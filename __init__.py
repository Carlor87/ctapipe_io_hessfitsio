'''
Event source for HESS data in fits format.
The data are passed with a single file for each telescope
Each telescope files in the name has a tag _CTX where X
is the telescope number

Author: Carlo Romoli - MPIK
Using the hessioeventsource module as base
'''

from astropy import units as u
from astropy.time import Time,TimeDelta
from fitsio import FITS      ## change all the code to use this library
from astropy.table import Table
from ctapipe.io.eventsource import EventSource
from ctapipe.io.containers import DataContainer
from ctapipe.instrument import CameraGeometry,OpticsDescription,TelescopeDescription, SubarrayDescription


__all__ = ['HESSfitsIOEventSource']


class HESSfitsIOEventSource(EventSource):
    """
    EventSource for the HESS data in fits file format.

    This class utilises astropy.io.fits to read the HESS fits files, and stores the
    information into the event containers.
    """
    _count = 0

    def __init__(self, config=None, tool=None, **kwargs):
        #super().__init__(config=config, tool=tool, **kwargs)
        super().__init__(config=config,**kwargs)
        try:
            import numpy as np
            import os as os
        except ImportError:
            msg = "The `numpy` python module is required to run the pipe"
            self.log.error(msg)
            raise
        self.np = np
        self.os = os


    @staticmethod
    def is_compatible(file_path):
        '''
        Check if the fits file is a valid one and
        if it has the right columns in it
        use build-in commands in astropy.io.fits
        '''
        listfiles = [f for f in file_path if
                     (f[-8:-6] == 'CT' and int(f[-6]))]
        for f in listfiles:
            with FITS(f) as hdul:
                try:
                    try:
                        hdul['TEVENTS'].has_data() #check if the table has data inside
                    except IOError:
                        msg = "Table TEVENTS not found!"
                        self.log.error(msg)
                except FileNotFoundError:
                    msg = "File not found, please provide existing file"
                    self.log.error(msg)
        return True
                

    def __exit__(self, exc_type, exc_val, exc_tb):
        '''
        !!!TO BE DEFINED!!!
        '''
        #hdul.close()
        #del table
        

    def _generator(self):
        # avoid extra dependences and check on CTX tag and with the list of allowed telescopes
        # WARNING, requires a rigid structure for the end of the file name: "*CTX.fits"
        # TO BE modified...it will break with more than 9 telescopes :(
        # option to import the regex package to make the selection more powerful
        listfiles = [self.input_url+f for f in self.os.listdir(self.input_url) if (f[-8:-6]=='CT' and int(f[-6]) in self.allowed_tels)]
        listfitstable = [] # store the table with the events for each file
        listpointingtable = [] # store the pointing correction table
        ntelescope = 0 # counter for the telescopes (maybe not needed)
        ev = 0 # event counter
        if len(listfiles)==0:
            print("File list empty, check the input directory or the list of allowed telescopes")
        for f in listfiles:
            hdul = FITS(f)
            # the container is initialized once, and data is replaced within
            # it after each yield
            eventtable = hdul['TEVENTS']
            listfitstable.append(eventtable)
        maxevent = max(table.read_header()['NAXIS2'] for table in listfitstable)
        data = DataContainer() # initialization of the container
        data.meta['origin'] = "HESSfits"
        data.meta['input_url'] = self.input_url
        data.meta['pointing_alt'] = listfitstable[0].read_header()['ALT_PNT']
        data.meta['pointing_az'] = listfitstable[0].read_header()['AZ_PNT']
        obs_id = listfitstable[0].read_header()['OBS_ID']  # same for all the tables, take only the first
        reftime = Time(listfitstable[0].read_header()['MJDREFI'] + listfitstable[0].read_header()['MJDREFF'],
                       format='mjd')
        data.dl0.obs_id = obs_id
        ct = self.np.zeros(len(listfiles)) # position of the event in each of the files
        while ev < maxevent:
            if ev == 0:
                data.inst.subarray = self._build_subarray_info(f) #based on the last file. The table is always the same
            data.count = ev
            evdata = []
            j = 0
            for table in listfitstable:
                if table['EVENT_ID_HESS'][ct[j]] == (ev + 1):
                    evdata.append(table[ct[j]])
                    ct[j] = ct[j]+1
                j = j + 1
            event_id = [evdata[x]['EVENT_ID_HESS'] for x in range(len(evdata))] # take the last element of the array
            data.pointing.azimuth = evdata[-1]['AZ_PNT']  # azimuth per event - needs to be a number
            data.pointing.altitude = evdata[-1]['ALT_PNT']  # altitude per event - needs to be a number
            tels_with_data = [evdata[x]['TEL_ID'] for x in range(len(evdata))]
            data.trig.gps_time = (reftime + TimeDelta(evdata[-1]['TIME'], format='sec'))
            data.trig.tels_with_trigger = tels_with_data
            data.count = ev
            data.r0.obs_id = obs_id
            data.r0.event_id = event_id
            data.r0.tels_with_data = tels_with_data
            data.r1.obs_id = obs_id
            data.r1.event_id = event_id
            data.r1.tels_with_data = tels_with_data
            data.dl0.obs_id = obs_id
            data.dl0.event_id = event_id
            data.dl0.tels_with_data = tels_with_data
            
            data.r0.tel.clear()
            data.r1.tel.clear()
            data.dl0.tel.clear()
            data.dl1.tel.clear()
            count = 0
            for tel_id in tels_with_data:
                npix = len(data.inst.subarray.tel[tel_id].camera.pix_id)
                data.dl1.tel[tel_id].image = self.np.zeros(npix) # initialize with empty image
                # write the image in the container
                data.dl1.tel[tel_id].image[evdata[count]['TEL_IMG_IPIX']] = evdata[count]['TEL_IMG_INT']
                count += 1
            ev += 1
            yield data
        return


    def _build_subarray_info(self,fitsfile):
        '''
        Added also the reading of the camera file "chercam.fits.gz"
        with the position of the pixels.
        The HESS-I camera in the ctapipe database does not corresespond to
        the real camera (it works only for the Montecarlo data).

        Parameters
        ----------
        fitsfile: fitsfile with the table of about the array details

        Returns
        -------
        SubarrayDescription :
            instrumental information
        '''
        pathtofile = self.os.path.split(self.input_url)[0]+'/'
        chercamfile = pathtofile+"chercam.fits.gz"
        pixtab = Table.read(chercamfile,format='fits') #read the table with the pixel position
        subarray = SubarrayDescription("HESS-I")
        try:
            hdu_array = FITS(fitsfile)[3] # open directly the table with the telarray
            teldata = hdu_array.read()
            telescope_ids = list(teldata['TELID'])

            for tel_id in telescope_ids:
                cam=pixtab[(tel_id-1)*960:tel_id*960]
                geom = CameraGeometry('CT%i'%tel_id,
                                      cam['PIX_ID'],
                                      self.np.array(cam['PIX_POSX'])*u.m,
                                      self.np.array(cam['PIX_POSY'])*u.m,
                                      self.np.array(cam['PIX_AREA'])*u.m*u.m,
                                      pix_type='hexagonal')
                foclen = teldata['FOCLEN'][tel_id-1] * u.m
                mirror_area = 108*u.m*u.m  # hard coded, NEED FIX! This column is empty in the original fits file!
                num_tiles = 382            # hard coded, missing in the original file
                optic = OpticsDescription(name = 'MST',
                                          num_mirrors = 1,
                                          equivalent_focal_length = foclen,
                                          mirror_area = mirror_area,
                                          num_mirror_tiles = num_tiles)
                tel_pos = [
                    teldata['POSX'][tel_id-1],
                    teldata['POSY'][tel_id-1],
                    teldata['POSZ'][tel_id-1]
                    ] * u.m
                tel = TelescopeDescription(name = 'CT%i'%tel_id, tel_type = 'MST',optics = optic,camera = geom)
                subarray.tels[tel_id] = tel
                subarray.positions[tel_id] = tel_pos
            
            return subarray
        except FileNotFoundError:
            msg = "Secondary file not found, check the presence in same folder!"
            self.log.error(msg)
            raise SystemExit

