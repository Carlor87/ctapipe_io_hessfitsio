import numpy as np
import matplotlib.pyplot as plt

from tqdm import tqdm
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from ctapipe.image import hillas, tailcuts_clean
from ctapipe.reco import HillasReconstructor
from ctapipe.reco.reco_algorithms import InvalidWidthException

from ctapipe_io_hessfitsio import HESSfitsIOEventSource as hfio
from ctapipe_io_hessfitsio.pointingcorrections import PointingCorrections as pc

if __name__ == "__main__":
    inputfile = hfio(input_url="/home/cromoli/Documents/hessfits/singlefiles_head/",allowed_tels=[1, 2, 3, 4])
    loc = EarthLocation.from_geodetic(lon = 16.5002222222222*u.deg,lat=-23.2717777777778*u.deg,height=1835*u.m) # HESS site location

    ShowerList = []
    ImageList = []
    MaskList = []
    HillasList = []
    CoordList = []
    
    # Initialize the pointing corrections at the beginning so to load the tables only once.
    pointcorr = pc.PointingCorrections("/home/cromoli/Documents/hessfits/singlefiles_head/")
    
    with inputfile as source:
        gsource = (x for x in source)
        hillasReco = HillasReconstructor()
        for i in tqdm(range(1000)):
            event = next(gsource)
            HillasIm = dict()
            NomHillasIm = dict()
            ExtendedImage = dict()  # only for presentation purposes
            mask0510 = dict()
            for j in event.dl0.tels_with_data:
                foc = event.inst.subarray.tel[j].optics.equivalent_focal_length
                ExtendedImage[j] = event.dl1.tel[j].image
                mask0510[j] = tailcuts_clean(event.inst.subarray.tel[j].camera,
                                             event.dl1.tel[j].image,
                                             10,5)  # implement a (5,10) cleaning
                try:
                    HillasIm[j] = hillas.hillas_parameters(event.inst.subarray.tel[j].camera,
                                                           event.dl1.tel[j].image * mask0510[j])
                except hillas.HillasParameterizationError:
                    continue
            ImageList.append(ExtendedImage)  # only for presentation purposes
            MaskList.append(mask0510)
            HillasList.append(HillasIm)
            if len(HillasIm) < 2:
                # at least 2 telescopes
                continue
            altaz = AltAz(location=loc, obstime=Time("2006-07-28 00:34:46")) # for testing, initial time of observation of the run
            array_pointing = SkyCoord(alt = 82.918*u.deg, az=191.869*u.deg, frame = altaz) # for testing, initial pointing position of the run
            obstime = event.trig.gps_time
            #altaz = AltAz(location=loc, obstime=obstime)
            #array_pointing = SkyCoord(alt = event.pointing.altitude*u.deg, az=event.pointing.azimuth*u.deg, frame = altaz)
            corrections = np.array([pointcorr.get_delta_azalt(t,obstime) for t in event.inst.subarray.tel_ids])
            #telescope_pointing = SkyCoord(alt=(event.pointing.altitude+corrections[:,1])*u.deg,
            #                              az=(event.pointing.azimuth+corrections[:,0])*u.deg,
            #                              frame=altaz)
            telescope_pointing = SkyCoord(alt=(array_pointing.alt+corrections[:,1]*u.deg),
                                          az=(array_pointing.az+corrections[:,0]*u.deg),
                                          frame=altaz) # for testing, corrections on initial pointing position of the run
            telescope_pointing = {x:y for x,y in zip(event.inst.subarray.tel_ids,telescope_pointing)} # create dict
            try:
                shower = hillasReco.predict(hillas_dict=HillasIm,
                                            inst=event.inst,
                                            array_pointing = array_pointing,
                                            telescopes_pointings = telescope_pointing
                                            )
            except InvalidWidthException:
                continue
            CoordList.append(SkyCoord(shower['az'],shower['alt'],frame=altaz))
            ShowerList.append(shower)
        print("End process!")


    # storing the main results of the shower reconstruction into lists
    # for better analysis and checks
    corx = [x['core_x'].value for x in ShowerList]
    cory = [x['core_y'].value for x in ShowerList]
    alt = [x['alt'].value for x in ShowerList]
    az = [x['az'].value for x in ShowerList]
    hmax = [x['h_max'].value for x in ShowerList]
    ra = [x.fk5.ra.value for x in CoordList]
    dec = [x.fk5.dec.value for x in CoordList]
    
    
    # produces a rough skymap according to the initial
    # pointing of the array and compare with the position
    # of PKS2155-304 and of the array pointing.
    plt.figure()
    plt.hist2d(ra,dec,bins=[np.linspace(324,334,50),np.linspace(-35,-25,50)])
    plt.xlabel("Ra")
    plt.ylabel("Dec")
    plt.plot( 329.716938 ,  -30.225589 ,'ro')
    plt.plot(array_pointing.fk5.ra.value,array_pointing.fk5.dec.value,'go')
    plt.show()
    
    
    
    