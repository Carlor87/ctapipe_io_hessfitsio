"""
Class to read and implement the
pointing corrections for HESS
analysis of data.
WARNING: It works only for the HESS-I telescopes

author: Carlo Romoli - MPIK

date: 2019-08-20

The file is based on the original transormations
present in the HAP software.
"""
import os
import numpy as np
from fitsio import FITS
from astropy.time import Time,TimeDelta

__all__ = ["PointingCorrections"]

class PointingCorrections:
    def __init__(self,url):
        """
        Initialization of the class.
        Initialize to 0 all the needed quantities and reads
        calls the function to read the fits table where the
        pointing corrections are stored.
        Arguments:
           - url:
               path to the folder with the fits files of the run
        """
        self.dx = 0
        self.dy = 0
        self.phi = 0
        self.focal = 0
        self.delta_az = 0  # horizontal direction
        self.delta_alt = 0 # vertical direction
        self.table = 0     # table to store the entire fits table
        self._readpc(url)  # read the pointing corrections
        
        
    def _interpolate(self,tel,evtime):
        """
        Interpolate the pointing correction table and fill the corresponding variable.
        Given the time of the event, linear interpolation between the values stored in the
        FITS table (which are every ~second).
        Arguments:
            - tel:
                telescope number (as int then converted in string to access the dictionary)
            - evtime:
                time of the event (in MET)
        """
        self.dx  = np.interp(evtime,self.table[str(tel)]['TIME'][:],self.table[str(tel)]['TRANSLATION_X'][:])
        self.dy  = np.interp(evtime,self.table[str(tel)]['TIME'][:],self.table[str(tel)]['TRANSLATION_Y'][:])
        self.phi = np.interp(evtime,self.table[str(tel)]['TIME'][:],self.table[str(tel)]['ROTATION_ANGLE'][:])
        
        
    def _readpc(self,url):
        """
        Called in the initialization process of the class
        Reads the FITS table with the pointing corrections and stores it in a
        dictionary of tables usinf the telescope number as key.
        Sets up the focal distance of the telescope and reads the reference time
        of the MET data format.
        Arguments:
           - url:
               path to the folder with the fits files of the run
        """
        listfiles = [url+f for f in os.listdir(url) if (f[-8:-6]=='CT')][::-1] # reverse the order of the lists
        tels = [t[-6:-5] for t in listfiles]
        self.table = dict()
        for i in range(len(listfiles)):
            pointingcorr = FITS(listfiles[i])[2]
            self.table[tels[i]] = pointingcorr
        self.focal = pointingcorr['FOCAL_DISTANCE'][0] #focal length does not change
        header = FITS(listfiles[0])[1].read_header()   # read header of the first files (has same information for the reftime)
        self._reftime = Time(header['MJDREFI'] + header['MJDREFF'], format='mjd') # sets up the reference time for the MET
        
        
    def _from_camera_to_nominal(self):
        """
        Converts the pointing corrections, which are in the camera frame,
        in an angular shift in azimuth and altitude.
        WARNING: works only for the HESS-I telescopes.
        It involves a rotation by the angle -phi, a scaling by the
        focal length (to obtain an angluar distance) of the translations
        in X and Y direction.
        The result is converted in a shift in units of degrees for
        the azimuth and altitude direction.
        """
        degtorad = np.pi/180.
        c = np.cos(-self.phi*degtorad)
        s = np.sin(-self.phi*degtorad)
        rot = np.array([[c,-s],[s,c]])
        delta = 1./self.focal*np.array([-self.dx,-self.dy])
        result = np.dot(rot,delta)
        self.delta_az = result[0]/degtorad
        self.delta_alt = result[1]/degtorad
    
    
    def get_delta_azalt(self,tel,evtime):
        """
        Return the pointing corrections for all the telescope.
        First construct the time of the event in MET
        (as a difference in seconds from self._reftime).
        Then interpolates the value to obtain the pointing corrections
        at the said time step. Computes the shift in azimuth and altitude.
        Arguments:
            - tel:
                Telescope number (as int)
            - evtime:
                Time of the event (as astropy.time.Time object)
            Returns:
                Tuple with shift in azimuth and altitude in units
                of degrees.
        """
        eventtime = (evtime - self._reftime).sec
        self._interpolate(tel,eventtime)
        self._from_camera_to_nominal()
        return self.delta_az,self.delta_alt
        
    
        
        
    
    
        
