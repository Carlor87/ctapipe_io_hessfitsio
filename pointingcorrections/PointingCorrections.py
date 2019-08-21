"""
Class to read and implement the
pointing corrections for HESS
analysis of data

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
        self.dx = 0
        self.dy = 0
        self.phi = 0
        self.focal = 0
        self.delta_az = 0  # horizontal direction
        self.delta_alt = 0 # vertical direction
        self.table = 0     # table to store the entire fits table
        self._readpc(url)  # read the pointing corrections
        
        
    def _interpolate(self,tel,evtime):
        self.dx  = np.interp(evtime,self.table[str(tel)]['TIME'][:],self.table[str(tel)]['TRANSLATION_X'][:])
        self.dy  = np.interp(evtime,self.table[str(tel)]['TIME'][:],self.table[str(tel)]['TRANSLATION_Y'][:])
        self.phi = np.interp(evtime,self.table[str(tel)]['TIME'][:],self.table[str(tel)]['ROTATION_ANGLE'][:])
        
        
    def _readpc(self,url):
        listfiles = [url+f for f in os.listdir(url) if (f[-8:-6]=='CT')][::-1] # reverse the order of the lists
        tels = [t[-6:-5] for t in listfiles]
        self.table = dict()
        for i in range(len(listfiles)):
            pointingcorr = FITS(listfiles[i])[2]
            self.table[tels[i]] = pointingcorr
        self.focal = pointingcorr['FOCAL_DISTANCE'][0] #focal length does not change
        header = FITS(listfiles[0])[1].read_header()
        self._reftime = Time(header['MJDREFI'] + header['MJDREFF'], format='mjd')
        
        
    def _from_camera_to_nominal(self):
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
        Return the pointing corrections for all the telescope
        """
        eventtime = (evtime - self._reftime).sec
        self._interpolate(tel,eventtime)
        self._from_camera_to_nominal()
        return self.delta_az,self.delta_alt
        
    
        
        
    
    
        
