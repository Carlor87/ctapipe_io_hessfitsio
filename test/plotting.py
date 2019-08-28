from ctapipe.instrument import CameraGeometry
from ctapipe.visualization import CameraDisplay, ArrayDisplay
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.table import Table
import numpy as np
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

from ctapipe.coordinates import (
    GroundFrame,
    TiltedGroundFrame,
    NominalFrame,
    TelescopeFrame,
    EngineeringCameraFrame,
    CameraFrame,
)


camdef=Table.read("/home/cromoli/Documents/hessfits/pks2155flarerun_wImages_test.hap/chercam.fits.gz",format='fits')
# separate the different cameras in different tables
camdef1=camdef[0:960]
camdef2=camdef[960:1920]
camdef3=camdef[1920:2880]
camdef4=camdef[2880:3840]
geom = CameraGeometry(
    'HESSI',
    camdef3['PIX_ID'],
    np.array(camdef3['PIX_POSX']) * u.m,
    np.array(camdef3['PIX_POSY']) * u.m,
    np.array(camdef3['PIX_AREA']) * u.m * u.m,
    pix_type='hexagonal')  # Real camera

def PlotEventExtended(extendedimage):
    """
    Plot the extended images of the events in the 
    camera frame.
    """
    camdef=Table.read("/home/cromoli/Documents/hessfits/pks2155flarerun_wImages_test.hap/chercam.fits.gz",format='fits')
    camdef1=camdef[0:960]
    camdef2=camdef[960:1920]
    camdef3=camdef[1920:2880]
    camdef4=camdef[2880:3840]
    geom = CameraGeometry(
    'HESSI',
    camdef3['PIX_ID'],
    np.array(camdef3['PIX_POSX']) * u.m,
    np.array(camdef3['PIX_POSY']) * u.m,
    np.array(camdef3['PIX_AREA']) * u.m * u.m,
    pix_type='hexagonal')  # Real camera
    #geom2 = CameraGeometry.from_name('HESS-I') # only for MC
    ntrig = len(extendedimage)
    fig,axs = plt.subplots(1, ntrig, constrained_layout=True, figsize=(ntrig*3, 3))
    for ax,ima in zip(axs,extendedimage):
        CameraDisplay(geom, ax=ax, image=extendedimage[ima])
        ax.set_title("CT%i"%ima)
    

def PlotSummedImages(event,mask,hillas):
    plt.figure()
    disp = CameraDisplay(geom,title="HESS I")
    disp.cmap = 'jet'
    image=0
    for ii in event.dl0.tels_with_data:
        image += event.dl1.tel[ii].image*mask[ii]
        disp.image = image
        disp.overlay_moments(hillas[ii], color='pink', lw=3, with_label=False,keep_old=True)
        

def PlotArrayDisplay(event,hillas,shower,azim):
    plt.figure()
    disp = ArrayDisplay(event.inst.subarray)
    time_gradients = {}
    for telescope_id, dl1 in event.dl1.tel.items():
        time_gradients[telescope_id] = hillas[telescope_id].skewness
    disp.set_vector_hillas(hillas,
                           time_gradient = time_gradients,
                           angle_offset = azim, # not clear the meaning of this parameter
                           length=200)
    #disp.set_line_hillas(hillas,200) # this function will not take into account the azimuth angle of the array
    if (not np.isnan(shower.core_uncert)):
        plt.errorbar(shower.core_x.value,
                 shower.core_y.value,
                 xerr = shower.core_uncert.value,
                 yerr = shower.core_uncert.value,
                 fmt = 'ro',
                 label = 'Reconstructed Core')
    else:
        plt.plot(shower.core_x,
                 shower.core_y,
                 'ro',
                 label = 'Reconstructed Core')
    plt.xlabel("Distance [m]")
    plt.ylabel("Distance [m]")
    
    
    
def CheckFrames(event):
    loc = EarthLocation.from_geodetic(lon = 16.5000000*u.deg,lat=-23.2716667*u.deg,height=1800*u.m)
    focal_length = event.inst.subarray.tel[1].optics.equivalent_focal_length
    obstime = event.trig.gps_time
    altaz = AltAz(location=loc, obstime=obstime)
    telescope_pointing = SkyCoord(alt=event.pointing.altitude * u.deg, az=event.pointing.azimuth * u.rad, frame=altaz)
    camera = CameraFrame(focal_length=focal_length,rotation=0.*u.deg,telescope_pointing=telescope_pointing,obstime = event.trig.gps_time,location = loc)
    camera2 = CameraFrame(focal_length=focal_length,rotation=0.25*u.deg,telescope_pointing=telescope_pointing,obstime = event.trig.gps_time,location = loc)
    cam_coord = SkyCoord(geom.pix_x,geom.pix_y,frame = camera)
    cam_coord2 = SkyCoord((geom.pix_x-0.025*u.m),(geom.pix_y-0.025*u.m),frame = camera)
    pks = SkyCoord.from_name('pks 2155-304')
    pks_cam = pks.transform_to(camera)
    plt.figure()
    plt.scatter(cam_coord.x, cam_coord.y)
    plt.scatter(cam_coord2.x, cam_coord2.y) 
    plt.xlabel(f'x / {cam_coord.x.unit}') 
    plt.ylabel(f'y / {cam_coord.y.unit}') 
    plt.axis('square') 
    plt.plot(pks_cam.x.to_value(u.m),pks_cam.y.to_value(u.m),'ro')
    
    telescope_frame = TelescopeFrame(telescope_pointing = telescope_pointing)
    telescope_coords = cam_coord.transform_to(telescope_frame)
    wrap_angle = telescope_pointing.az + 180* u.deg
    plt.figure()
    plt.axis('equal')
    plt.scatter(
        telescope_coords.altaz.az.wrap_at(wrap_angle).deg,
        telescope_coords.altaz.alt.deg
    )
    plt.xlabel('x / {}'.format(telescope_coords.altaz.az.unit))
    plt.ylabel('y / {}'.format(telescope_coords.altaz.alt.unit))
    
    