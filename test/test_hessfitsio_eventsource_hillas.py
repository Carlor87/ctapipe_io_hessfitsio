from tqdm import tqdm
from ctapipe.image import hillas, tailcuts_clean

from ctapipe_io_hessfitsio import HESSfitsIOEventSource as hfio


if __name__ == "__main__":
    inputfile = hfio(input_url="/home/cromoli/Documents/hessfits/singlefiles_head/",allowed_tels=[1, 2, 3, 4])
    
    ImageList = []
    MaskList = []
    HillasList = []
    with inputfile as source:
        gsource = (x for x in source)
        for i in tqdm(range(1000)):
            event = next(gsource)
            HillasIm = dict()
            ExtendedImage = dict()  # only for presentation purposes
            mask0510 = dict()
            for j in event.dl0.tels_with_data:
                ExtendedImage[j] = event.dl1.tel[j].image
                mask0510[j] = tailcuts_clean(event.inst.subarray.tel[j].camera,
                                             event.dl1.tel[j].image,
                                             10,5)  # implement a (5,10) cleaning
                try:
                    HillasIm[j] = hillas.hillas_parameters(event.inst.subarray.tel[j].camera,
                                                           event.dl1.tel[j].image * mask0510[j])
                except hillas.HillasParameterizationError:
                    continue
    print("End process!")

