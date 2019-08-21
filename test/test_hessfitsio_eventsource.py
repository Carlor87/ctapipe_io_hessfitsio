"""
Testing script for ctapipe_io_hessfitsio module

Runs through the files of a single run
reading the first 100 events

Only for testing
No reconstruction involved, only reading and
filling the DataContainer with dl1 infos for
the events
"""
from tqdm import tqdm

from ctapipe_io_hessfitsio import HESSfitsIOEventSource as hfio


if __name__ == "__main__":
    inputfile = hfio(input_url="/home/cromoli/Documents/hessfits/singlefiles_head/",allowed_tels=[1, 2, 3, 4])

    with inputfile as source:
        gsource = (x for x in source)
        for i in tqdm(range(1000)):
            event = next(gsource)
            
    print("End process!")

