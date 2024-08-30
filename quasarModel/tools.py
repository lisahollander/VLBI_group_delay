from datetime import datetime
import time
import gzip
import os
import logging
import logging.handlers
import math

from gauss import Gaussian, GaussList
from numpy import radians, meshgrid,linspace,arange, empty,ndarray,pi,sqrt
from astropy.io import fits
from astropy.table import Table
from constants import one_source1, one_source2, one_source3, one_source4, one_source5, one_source6, \
    two_sources

import csv
from logging import getLogger
from source_model import SourceModel
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
import numpy as np



"""
##################################################
Tools for logging 
##################################################
"""


# Custom filter use to format records
class ContextFilter(logging.Filter):
    def filter(self, record):
        setattr(record, 'utc', datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S.%f')[:-3])
        return True


# Set default logger
def set_logger(log_path='', prefix='', console=False, size=1000000, debug=False):
    # Functions needed to provide name of new compress file
    def namer(name):
        folder = os.path.dirname(name)
        return os.path.join(folder, datetime.utcnow().strftime(f'{prefix}%Y-%m-%d.%H%M%S.gz'))

    # Functions needed to created file rotator with gzip compression
    def rotator(source, destination):
        with open(source, "rb") as sf, open(destination, "wb") as df:
            df.write(gzip.compress(sf.read(), 9))
        os.remove(source)

    logger = logging.getLogger('QuasarModelLog')
    logger.setLevel(logging.DEBUG)
    logger.addFilter(ContextFilter())
    formatter = logging.Formatter('%(message)s')
    # Add File handler
    if log_path:
        fh = logging.FileHandler(log_path, 'w')
        fh.setLevel(logging.INFO)
        fh.setFormatter(formatter)
        fh.rotator = rotator
        fh.namer = namer
        logger.addHandler(fh)
    # Add console filter
    if console:
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG if debug else logging.INFO)
        ch.setFormatter(formatter)
        logger.addHandler(ch)
    return logger

# decorator to calculate duration taken by any function.
def time_it(func):
    # added arguments inside the inner1,
    # if function takes any arguments,
    # can be added like this.
    def inner_fnc(*args, **kwargs):
        # storing time before function execution
        begin = time.time()

        results = func(*args, **kwargs)

        # storing time after function execution
        end = time.time()
        logger = logging.getLogger('QuasarModelLog')
        logger.info(f'{func.__name__} took {end - begin:.2f} seconds')
        return results

    return inner_fnc

"""
##################################################
Data class for image data acquisition 
##################################################
"""

@dataclass
class Header:
    grid_RA: ndarray = field(default_factory=lambda: np.array([]))
    grid_DEC: ndarray = field(default_factory=lambda: np.array([]))
    u : ndarray = field(default_factory=lambda: np.array([]))
    v : ndarray = field(default_factory=lambda: np.array([]))
    pixel_increment_RA: float = 0
    pixel_increment_DEC: float = 0
    reference_pixel_RA: float = 0
    reference_pixel_DEC: float = 0
    frequency: float = 0
    w : float = 0
    size: int = 0

    def initialize_coord_system(self):
        #Initializing x,y (RA,DEC) grid
        ascStart =  -(self.reference_pixel_RA * self.pixel_increment_RA)
        ascEnd =  (self.size - self.reference_pixel_RA) * self.pixel_increment_RA
                
        decStart = -(self.reference_pixel_DEC * self.pixel_increment_DEC)
        decEnd =  (self.size - self.reference_pixel_DEC) * self.pixel_increment_DEC

        right_ascension = arange(ascStart, ascEnd, self.pixel_increment_RA) 
        declination = arange(decStart, decEnd, self.pixel_increment_DEC)

        self.grid_RA, self.grid_DEC = meshgrid(right_ascension, declination)

        increment_size = 2*6371000/self.size

        uStart = -(self.reference_pixel_RA * increment_size)
        uEnd =  (self.size - self.reference_pixel_RA) * increment_size

        vStart = -(self.reference_pixel_RA * increment_size)
        vEnd =  (self.size - self.reference_pixel_RA) * increment_size

        u = arange(uStart,uEnd, increment_size) 
        v = arange(vStart,vEnd, increment_size)  
        self.u, self.v = meshgrid(u,v)


    def get_image_from_path_2(self,source_path):
        with fits.open(source_path) as data:
            image = data[0].data
            image.shape = image.shape[2:]

            header_file = data[0].header
            
            self.reference_pixel_RA = header_file['CRPIX1']
            self.reference_pixel_DEC = header_file['CRPIX2']
            self.pixel_increment_RA = abs(header_file['CDELT1'])
            self.pixel_increment_DEC = abs(header_file['CDELT2'])
            print('Pixel in degrees', self.pixel_increment_RA, self.pixel_increment_DEC)
            self.frequency = header_file['CRVAL3']
            self.size = header_file['NAXIS1']
            self.w = 2 * pi * self.frequency 
            self.wavelength = 3*10**8 / self.frequency

            clean_components = Table(data[1].data).to_pandas() 
            clean_components['X'] = (clean_components['DELTAX'] / self.pixel_increment_RA + self.reference_pixel_RA)
            clean_components['Y'] = (clean_components['DELTAY'] / self.pixel_increment_DEC + self.reference_pixel_DEC)

            for i in range(len(clean_components['X'])):
                if clean_components['X'][i] >= self.size:
                    clean_components['X'][i] = clean_components['X'][i] - 1
                    
                if clean_components['Y'][i] >= self.size:
                    clean_components['Y'][i] = clean_components['Y'][i] - 1

            self.pixel_increment_RA = self.pixel_increment_RA * pi/180
            self.pixel_increment_DEC = self.pixel_increment_DEC * pi/180

            self.initialize_coord_system()

        return image,clean_components

def check_if_test_source(source_path):
    test_sources = {'1': one_source1, '2': one_source2, '3': one_source3, '4': one_source4, '5': one_source5,
                    '6': one_source6, '7': two_sources}
    if source_path in test_sources:
        source_path = test_sources[source_path]
    return source_path


def get_image_from_path(source_path):
    image = fits.getdata(source_path, ext=0)
    image.shape = image.shape[2:]

    return image


"""
##################################################
Tools for running and logging several files
##################################################
"""

#TO DO: Fix this part of the code after implementation 
def model_several_sources(base_path):
    """
    Model several sources by providing the directory path to where they are located. Will create two .csv files.
    precision.csv stores a summary of the modelled precision of all sources.
    parameters.csv stores the parameters of the modelled gaussians
    :param base_path: path to directory with source images.
    directory structure should be: base_path / source_dir / source.fits
    :return: None
    """
    logger = getLogger('QuasarModelLog')

    header1 = ['Source', 'Mean error', 'Total error', 'mean error as % of max val', 'total error as % of max val']
    header2 = ['Source', 'Amplitude', 'X0', 'Y0', 'Sigma_x', 'Sigma_y', 'theta']

    source_model = SourceModel()
    with open('precisions.csv', 'w', encoding='UTF8') as f1, open('parameters.csv', 'w', encoding='UTF8') as f2:
        writer1 = csv.writer(f1)
        writer1.writerow(header1)

        writer2 = csv.writer(f2)
        writer2.writerow(header2)

        for index, source_path in enumerate([os.path.join(path, name)
                                             for path, _, files in os.walk(base_path) for name in files if name.endswith('map.fits')], 1):
            try:
                header = Header()
                image,_ = header.get_image_from_path_2(source_path) 
                radius_earth = 6378.1*10**3
                diameter_earth = 12756000
                u = np.linspace(-diameter_earth,diameter_earth, header.size)
                v = np.linspace(-diameter_earth,diameter_earth, header.size) 
                U, V = np.meshgrid(u,v)

                org, mdl, _, gauss_fnd = source_model.process(image,header,U,V)
                precision = source_model.precision(org, mdl)

                data1 = [source_path, precision['mean error'],
                         precision['total error'],
                         precision['mean error as % of max val'],
                         precision['total error as % of max val']]
                writer1.writerow(data1)

                writer2.writerow([source_path])
                for index, gauss in enumerate(gauss_fnd):
                    data2 = [f"Gauss #{index+1}", gauss.amp, gauss.x0, gauss.y0, 1/sqrt(gauss.a), 1/sqrt(gauss.b), gauss.theta]
                    writer2.writerow(data2)

                logger.info(f'{source_path} ok')
            except Exception as exc:
                logger.error(f'{source_path} {str(exc)}')
