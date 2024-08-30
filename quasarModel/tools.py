from datetime import datetime
import time
import gzip
import os
import logging
import logging.handlers

import csv
from logging import getLogger
from source_model import SourceModel
from read_file import FileData
from numpy import sqrt

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
Tools for running and logging several files
##################################################
"""

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
                file_data = FileData()
                image = file_data.get_image_from_path(source_path) 

                org, mdl, _, _, _, _,_,gauss_fnd = source_model.process(image,file_data)
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
