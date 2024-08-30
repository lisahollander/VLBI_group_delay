
from dataclasses import dataclass, field

from numpy import meshgrid,arange,ndarray,pi,array
from astropy.io import fits
from constants import one_source1, one_source2, one_source3, one_source4, one_source5, one_source6, \
    two_sources


"""
##############################################################
Data class for reading fits files for image data acquisition 
and initialization of coordinate system
##############################################################
"""

@dataclass
class FileData:
    grid_RA: ndarray = field(default_factory=lambda: array([]))
    grid_DEC: ndarray = field(default_factory=lambda: array([]))
    u : ndarray = field(default_factory=lambda: array([]))
    v : ndarray = field(default_factory=lambda: array([]))
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


    def get_image_from_path(self,source_path):
        with fits.open(source_path) as data:
            image = data[0].data
            image.shape = image.shape[2:]

            header_file = data[0].header
            
            self.reference_pixel_RA = header_file['CRPIX1']
            self.reference_pixel_DEC = header_file['CRPIX2']
            self.pixel_increment_RA = abs(header_file['CDELT1'])
            self.pixel_increment_DEC = abs(header_file['CDELT2'])

            self.frequency = header_file['CRVAL3']
            self.size = header_file['NAXIS1']
            self.w = 2 * pi * self.frequency 
            self.wavelength = 3*10**8 / self.frequency

            self.pixel_increment_RA = self.pixel_increment_RA * pi/180
            self.pixel_increment_DEC = self.pixel_increment_DEC * pi/180

            self.initialize_coord_system()

        return image

def check_if_test_source(source_path):
    test_sources = {'1': one_source1, '2': one_source2, '3': one_source3, '4': one_source4, '5': one_source5,
                    '6': one_source6, '7': two_sources}
    if source_path in test_sources:
        source_path = test_sources[source_path]
    return source_path


