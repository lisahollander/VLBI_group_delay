from dataclasses import dataclass
from numpy import exp, cos, sin, pi, sqrt, degrees
import numpy as np

@dataclass(order=True)
class Gaussian:
    amp: float = 0
    x0: float = 0           
    y0: float = 0           
    a: float = 0
    b: float = 0
    theta: float = 0
    a_scaled: float = 0
    b_scaled: float = 0
    delta_RA: float = 0
    delta_DEC: float = 0

    def __add__(self, gaussian):
        return Gaussian(gaussian.amp + self.amp, gaussian.x0 + self.x0, gaussian.y0 + self.y0, gaussian.a + self.a,
                        gaussian.b + self.b, gaussian.theta + self.theta)

    def __mul__(self, number):
        return Gaussian(number * self.amp, number * self.x0, number * self.y0, number * self.a,
                        number * self.b, number * self.theta)

    def __iter__(self):
        return self

    def __str__(self):
        return f"{round(self.amp, 5):<15}{round(self.x0, 2):<15}"\
                               f"{round(self.y0, 2):<15}"\
                               f"{round(sqrt(1 / self.a), 2):<15}{round(sqrt(1 / self.b), 2):<15}"\
                               f"{round(float(degrees(self.theta % pi)), 2)}\n"

    def __format__(self, format_spec):
        return f"{round(self.amp, 5):<15}{round(self.x0, 2):<15}" \
               f"{round(self.y0, 2):<15}" \
               f"{round(sqrt(1 / self.a), 2):<15}{round(sqrt(1 / self.b), 2):<15}" \
               f"{round(float(degrees(self.theta % pi)), 2)}\n"

    def get_gauss(self, y, x):
        """
        Returns value of gaussian for coordinates x and y
        """
        cos_2_theta, sin_2_theta = cos(2 * self.theta), sin(2 * self.theta)
        x -= self.x0
        y -= self.y0

        x2, y2, xy = x * x, y * y, x * y

        a_term = (x2 * (1.0 + cos_2_theta) / 2 + y2 * (1.0 - cos_2_theta) / 2 + xy * sin_2_theta)
        b_term = (x2 * (1.0 - cos_2_theta) / 2 + y2 * (1.0 + cos_2_theta) / 2 - xy * sin_2_theta)
        g = abs(self.a * a_term + self.b * b_term)
        g[g>50] = 50    #Any value above 50 is set to 50

        return exp(-g) * self.amp

    def get_terms(self, y, x):
        """
        Returns terms to build up partials for normal vector to the least square method in a specific coordinate
        """
        cos_2_theta, sin_2_theta = cos(2 * self.theta), sin(2 * self.theta)
        x -= self.x0
        y -= self.y0

        x2, y2, xy = x * x, y * y, x * y

        return [(x2 * (1.0 + cos_2_theta) / 2 + y2 * (1.0 - cos_2_theta) / 2 + xy * sin_2_theta),
                (x2 * (1.0 - cos_2_theta) / 2 + y2 * (1.0 + cos_2_theta) / 2 - xy * sin_2_theta),
                x + x * cos_2_theta + y * sin_2_theta,
                x - x * cos_2_theta - y * sin_2_theta,
                y - y * cos_2_theta + x * sin_2_theta,
                y + y * cos_2_theta - x * sin_2_theta,
                x2 * sin_2_theta - y2 * sin_2_theta - 2 * xy * cos_2_theta
                ]
    
    def get_fourier_transform_value(self, y,x, scale_factor=1):
        x_bar = (x - scale_factor / 2) / scale_factor
        y_bar = (y - scale_factor / 2) / scale_factor 

        return self.amp * pi / sqrt(self.a * self.b) * exp(-2 * pi * 1j * (self.x0 * x_bar + self.y0 * y_bar)) * \
            exp(-pi**2 * ((x_bar * cos(self.theta) + y_bar * sin(self.theta)) ** 2 / self.a + (y_bar * cos(self.theta) - x_bar * sin(self.theta)) ** 2 / self.b))
    
    
    def get_gaussian_value_2(self,header):
        """
        Returns value of gaussian for coordinates according to the fits files in ascention and declination
        """
        cos_2_theta, sin_2_theta = cos(2 * (self.theta)), sin(2 * self.theta)

        self.delta_RA = (self.x0 - header.reference_pixel_RA) * header.pixel_increment_RA          
        self.delta_DEC = (self.y0 - header.reference_pixel_DEC) * header.pixel_increment_DEC       

        x2, y2, xy = (header.grid_RA - self.delta_RA)**2, (header.grid_DEC - self.delta_DEC)**2, (header.grid_DEC - self.delta_DEC) * (header.grid_RA - self.delta_RA)

        self.a_scaled = self.a / header.pixel_increment_RA**2            # Note: We assume pixel increment is the same for a and b, pixel_increment_RA = pixel_increment_DEC
        self.b_scaled = self.b / header.pixel_increment_DEC**2           # Note: We assume pixel increment is the same for a and b, pixel_increment_RA = pixel_increment_DEC

        a_term = x2 * (1.0 + cos_2_theta) / 2 + y2 * (1.0 - cos_2_theta) / 2 + xy * sin_2_theta 
        b_term = x2 * (1.0 - cos_2_theta) / 2 + y2 * (1.0 + cos_2_theta) / 2 - xy * sin_2_theta

        g = abs(self.a_scaled * a_term + self.b_scaled * b_term)
        g[g > 50] = 50   #Any value above 50 is set to 50

        return self.amp * exp(-g)

    def get_fourier_transform_value_2(self,U,V,header):
        c = 3*10**8
        w = header.w

        A = self.amp * pi / sqrt(self.a_scaled * self.b_scaled) 

        dTheta= 1/c * (self.delta_RA * U + self.delta_DEC * V)
        dExp = w/c**2 * ((U * cos(self.theta) + V * sin(self.theta))**2/(4*self.a_scaled) + (V * cos(self.theta) - U * sin(self.theta))** 2/(4*self.b_scaled)) 

        visibility = A * exp(-1j* w * dTheta) * exp(-w * dExp)
        dV = visibility * (-1j * dTheta - 2*dExp)

        return visibility, dV

@dataclass
class GaussList:
    """
    Class for storing several gaussians to be depicted. Define size for the class to define the size of the image
    of the gaussians that the build_image function returns.
    """

    def __init__(self, gaussians=None, size=256):
        self.size = size
        self.gaussians = list(gaussians) if gaussians else []

    def __getitem__(self, num):
        return self.gaussians[num]

    def __str__(self):
        return f"{[str(gauss) for gauss in self.gaussians]}"

    def __len__(self):
        return len(self.gaussians)

    def append(self, gaussian):
        if isinstance(gaussian, Gaussian):
            self.gaussians.append(gaussian)
        elif isinstance(gaussian, GaussList):
            self.gaussians.extend(gaussian.gaussians)
        else:
            raise ValueError(f"tried to append {type(gaussian)} to GaussList")
        
    def build_image_2(self,header):
        result = np.zeros((header.size,header.size))
        for gauss in self.gaussians:
            result = result + gauss.get_gaussian_value_2(header)

        return result
    
    def build_image(self):
        return sum([np.fromfunction(gauss.get_gauss, (self.size, self.size), dtype=float) for gauss in self.gaussians])
    
    def get_analytical_results(self,U,V,gauss_fnd,header):
        anl2 = np.zeros_like(U,dtype=np.complex128)
        anl2_dV = np.zeros_like(U,dtype=np.complex128)

        for gauss in gauss_fnd:
            visibility, dV = gauss.get_fourier_transform_value_2(U,V,header)
            anl2 += visibility
            anl2_dV += dV
        
        anl2_derivative = (np.real(anl2) * np.imag(anl2_dV) - np.imag(anl2) * np.real(anl2_dV)) / (np.real(anl2)**2 + np.imag(anl2)**2)

        return anl2, anl2_derivative

    def add(self, delta, step):
        return GaussList([gauss + d_gauss * step for gauss, d_gauss in zip(self.gaussians, delta)], size=self.size)

    def sort(self):
        """
        Sorts GaussList in ascending order because that's the order that peak_local_max function returns the peaks.
        Just so the order of final components and original components match when modelling analytical gaussians.
        """
        self.gaussians.sort(key=lambda x: x.amp, reverse=True)
