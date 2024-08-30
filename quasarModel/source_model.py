import numpy as np
import matplotlib.pyplot as plt
from numpy import log, amax, exp, sqrt
from skimage.feature import peak_local_max
from gauss import Gaussian, GaussList
import copy
import logging
import scipy
from plot import *
from constants import stars, line


class SourceModel:
    stop: bool = False

    def process(self,image,file_data):
        """
        Model a quasar image with gaussian functions.
        :param image: np.array of image data
        :param file_data: object containig data from fits file header, RA/DEC grid and u/v grid
        :return: org - original input image, mdl - modelled image, anl - analytical fourier transform of modelled image
        """
        logger = logging.getLogger('QuasarModelLog')

        size = file_data.size

        org = copy.deepcopy(image)       #Creates copy of source image array

        gauss_fnd = GaussList(size=size) #Initializes list of gaussians (empty here)
        for guess in range(10):
            if self.stop:
                break

            #skimage.feature.peak_local_max built in function that extracts maximums in array
            peaks_position = peak_local_max(image, min_distance=3, threshold_abs=amax(org)*0.1,     
                                            exclude_border=True, num_peaks=10)

            num_peaks = len(peaks_position)
            if num_peaks == 0:
                break
                
            gauss_now = GaussList(size=size)
            for i in range(num_peaks):
                gauss_now.append(self.initial_guess(peaks_position[i][1], peaks_position[i][0], image))    
            
            mean_sq_dof = self.find_mean_sq(image, gauss_now)                                             
            mean_sq_dof_old = mean_sq_dof
        
            self.log_initial_guess(logger, num_peaks, mean_sq_dof, guess, gauss_now,file_data)              #Adds initial guess to logger                       

            gauss_now, mean_sq_dof, point_min = self.model_fitting(gauss_now, image)
            delta_mean = mean_sq_dof_old - mean_sq_dof
            mean_sq_dof_old = mean_sq_dof

            self.log_iteration(logger, mean_sq_dof, delta_mean, point_min, 1, gauss_now,file_data)

            for iteration in range(2, 12):
                if self.stop:
                    break

                if point_min == 0 or abs(delta_mean) < 10 ** (-8) or iteration == 11:
                    gauss_fnd.append(gauss_now)
                    image_comp = gauss_fnd.build_image()                    #Brightness distribution using image coordinate system
                    image = image - image_comp
                    break
                else:
                    gauss_now, mean_sq_dof, point_min = self.model_fitting(gauss_now, image)
                    delta_mean = mean_sq_dof_old - mean_sq_dof
                    mean_sq_dof_old = mean_sq_dof
                    self.log_iteration(logger, mean_sq_dof, delta_mean, point_min, iteration, gauss_now,file_data)

        mdl_fitted = gauss_fnd.build_image()
        precision = self.precision(org, mdl_fitted)

        mdl = gauss_fnd.build_image_scaled(file_data)                        #Brightness distribution using physical coordinate system and scaled Gaussians
        self.log_final_iteration(logger, gauss_fnd,file_data)

        anl, anlDerivative = gauss_fnd.get_analytical_results(gauss_fnd,file_data)     

        residuals, rms_value = self.calculate_group_delay_rms(file_data, anlDerivative)

        logger.info(f"mean error: {precision['mean error']}")
        logger.info(f"mean error as % of max val: {precision['mean error as % of max val']}")
        logger.info(f"total error: {precision['total error']}")
        logger.info(f"total error as % of max val: {precision['total error as % of max val']}")
        logger.info(f"RMS-error for group delay for baselines smaller than earth radius: {rms_value}")
  
        return org, mdl, anl, anlDerivative,residuals,gauss_fnd
    
    def model_fitting(self, gauss_now, image):
        """Use least squares to find the optimal direction for the next guess.."""
        delta_gauss = self.least_squares(gauss_now, image)
        """start of parabolic interpolation, used for refining an initial model fitting ie least_sq"""
        parabola_error_values = []
        """
        Want to compute the error function ie mean sq at three different points around the 
        current estimate. Decided from analyzing that this should probably not be hardcoded, 
        since the idea of parabola fitting is to use the CURRENT ESTIMATE given by least squares.
        Did not have time to fix it, but keep in mind that the choice of these three points
        majorly affect the fitting. 
        """

        three_different_points = [0, 0.6, 1.2] 
        
        """Compute the error function ie mean squared error to get the error values for
        the three different points"""
        for one_point in three_different_points:
            parabola_error_values.append(self.find_mean_sq(image, gauss_now.add(delta_gauss, one_point)))

        mean_sq_dof_one, mean_sq_dof_two, mean_sq_dof_three = parabola_error_values
        point_one, point_two, point_three = three_different_points

        """
        For a parabola on the form: f(x) = a*x^2 + b*x + c, fitted through three points, a and b can be calculated with 
        the following generalised formulas. 
        """
        a = (point_three * (mean_sq_dof_two - mean_sq_dof_one) + point_two * (
                mean_sq_dof_one - mean_sq_dof_three) + point_one * (mean_sq_dof_three - mean_sq_dof_two)) / (
                    (point_one - point_two) * (point_one - point_three) * (point_two - point_three))
        b = (point_one ** 2 * (mean_sq_dof_two - mean_sq_dof_three) + point_three ** 2 * (
                mean_sq_dof_one - mean_sq_dof_two) + point_two ** 2 * (
                     mean_sq_dof_three - mean_sq_dof_one)) / ((point_one - point_two) * (point_one - point_three) * (
                point_two - point_three))
        "Finding the vertex aka the extreme (minimum value) of a parabola is situated at this point:"
        point_min = np.divide(- b, (2 * a))  
               
        if point_min > 2:
            point_min = 2
        elif point_min < -1:
            point_min = -1

        gauss_min = gauss_now.add(delta_gauss, point_min)
        min_mean = self.find_mean_sq(image, gauss_min)            #Mean squared error

        #a and b should not be negative since that gives a non elliptical gaussian distribution function
        for gauss, gauss_n in zip(gauss_min, gauss_now):
            if gauss.a <= 0:
                gauss.a = gauss_n.a * 0.5
            if gauss.b <= 0:
                gauss.b = gauss_n.b * 0.5

        gauss_now = gauss_min
        mean_sq_dof = min_mean

        return gauss_now, mean_sq_dof, point_min

    def initial_guess(self, x0, y0, image):
        gauss_now = Gaussian()
        gauss_now.theta = 0
        gauss_now.amp = image[y0][x0]                      
        gauss_now.x0 = x0                                  #These are pixel values, offset from origin of image bottom left corner
        gauss_now.y0 = y0                                  #These are pixel values, offset from origin of image bottom left corner
        a_nominator = image[gauss_now.y0][gauss_now.x0 + 1] + image[gauss_now.y0][gauss_now.x0 - 1]
        b_nominator = image[gauss_now.y0 + 1][gauss_now.x0] + image[gauss_now.y0 - 1][gauss_now.x0]

        """
        In certain cases some of the values next to the max point (x0, y0), can be zero. The algebraic estimate of a and
        b uses log() so input shouldn't be zero. The value of a or b is then instead set to 10000, equating to very 
        rapid decline of intensity in the direction of a or b.
        """
        gauss_now.a = 10000 if a_nominator <= 0 else -log(a_nominator / (2 * gauss_now.amp))
        gauss_now.b = 10000 if b_nominator <= 0 else -log(b_nominator / (2 * gauss_now.amp))

        return gauss_now

    def find_mean_sq(self, image_data, gauss_now):
        """
        Finds and returns the mean squared error of the current gaussian model evaluated in every point in the image.
        """
        full_size = image_data.shape[0]
        num_points = full_size * full_size #tot number of pixels in the image
        mdl_data = gauss_now.build_image()
        squared_error= np.sum(np.sum(np.power((np.subtract(image_data, mdl_data)), 2)))
        #returning the mean squared error aka squared error/tot num of pixels in the image
        return squared_error / num_points

    def least_squares(self, gauss_in, image):
        """
        Calculates and sets delta_gauss which contains the adjustment for gauss_now based on its partial derivatives.
        least_squares minimizes and finds the next optimal path with the function ((J^T)*J)^-1*(J^T)*r
        """
        size = image.shape[0]
        num_peaks = len(gauss_in)
        delta_gauss = GaussList(size=size)
        """Get the normal vector (J^T)*r and the normal matrix (J^T)*J from make_norm_np"""
        norm_vec, norm_mat = self.make_norm_np(gauss_in, image, num_peaks)

        if scipy.linalg.det(norm_mat) == 0:
            for i in range(len(norm_mat)):
                if norm_mat[i][i] < 10e-8:
                    norm_mat[i][i] = 1
        """Get the inverse of the normal matrix in order to compute the minimazing function ((J^T)*J)^-1*(J^T)*r"""
        norm_inv = np.linalg.inv(norm_mat)
        adjust = np.matmul(norm_inv, norm_vec)
        """Now that we have what we need to adjust with in order to move toward the minimum of the function, 
        we adjust/update all of the current parameter values with them. We do this for all gaussian models,
        aka range(num_peaks) """
        for i in range(num_peaks):
            d_gauss = Gaussian()
            d_gauss.amp = adjust[0 + 6 * i]
            # Keep from making large adjustments in amplitude
            if abs(d_gauss.amp / gauss_in[i].amp) > 0.05:
                d_gauss.amp = 0.1 * gauss_in[i].amp if d_gauss.amp > 0 else -0.1 * gauss_in[i].amp
            d_gauss.a = adjust[1 + 6 * i]
            d_gauss.b = adjust[2 + 6 * i]
            d_gauss.x0 = adjust[3 + 6 * i]
            d_gauss.y0 = adjust[4 + 6 * i]
            # Keep from making large adjustments in offset
            if abs(d_gauss.x0 / image.shape[0]) > 0.2:
                d_gauss.x0 = 0.2 * gauss_in[i].x0 if d_gauss.x0 > 0 else -0.2 * gauss_in[i].x0
            if abs(d_gauss.y0 / image.shape[0]) > 0.2:
                d_gauss.y0 = 0.2 * gauss_in[i].y0 if d_gauss.y0 > 0 else -0.2 * gauss_in[i].y0
            d_gauss.theta = adjust[5 + 6 * i]
            delta_gauss.append(d_gauss)
        return delta_gauss
    
    @staticmethod
    def calculate_group_delay_rms(header, anlDerivative):
        """
        Function that fits a plane to the Group delay and returns the residual between the plane and the group delay.
        The root mean square error (RMS-error) is calculated for the residual with baselines smaller than the radius of the earth.
        """
        u = header.u.flatten()
        v = header.v.flatten()

        A = np.c_[u, v, np.ones(v.shape[0])]

        params, _, _, _ = np.linalg.lstsq(A, anlDerivative.flatten(), rcond=None) #Fitting plane to anlDerivative
        a, b, c = params

        fitted_plane = a*header.u + b * header.v + c

        scale_factor_residual = 1e12
        residuals = anlDerivative*scale_factor_residual - fitted_plane*scale_factor_residual

        #Extracting baselines smaller than radius of earth
        distance = sqrt((header.u - header.reference_pixel_RA) ** 2 + (header.v - header.reference_pixel_DEC) ** 2)
        mask = distance <= header.v[-1,0]

        rms_value = sqrt(np.mean(np.square(residuals[mask])))    #Calculating rms-error for baslines smaller than radius of earth

        return residuals, rms_value

    @staticmethod
    def precision(img_data, mdl_data):
        """
        Calculate the difference between the original and modelled image.
        :param img_data: np.array of original image
        :param mdl_data: np.array of modelled image
        :return: precisions - the different measurements of error as a dictionary where the type of precision measure
        is the key.
        """
        mean_sq = np.sum(np.sum(np.power(np.subtract(img_data, mdl_data), 2)))
        num_points = img_data.shape[0] ** 2

        precisions = {'mean error': sqrt(mean_sq) / num_points, 'mean error as % of max val':
            (sqrt(mean_sq) / num_points) / amax(img_data) * 100, 'total error': sqrt(mean_sq),
                      'total error as % of max val': sqrt(mean_sq) / amax(img_data) * 100}

        return precisions

    def make_norm_np(self, gauss_in, image, num_peak):
        """make_norm_np uses get_terms to get the partial derivatives needed in order to
        create the normal vector (J*T)*r where J*T is the transposed jacobian and r is the residual, as well
        as the normal matrix (J*T)*J. These are needed to compute the least squares, where the inverse of the 
        normal matrix times the normal vector gives us what we will add to the current parameters to move
        towards the minimum of the objective function"""
        norm_vec = np.zeros(6 * num_peak, dtype=float)
        norm_mat = np.zeros((6 * num_peak, 6 * num_peak))
        for index, gauss in enumerate(gauss_in, 0):
            terms = np.fromfunction(gauss.get_terms, (image.shape[0], image.shape[1]), dtype=float) #
            term_a, term_b = terms[::][0], terms[::][1] 
            g = abs(gauss.a * term_a + gauss.b * term_b)
            g[g > 50.] = 50.                                #Any value in array larger than 50 is set to 50
            exp_minus_g = exp(-g)

            apriori = exp_minus_g * gauss.amp
            residuals = image - apriori
            """These are the partial derivatives based on all of the parameters, beginning with A"""
            partials = [exp_minus_g,
                        -apriori * term_a,
                        -apriori * term_b,
                        apriori * (gauss.a * terms[::][2] + gauss.b * terms[::][3]),
                        apriori * (gauss.a * terms[::][4] + gauss.b * terms[::][5]),
                        -apriori * (gauss.b - gauss.a) * terms[::][6]
                        ]

            for i in range(6):
                norm_vec[index * 6 + i] = np.sum(residuals * partials[i])
                for j in range(6):
                    norm_mat[index * 6 + i][index * 6 + j] = np.sum(partials[i] * partials[j])

        return norm_vec, norm_mat

    def log_initial_guess(self, logger, mean_sq_dof, num_peak, guess, gauss_now,header):
        logger.info(f'Initial Guess #{guess + 1} \n')
        logger.info(f"Mean_sq_DOF: {round(mean_sq_dof, 10):<15} #Found peaks: {num_peak} \n")
        logger.info(f"{'#':<15}{'Amplitude':<15}{'X0':<15}{'Y0':<15}{'Sigma_X':<15}{'Sigma_Y':<15}"
                    f"{'Rotation':<15}\n")
        for index,gauss in enumerate(gauss_now):
            logger.info(f"{index:<15}{str(gauss)}")
        logger.info(f'{line}\n')

    def log_iteration(self, logger, mean_sq_dof, delta_mean, point_min, iteration, gauss_now,header):
        logger.info(f'Iteration #{iteration} \n')
        logger.info(f'Mean_sq_DOF: {round(mean_sq_dof, 10):<15} Delta_mean: {round(delta_mean, 10):<15}'
                     f'Point: {round(point_min, 5):<15} #Found peaks: {len(gauss_now)} \n')
        logger.info(f"{'#':<15}{'Amplitude':<15}{'X0':<15}{'Y0':<15}{'a':<15}{'b':<15}"
                     f"{'Rotation':<15}\n")
        logger.info("Found components: \n")
        for index, gauss in enumerate(gauss_now):
            logger.info(f"{index:<15}{str(gauss)}")
        logger.info(f'{line}\n')

    def log_final_iteration(self, logger, gauss_fnd,header):
        logger.info(f'{stars}FINAL{stars}\n')
        logger.info("Found components: \n")

        logger.info(f"{'#':<15}{'Amplitude':<15}{'X0':<15}{'Y0':<15}{'a':<15}{'b':<15}"
                     f"{'Rotation':<15}\n")

        for index, gauss in enumerate(gauss_fnd):
            #Code to print scaled components in log
            logger.info(f"{index:<15}{str(gauss)}")

        logger.info("Found components scaled according to physical coordinate system: \n")
        for index, gauss in enumerate(gauss_fnd):
            gauss_temp = [gauss.amp,(gauss.x0 - header.reference_pixel_RA) * header.pixel_increment_RA,(gauss.y0 - header.reference_pixel_DEC) * header.pixel_increment_DEC ,gauss.a/header.pixel_increment_RA**2, gauss.b/header.pixel_increment_RA**2, gauss.theta]
            formatted_values = "\t".join(f"{value:<15.3e}" if (isinstance(value, float) and (abs(value) < 1e-3 or abs(value) >= 1e6)) else f"{value:<15.5f}"for value in gauss_temp)            
            logger.info(f"{index:<15}{formatted_values}")

        logger.info("\n")