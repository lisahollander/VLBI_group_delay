#A script that you can run if you want to create a csv file from a fits file. Useful if you want to look closer on the values without printing. 
import os
from astropy.io import fits
import pandas as pd

# Function to process FITS file and save data to CSV for easier viewing and analysis
def fits_to_csv(fits_file, output_dir):
    with fits.open(fits_file) as hdulist:
        hdulist.verify('fix') #need to use fix or else it complains about RA and DEC, built in astropy function
        for i, hdu in enumerate(hdulist): #loops through all hdus in the fits file, usually 2 (PrimaryHDU and HeaderHDU) but could be more
            # Create dictionaries to store header information
            header_dict = {key: value for key, value in hdu.header.items()}

            # Save header information to a DataFrame and then to a CSV file
            df_header = pd.DataFrame(list(header_dict.items()), columns=['Key', 'Value'])
            header_csv_file = os.path.join(output_dir, f'header_values_{i}.csv')
            df_header.to_csv(header_csv_file, index=False)
            print(f"Saved header information to {header_csv_file}")

            # Process and save data if available
            data = hdu.data
            if data is not None:
                if isinstance(hdu, (fits.BinTableHDU, fits.TableHDU)):
                    # Convert table data to a DataFrame with the column names specified in the fits file, if there are some
                    if hasattr(data, 'columns'):
                        df_data = pd.DataFrame(data, columns=data.columns.names)
                    else:
                        # If there is no column names defined in the fits file, just pick generic ones
                        df_data = pd.DataFrame(data)
                    data_csv_file = os.path.join(output_dir, f'data_{i}.csv')
                    df_data.to_csv(data_csv_file, index=False)
                    print(f"Saved table data to {data_csv_file}")

                else:
                    # Handle different shapes of image data
                    if data.ndim == 4:
                        # Slice off the first two dimensions and create a DataFrame from the resulting 2D array
                        first_two_dims_data = data[0, 0, :, :] #This is usually 512, 512 aka NAXIS1 and NAXIS2
                        df_first_two_dims = pd.DataFrame(first_two_dims_data)
                        first_two_dims_csv_file = os.path.join(output_dir, f'first_two_dims_data_{i}.csv')
                        df_first_two_dims.to_csv(first_two_dims_csv_file, index=False)
                        print(f"Saved first two dimensions data to {first_two_dims_csv_file}")

                        # Slice off the second two dimensions and create a DataFrame from the resulting 2D array
                        second_two_dims_data = data[:, :, 0, 0] #this is usually 1,1 aka NAXIS3 and NAXIS4
                        df_second_two_dims = pd.DataFrame(second_two_dims_data)
                        second_two_dims_csv_file = os.path.join(output_dir, f'second_two_dims_data_{i}.csv')
                        df_second_two_dims.to_csv(second_two_dims_csv_file, index=False)
                        print(f"Saved second two dimensions data to {second_two_dims_csv_file}")

                    elif data.ndim == 3:
                        # Handle 3D data if needed (e.g., slice or flatten)
                        df_data = pd.DataFrame(data.flatten())
                        data_csv_file = os.path.join(output_dir, f'data_{i}.csv')
                        df_data.to_csv(data_csv_file, index=False)
                        print(f"Saved 3D data to {data_csv_file}")

                    else:
                        # For other shapes or 2D data#
                        df_data = pd.DataFrame(data.flatten())
                        data_csv_file = os.path.join(output_dir, f'data_{i}.csv')
                        df_data.to_csv(data_csv_file, index=False)
                        print(f"Saved 2D data to {data_csv_file}")

            else:
                print(f"No data found in HDU {i}")

# Specify the path to the FITS file and the output directory
fits_file = r'fits\J0555+3922\J0555+3922_S_2018_07_05_pet_map.fits'
output_dir = 'csv_files'

#Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Call the function to convert FITS to CSV
fits_to_csv(fits_file, output_dir)

# Construct CSV file paths, 1/2 hardcoded, update as needed
primary_data_csv = os.path.join(output_dir, 'first_two_dims_data_0.csv')  # Update if needed
header_csv = os.path.join(output_dir, 'header_values_0.csv')  # Update if needed

# Debug prints
print(f"Primary Data CSV Path: {primary_data_csv}")
print(f"Header CSV Path: {header_csv}")

# Check if files exist
if not os.path.isfile(primary_data_csv):
    print(f"File not found: {primary_data_csv}")

if not os.path.isfile(header_csv):
    print(f"File not found: {header_csv}")
