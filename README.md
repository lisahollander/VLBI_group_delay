# VLBI_group_delay
Summer internship project for NVI inc

## Downloading the Program
Download the github repository and then install the requirements. 

```bash
git clone https://github.com/lisahollander/VLBI_group_delay.git
pip install -r VLBI_group_delay/requirements.txt
```

## Running the Program
The program was constructed using visual studio code, but it can run from any editor or from the terminal. When running the
program from an editor, run the `__main__.py` file if you want to run the GUI. If not, use the terminal to run the script.

Example for using the terminal with a fits file from the fits folder (Note: You must be in the quasarModel directory for this command)
```bash
python quasarModel -r fits/fits_files/J0004-4736/J0004-4736_S_2017_07_17_pet_map.fits
```

Example for using the terminal with the sourcss predefined in constants.py, this can be useful if you want to run the script quickly (Note: You must be in the QuasarModel directory for this)
```bash
python quasarModel -r 1
```
will run source 1 in constants.py. 
