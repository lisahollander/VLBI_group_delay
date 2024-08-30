# VLBI_group_delay
Summer internship project for NVI inc

## Downloading the Program
Downloading via pip is possible with 

```bash
pip install QuasarModel
```

and gives easy access to the special classes to build your own programs with. However, if you want to run the finished program with the gui or from the terminal it's better to download the github repository and then installing the requirements. 

```bash
git clone https://github.com/lisahollander/VLBI_group_delay.git
pip install -r QuasarModel/requirements.txt
```

## Running the Program
The program can run from an editor or from the terminal. When running the
program from an editor, run the `__main__.py` file if you want to run the GUI. If not, use the terminal to run the script and specify the fits file you would like to run. 

Example for using the terminal, not the GUI (Note: You must be in the quasarModel directory for this command)
```bash
python quasarModel -r fits/fits_files/J0004-4736/J0004-4736_S_2017_07_17_pet_map.fits
```

Example for using the terminal with the sourcss predefined in constants.py (Note: You must be in the QuasarModel directory for this)
```bash
python quasarModel -r 1
```
will run source 1 in constants.py. 
