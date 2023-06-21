# pyppe
Python code computing the probability of planet-planet eclipses (PPE)

This code computes the probability of planet-planet eclipses from a population of objects defined in a CSV file.
The format of the CSV file is the same as the one provided by the NASA Exoplanet Archive (https://exoplanetarchive.ipac.caltech.edu/) when downloading the table of planetary systems (e.g. confirmed planets - https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=PS).

The code requires the following packages: pickle, numpy, astropy, pandas, multiprocessing

To run the code, one must modify the parameter file 'PARAMS.txt' accordingly and then type 'python main.py PARAMS.txt'.
Note that the typical computation time for a single planet pair is of several minutes.
