# pyppe
Python code computing the probability of planet-planet eclipses (PPE)

### Description
This code computes the probability of planet-planet eclipses from a population of objects defined in a CSV file.
The format of the CSV file is the same as the one provided by the [NASA Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu/) when downloading the table of planetary systems (e.g. [confirmed planets](https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=PS)).

### Dependencies
The code requires the following packages: `pickle`, `numpy`, `astropy`, `pandas`, `multiprocessing`

### How to run pyppe
To run the code, one must first modify the parameter file **PARAMS.txt** accordingly, and then type `python3 main.py PARAMS.txt` in a terminal window.
*Note that the typical computation time for a single planet pair is of several minutes.*

To create graphs showing the probability of PPE events as a function of time, type `python3 createGraphs.py <inputDir> <outputDir>` in a terminal window.
