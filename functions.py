import numpy as np
import pickle
import sys
import os
from astropy import time as tm
from astropy import constants as const


def get_available_stars(data):
    """
    Returns a list containing all the stars in the dataset

    PARAMETERS
        data : The dataset (list of dictionnaries)
    """
    return np.unique([_['st_name'] for _ in data])


def get_number_of_stars(data):
    """
    Returns the number of stars in the dataset

    PARAMETERS
        data : The dataset (list of dictionnaries)
    """
    return len(get_available_stars(data))


def get_planets_set(data, star):
    """
    Returns a list of dictionnaries containing all the planets orbiting
    around a given star

    PARAMETERS
        data : The dataset (list of dictionnaries)
        star : The name of the star
    """
    output = []

    for pl in data:
        if (star == pl['st_name']):
            output.append(pl)

    return output


def get_planet(data, pl_name):
    """
    Returns the information of a planet (dictionnary), given its name

    PARAMETERS
        data    : The dataset (list of dictionnaries)
        pl_name : The name of the planet

    RAISE
        Raises a ValueError if the planet is not found in the dataset
    """
    for pl in data:
        if pl['pl_name'] == pl_name:
            return pl

    raise ValueError("'{}' not found in database".format(pl_name))


def compute_mean_error(data, pl_name, key):  # OK
    """
    Return the average between sup and inf errors :
        sqrt( (err+^2 + err-^2) / 2 )

    PARAMETERS
        data    : The dataset (list of dictionnaries)
        pl_name : The name of the planet
        key     : The quantity on which to compute the error
    """

    pl = get_planet(data, pl_name)
    return np.sqrt(0.5 * (pl[key + '_err+']**2 + pl[key + '_err-']**2))


def compute_T14(pl):
    """
    Return the value of T14 according to the appoximation formula

    PARAMETERS
        pl : The planet (dictionnary)
    """

    P   = pl['period'] * 24           # Hours
    Rp  = pl['pl_radius'] / const.M_jup.value * const.M_earth.value
    Rs  = pl['st_radius'] * const.M_sun.value
    aRs = pl['aRs']
    i   = pl['incl']

    frac = (1 + Rp / Rs) / (aRs)
    return ((P / np.pi)
            * np.arcsin(np.sqrt(np.power(frac, 2)
                                - np.power(np.cos(np.radians(i)), 2))))


def get_planets_with_lambda(data, allWithLambda=False):
    """
    Returns the list of planets where (all or at least one of) the obliquity
    is known

    PARAMETERS
        data          : The dataset (list of dictionnaries)
        allwithLambda : A boolean (default False) to refine the selection :
                        True  will return the systems where the obliquity is
                            known for every planets
                        False will return the systems where at least one
                            obliquity is known (but not all)
    """
    output = []
    av_stars = get_available_stars(data)

    for i in range(len(av_stars)):
        star = av_stars[i]                    # Star name
        pl_set = get_planets_set(data, star)  # Planets orbiting around star
        pl_lambda = np.array([_['lambda'] for _ in pl_set])  # Obliquities

        if allWithLambda and np.all(np.isfinite(pl_lambda)):
            for pl in pl_set:
                output.append(pl)

        elif (not allWithLambda and np.any(np.isfinite(pl_lambda))
              and not np.all(np.isfinite(pl_lambda))):

            for pl in pl_set:
                output.append(pl)

    return output


def get_planets_without_lambda(data):
    """
    Returns the list of planets where the obliquity is not known for every
    planet

    PARAMETERS
        data        : The dataset (list of dictionnaries)
    """
    output = []
    av_stars = get_available_stars(data)

    for i in range(len(av_stars)):
        star = av_stars[i]
        pl_set = get_planets_set(data, star)
        pl_lambda = np.array([_['lambda'] for _ in pl_set])

        if all(np.isnan(pl_lambda)):
            for pl in pl_set:
                output.append(pl)

    return output


def get_transit_multiples(T0, P, Tstart, Tend):
    """
    Returns a list of integers n_i such that T0 + n_i * P is included between
    starting and ending times.

    PARAMETERS                                         UNIT    |  FORMAT
        T0     : Mid-transit time of the planet         JD     |  float  
        P      : Period of the planet                  days    |  float  
        Tstart : Starting time                      yyyy-mm-dd |  string 
        Tend   : Ending time                        yyyy-mm-dd |  string 
    """

    if P <= 0.0:
        raise ValueError("Period should be > 0.0.")

    output = []

    # Convert Tstart and Tend in Julian days (float numbers)
    Tstart = float(tm.Time(tm.Time(Tstart), format='jd').value)
    Tend   = float(tm.Time(tm.Time(Tend), format='jd').value)

    # Invert if Tstart < Tend
    if Tstart > Tend:
        (Tstart, Tend) = (Tend, Tstart)

    multiple = 0

    while T0 + multiple * P <= Tend:
        if T0 + multiple * P >= Tstart:
            output.append(multiple)
        multiple += 1

    return output


def compute_transit_probability_fake(data, pl1, pl2, Tstart, Tend, lambdas,
                                     numsamples, ofName, saveInFile=True,
                                     overwrite=False, verbose=True):
    """
    Only for debugging purpose.
    """
    if verbose:
        print("    Computing transit probability for", pl1['pl_name'],
              "and", pl2['pl_name'], "lambda1 =", lambdas['l1'],
              ", lambda2 =", lambdas['l2'])
        print("    Saving data in file", ofName)

    if (saveInFile and os.path.isfile(ofName.replace('.pkl', '.txt'))
            and not overwrite):
        print("    ERROR : file already exists and overwriting is not "
              "allowed :", ofName.replace('.pkl', '.txt'))

    if saveInFile:
        with open(ofName.replace('.pkl', '.txt'), 'w') as outFile:
            outFile.write("This is file {}".format(ofName.replace('.pkl',
                                                                  '.txt')))
    print()


def compute_transit_probability(data, pl1, pl2, Tstart, Tend, lambdas,
                                numsamples, ofName, saveInFile=True,
                                overwrite=False, verbose=True):

    if verbose:
        print("Computing transit probability for {} & {} : l1 = {}, l2 = {}, "
              "ofName = {}".format(pl1['pl_name'], pl2['pl_name'],
                                   lambdas['l1'], lambdas['l2'], ofName))

    if saveInFile and os.path.isfile(ofName) and not overwrite:
        print("WARNING (from compute_transit_probability()) : Overwriting is "
              "not allowed and file already exists :", ofName, file=sys.stderr)
        return

    P1T14 = pl1['T14']

    if np.isnan(P1T14):
        P1T14 = compute_T14(pl1)

    # Ephemeris for choosing possible observing windows
    P10 = pl1['period']    # Period of planet c
    # Mid-transit time of planet c (Julian date = JD)
    T10 = pl1['T0']

    if np.isnan(pl1['period']) or np.isnan(pl1['T0']) or np.isnan(P1T14):
        print("Error (from compute_transit_probability()) :"
              " One of this value isn't defined. period :", pl1['period'],
              "T0 :", pl1['T0'], "T14 :", P1T14)
        print("The function will return None")
        return None

    # Minimum duration of the conjonction for detection
    # threshold above which one considers there was an eclipse [minutes]
    min_duration_m = 1.0

    # Observing window around mid-transit of planet c (longest period)

    Dt   = 0.6 * P1T14  # garder en jours
    step = min_duration_m/(2*24*60)  # Steps of 30 seconds
    dt   = np.arange(-Dt, Dt, step)
    nt   = len(dt)

    # Selected transits
    # All transits within a given time range
    transits = np.array(get_transit_multiples(pl1['T0'], pl1['period'],
                                              Tstart, Tend))

    ntransit = len(transits)
    Tmid = T10 + transits * P10

    # Full array of observing times
    # Each row is the transit time in the time interval dt
    t = Tmid[:, None] + dt[None, :]

    # Number of samples from fake MCMC
    nsample = numsamples

    # Array with 1 for conjonction, 0 otherwise
    # Transit probabilities
    conj = np.zeros((nsample, ntransit, nt), dtype=bool)

    # IAU convention - Resolution B3
    # http://arxiv.org/abs/1510.07674
    # No need to change.
    rSun_m      = 6.957e8
    AU_m        = 149597870700.
    rSun_AU     = rSun_m / AU_m
    GmSun_m3s2  = 1.3271244e20  # m^3/s2
    GmSun_AU3d2 = GmSun_m3s2 / AU_m**3 * (3600 * 24)**2

    # Compute errors ----------------------------------------------------------
    st_radius_err  = compute_mean_error(data, pl1['pl_name'], 'st_radius')
    st_mass_err    = compute_mean_error(data, pl1['pl_name'], 'st_mass')
    pl2_period_err = compute_mean_error(data, pl2['pl_name'], 'period')
    pl2_Tb_err     = compute_mean_error(data, pl2['pl_name'], 'T0')
    pl2_rRs_err    = compute_mean_error(data, pl2['pl_name'], 'RpRs')
    pl2_ib_err     = compute_mean_error(data, pl2['pl_name'], 'incl')
    pl2_Wb_err     = lambdas['l2e']
    pl1_period_err = compute_mean_error(data, pl1['pl_name'], 'period')
    pl1_Tc_err     = compute_mean_error(data, pl1['pl_name'], 'T0')
    pl1_Rcs_err    = compute_mean_error(data, pl1['pl_name'], 'RpRs')
    pl1_ic_err     = compute_mean_error(data, pl1['pl_name'], 'incl')
    pl1_Wc_err     = lambdas['l1e']
    # -------------------------------------------------------------------------

    # Draw samples from fake MCMC ---------------------------------------------
    # Star
    Rs  = np.random.normal(pl1['st_radius'], st_radius_err,
                           size=nsample) * rSun_AU
    ms  = np.random.normal(pl1['st_mass'], st_mass_err, size=nsample)
    # Planet 2
    Pb  = np.random.normal(pl2['period'], pl2_period_err, size=nsample)
    # Temps mi-transit
    Tb  = np.random.normal(pl2['T0'], pl2_Tb_err, size=nsample)
    # Rayon de la planète
    Rbs = np.random.normal(pl2['RpRs'], pl2_rRs_err, size=nsample)
    ib  = np.random.normal(pl2['incl'], pl2_ib_err, size=nsample) * np.pi/180
    Wb  = np.random.normal(lambdas['l2'], pl2_Wb_err, size=nsample) * np.pi/180
    # Planet 1
    Pc  = np.random.normal(pl1['period'], pl1_period_err, size=nsample)
    Tc  = np.random.normal(pl1['T0'], pl1_Tc_err, size=nsample)
    Rcs = np.random.normal(pl1['RpRs'], pl1_Rcs_err, size=nsample)
    Wc  = np.random.normal(lambdas['l1'], pl1_Wc_err,size=nsample) * np.pi/180
    ic  = np.random.normal(pl1['incl'], pl1_ic_err, size=nsample) * np.pi/180
    # -------------------------------------------------------------------------

    # Iterate over MCMC samples
    for i in range(nsample):

        # Compute cartesian coordinates
        Rb = Rs[i] * Rbs[i]
        ab = (GmSun_AU3d2 * ms[i] * (Pb[i] / (2 * np.pi)) ** 2) ** (1. / 3)
        Mb = 2 * np.pi * (t - Tb[i]) / Pb[i] + np.pi / 2
        Xb = ab * np.cos(Mb)
        Yb = ab * np.sin(Mb)
        xb = np.cos(Wb[i]) * Xb - np.sin(Wb[i]) * np.cos(ib[i]) * Yb
        yb = np.sin(Wb[i]) * Xb + np.cos(Wb[i]) * np.cos(ib[i]) * Yb
        zb = np.sin(ib[i]) * Yb

        Rc = Rs[i] * Rcs[i]
        ac = (GmSun_AU3d2 * ms[i] * (Pc[i] / (2 * np.pi)) ** 2) ** (1. / 3)
        Mc = 2 * np.pi * (t - Tc[i]) / Pc[i] + np.pi / 2
        Xc = ac * np.cos(Mc)
        Yc = ac * np.sin(Mc)
        xc = np.cos(Wc[i]) * Xc - np.sin(Wc[i]) * np.cos(ic[i]) * Yc
        yc = np.sin(Wc[i]) * Xc + np.cos(Wc[i]) * np.cos(ic[i]) * Yc
        zc = np.sin(ic[i]) * Yc

        # Compute projected distances between the 2 objects
        db  = np.sqrt(xb**2 + yb**2)
        dc  = np.sqrt(xc**2 + yc**2)
        dbc = np.sqrt((xb - xc)**2 + (yb - yc)**2)

        # Check for conjonction
        conj[i] = ((db < Rs[i]) & (dc < Rs[i])  # (1)
                   & (dbc < Rb + Rc)            # (2)
                   & (zb > 0) & (zc > 0))       # (3)
        # (1) - both planets on the stellar disc
        # (2) - Distance between 2 planets < radii sum  ->  overlap
        # (3) - both planets in front of the star (and not behind)

    # Duration of the conjonction event
    duration_m = np.sum(conj, axis=2) * (dt[1] - dt[0]) * 24 * 60

    # Detection if duration > min_duration
    detect = (duration_m > min_duration_m)

    detect_sample = np.sum(detect, axis=0)
    detect_proba = detect_sample / nsample
    total_dur = np.sum(duration_m, axis=0)
    av_dur = np.zeros(np.shape(detect_sample))
    idx_detect = detect_sample > 0
    av_dur[idx_detect] = total_dur[idx_detect]/detect_sample[idx_detect]

    output = {'Tmid'   : Tmid,   'prob'     : detect_proba,
              'av_dur' : av_dur, 'detect'   : detect, 'ntransit' : ntransit}

    if saveInFile:
        # Check if file already exists
        if os.path.isfile(ofName) and not overwrite:
            print("WARNING (from compute_transit_probability()) : "
                  "Overwriting is not allowed and file already exists :",
                  ofName, file=sys.stderr)
            return

        with open(ofName, 'wb') as outFile:
            pickle.dump(output, outFile)


def create_outFile_name(pl1, pl2, l1, l2, outputDir, fileExtension='.pkl'):
    """
    Returns the file name in the format <st_name>_<pl1>_<pl2>_<l1>_<l2>.txt
    Where <pl1> and <pl2> are the corresponding letters (b, c, ...)
    and <l1>, <l2> are the obliquities. If an obliquity is negative,
    the letter m is added in front of the value.
    """
    if pl1['st_name'] != pl2['st_name']:
        raise RuntimeError("Planets don't seem to belong to the same system :",
                           pl1['pl_name'], pl2['pl_name'], ".")

    st_name = pl1['st_name'].replace(' ', '-')
    pl1Letter = None
    pl2Letter = None
    l1_str    = None
    l2_str    = None

    if '.' in st_name:
        pl1Letter = pl1['pl_name'].split('.')[1]
        pl2Letter = pl2['pl_name'].split('.')[1]
    else:
        pl1Letter = pl1['pl_name'].split(' ')[-1]
        pl2Letter = pl2['pl_name'].split(' ')[-1]

    if l1 < 0:
        l1_str = 'm' + str(l1).replace('-', '')
    else:
        l1_str = str(l1)

    if l2 < 0:
        l2_str = 'm' + str(l2).replace('-', '')
    else:
        l2_str = str(l2)

    if fileExtension[0] == '.':
        filename = "{}_{}_{}_{}_{}{}".format(st_name, pl1Letter,
                                             pl2Letter, l1_str, l2_str,
                                             fileExtension)
    else:
        filename = "{}_{}_{}_{}_{}.{}".format(st_name, pl1Letter,
                                              pl2Letter, l1_str, l2_str,
                                              fileExtension)
    return os.path.join(outputDir, filename)
