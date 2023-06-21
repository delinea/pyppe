import os
import numpy as np
import pandas as pd
import sys
from multiprocessing import Pool
import functions as func

if __name__ == '__main__':

    # =========================================================================
    # Parsing file PARAMS.txt
    # =========================================================================

    # Checking number of parameters
    if len(sys.argv) != 2:
        print("Usage : python3", sys.argv[0], "<parameter_file_name>")
        quit(1)

    paramFileName = sys.argv[1]
    params = {}

    print("\n----- Parsing file", paramFileName, "...")

    with open(paramFileName, 'r') as inFile:

        content = inFile.readlines()

        for line in content:

            if line.replace(' ', '')[0] != '#' and line not in ['', '\n']:
                line = line.replace(' ', '').split("#")[0].split("=")

                if line[1].isnumeric():
                    params[line[0]] = int(line[1])

                elif line[1].replace('.', '').isnumeric():
                    params[line[0]] = float(line[1])

                else:
                    params[line[0]] = line[1]

    # Checking SELECTION entry
    if params['SELECTION'] not in ['LAMBDA_ALL', 'LAMBDA_NOTALL',
                                   'LAMBDA_NONE']:
        raise ValueError("SELECTION entry is not valid in parameter file ({})"
                         " : {} (LAMBDA_ALL / LAMBDA_NOTALL / LAMBDA_NONE)"
                         .format(paramFileName, params['SELECTION']))

    # Checking SAVE_IN_FILE entry
    if params['SAVE_IN_FILE'] in ['True', 'TRUE']:
        params['SAVE_IN_FILE'] = True
    elif params['SAVE_IN_FILE'] in ['False', 'FALSE']:
        params['SAVE_IN_FILE'] = False
    else:
        raise ValueError("SAVE_IN_FILE entry is not valid in parameter file "
                         "({}) : {} (TRUE / FALSE)"
                         .format(paramFileName, params['SAVE_IN_FILE']))

    # Checking OVERWRITE entry
    if params['OVERWRITE'] in ['True', 'TRUE']:
        params['OVERWRITE'] = True
    elif params['OVERWRITE'] in ['False', 'FALSE']:
        params['OVERWRITE'] = False
    else:
        raise ValueError("OVERWRITE entry is not valid in parameter file ({})"
                         " : {} (TRUE / FALSE)"
                         .format(paramFileName, params['OVERWRITE']))

    # Checking VERBOSE entry
    if params['VERBOSE'] in ['True', 'TRUE']:
        params['VERBOSE'] = True
    elif params['VERBOSE'] in ['False', 'FALSE']:
        params['VERBOSE'] = False
    else:
        raise ValueError("VERBOSE entry is not valid in parameter file ({}) :"
                         " {} (TRUE / FALSE)"
                         .format(paramFileName, params['VERBOSE']))

    # Checking if output directory exists
    if params['SAVE_IN_FILE'] and not os.path.isdir(params['OUTPUT_DIR']):
        os.mkdir(params['OUTPUT_DIR'])

    # Printing parameters
    print("\nParameters :\n")

    for _ in zip(params.keys(), params.values()):
        print("{:<15} | {:<50}".format(_[0], _[1]))

    # =========================================================================
    # Initializing data
    # =========================================================================

    print("\n----- Collecting data ...\n")

    # Creating a list of dictionnaries
    data = []
    df = pd.read_csv(params['DATAFILE_NAME'])

    for idx in df.index:
        data.append(dict(df.iloc[idx]))

    # We add 'st_name' (name of the star) for each planet
    for pl in data:
        if '.' in pl['pl_name']:
            pl['st_name'] = pl['pl_name'].split('.')[0]
        else:
            pl['st_name'] = ' '.join(pl['pl_name'].split(' ')[0:-1])

    print("Collected", len(data), "planets (",
          func.get_number_of_stars(data), "stars )\n")

    # =========================================================================
    # Filtering data
    # =========================================================================

    print("----- Filtering data ...\n")

    numPlanetsInit = len(data)

    print("\nInitial number of planets :", numPlanetsInit,
          "({} stars)".format(func.get_number_of_stars(data)))

    print("\nRemoving planets with no radius or error on radius ...")
    counter = 0

    pl_to_remove = []
    for pl in data:
        if ((np.isnan(pl['pl_radius'])) or (np.isnan(pl['pl_radius_err+']))
                or (np.isnan(pl['pl_radius_err-']))):
            pl_to_remove.append(pl)
            counter += 1

    print("    ->", counter, "/", numPlanetsInit,
          "planets where the radius isn't known")

    counter = 0

    print("Removing planets with eccentricity >", params['ECC_MAX'], "...")
    for pl in data:
        if pl['ecc'] > params['ECC_MAX']:
            pl_to_remove.append(pl)
            counter += 1
    print("    ->", counter, "/", numPlanetsInit, "planets where ecc >",
          params['ECC_MAX'])

    # Remove duplicates in list, and remove planets
    pl_to_remove_noDouble = []
    for i in range(len(pl_to_remove)):
        if pl_to_remove[i] not in pl_to_remove[i + 1:]:
            pl_to_remove_noDouble.append(pl_to_remove[i])

    print("->", len(pl_to_remove_noDouble), "planets to remove.")

    for pl in pl_to_remove_noDouble:
        data.remove(pl)

    counter = 0
    counterStars = 0

    print("Removing systems with < 2 planets ...")
    av_stars = func.get_available_stars(data)

    for star in av_stars:
        pl_set = func.get_planets_set(data, star)

        if len(pl_set) < 2:
            counterStars += 1
            for pl in pl_set:
                data.remove(pl)
                counter += 1

    print("    ->", counterStars, "systems with < 2 planets (",
          counter, "planets )")

    print("\nPlanets remaining :", len(data), "(",
          func.get_number_of_stars(data), "stars )\n")

    # =========================================================================
    # Distributing data according to the obliquity (lambda)
    # =========================================================================

    lambdaall    = func.get_planets_with_lambda(data, True)
    lambdanotall = func.get_planets_with_lambda(data, False)
    lambdanone   = func.get_planets_without_lambda(data)

    print("\n--- Repartition of the system wrt obliquity :\n")
    print("Obliquity known for every planet    :",
          func.get_number_of_stars(lambdaall), "  systems,",
          len(lambdaall), "planets")
    print("Obliquity unknown ...")
    print("\t... for at least one planet :",
          func.get_number_of_stars(lambdanotall), " systems,",
          len(lambdanotall), "planets")
    print("\t... for every planet        :",
          func.get_number_of_stars(lambdanone), "systems,",
          len(lambdanone), "planets")
    print()
    print("Total             :", len(data))
    print("Sum of subsystems :",
          len(lambdaall) + len(lambdanotall) + len(lambdanone), "\n")

    if len(data) != len(lambdaall) + len(lambdanotall) + len(lambdanone):
        raise RuntimeError("Sum of planets in subsystems is not equal "
                           "to the total number of planets.")

    # =========================================================================
    # Displaying data information (known values)
    # =========================================================================

    numplanets = [len(func.get_planets_set(data, star))
                  for star in func.get_available_stars(data)]

    distrib = {}

    for n in numplanets:
        if str(n) not in distrib.keys():
            distrib[str(n)] = 1
        else:
            distrib[str(n)] += 1

    print("\n----- Repartition of the umber of planets around stars\n")

    for k in sorted(distrib.keys()):
        print("{:>2} planets : {:>5} / {}"
              .format(k, distrib[k], func.get_number_of_stars(data)))

    print("\n----- Known values in data :\n")

    histWidth = 80
    keys = list(data[0].keys())

    for k in keys:
        if (k != 'pl_name') and (k != 'detection') and (k != 'st_name'):
            vec = np.array([_[k] for _ in data])
            sel = np.isnan(vec)
            numKnownValues = len(vec[~sel])

            numHashtags = int(numKnownValues * histWidth / len(vec))

            print("{:>14} : {:>4} / {:<5}".format(k, numKnownValues,
                                                  len(vec)), end=' |')
            for _ in range(numHashtags):
                print('#', end='')
            for _ in range(histWidth - numHashtags):
                print('.', end='')
            print('|')

    print("\n\n")

    # =========================================================================
    # Creating parameters list for function compute_transit_probability()
    # =========================================================================

    Tstart     = params['T_START']
    Tend       = params['T_END']
    outputDir  = params['OUTPUT_DIR']
    saveFile   = params['SAVE_IN_FILE']
    overwrite  = params['OVERWRITE']
    verbose    = params['VERBOSE']
    outFileExt = '.pkl'

    selected_planets = []

    if params['SELECTION'] == 'LAMBDA_ALL':
        selected_planets = lambdaall
    elif params['SELECTION'] == 'LAMBDA_NOTALL':
        selected_planets = lambdanotall
    elif params['SELECTION'] == 'LAMBDA_NONE':
        selected_planets = lambdanone

    # We store the systems where the computation is not possible
    problematicSystems = []

    # arguments is a list, where each element is a tuple in the format
    # [..., (data, pl1, pl2, Tstart, Tend, lambdas, numsamples, ofName,
    #        saveInFile=True) , ...]
    arguments = []

    for star in func.get_available_stars(selected_planets):
        pl_set = func.get_planets_set(selected_planets, star)

        for i in range(len(pl_set)):
            pl1 = pl_set[i]

            for j in range(i + 1, len(pl_set)):
                pl2 = pl_set[j]

                if (np.isnan(pl1['period']) or np.isnan(pl1['T0'])
                        or np.isnan(pl1['T14']) or np.isnan(pl2['period'])
                        or np.isnan(pl2['T0']) or np.isnan(pl2['T14'])):
                    ps_j = {'{}_{}'.format(pl1['pl_name'], pl2['pl_name']):
                            "Unknown value"}
                    problematicSystems.append(ps_j)
                    continue

                planetsInverted = False

                if pl2['period'] > pl1['period']:
                    (pl1, pl2) = (pl2, pl1)
                    planetsInverted = True

                if np.isfinite(pl1['lambda']) and np.isfinite(pl2['lambda']):

                    l1      = pl1['lambda']
                    l1e     = func.compute_mean_error(data, pl1['pl_name'],
                                                      'lambda')
                    l2      = pl2['lambda']
                    l2e     = func.compute_mean_error(data, pl2['pl_name'],
                                                      'lambda')
                    lambdas = {'l1' : l1, 'l1e' : l1e, 'l2' : l2, 'l2e' : l2e}

                    outFileName = func.create_outFile_name(
                        pl1, pl2, l1, l2, outputDir, fileExtension=outFileExt)

                    arguments.append(tuple([data, pl1, pl2, Tstart, Tend,
                                            lambdas, params['NUM_SAMPLES'],
                                            outFileName, saveFile, overwrite,
                                            verbose]))

                elif np.isfinite(pl1['lambda']) and np.isnan(pl2['lambda']):

                    vec_lambda2 = pl1['lambda'] + np.array([-90, -45, 0,
                                                            45, 90])

                    for lambda2 in vec_lambda2:

                        l1      = pl1['lambda']
                        l1e     = func.compute_mean_error(data, pl1['pl_name'],
                                                          'lambda')
                        l2      = lambda2
                        l2e     = l1e
                        lambdas = {'l1' : l1, 'l1e' : l1e, 'l2' : l2,
                                   'l2e' : l2e}

                        outFileName = func.create_outFile_name(
                            pl1, pl2, l1, l2, outputDir,
                            fileExtension=outFileExt)

                        arguments.append(tuple([data, pl1, pl2, Tstart, Tend,
                                                lambdas, params['NUM_SAMPLES'],
                                                outFileName, saveFile,
                                                overwrite, verbose]))

                elif np.isnan(pl1['lambda']) and np.isfinite(pl2['lambda']):

                    vec_lambda1 = pl2['lambda'] + np.array([-90, -45, 0,
                                                            45, 90])

                    for lambda1 in vec_lambda1:

                        l2      = pl2['lambda']
                        l2e     = func.compute_mean_error(data, pl2['pl_name'],
                                                          'lambda')
                        l1      = lambda1
                        l1e     = l2e
                        lambdas = {'l1' : l1, 'l1e' : l1e, 'l2' : l2,
                                   'l2e' : l2e}

                        outFileName = func.create_outFile_name(
                            pl1, pl2, l1, l2, outputDir,
                            fileExtension=outFileExt)

                        arguments.append(tuple([data, pl1, pl2, Tstart, Tend,
                                                lambdas, params['NUM_SAMPLES'],
                                                outFileName, saveFile,
                                                overwrite, verbose]))

                elif np.isnan(pl1['lambda']) and np.isnan(pl2['lambda']):

                    lambda1 = 0.0
                    vec_lambda2 = lambda1 + np.array([-90, -45, 0, 45, 90])

                    for lambda2 in vec_lambda2:

                        l1      = lambda1
                        l1e     = 10
                        l2      = lambda2
                        l2e     = l1e
                        lambdas = {'l1' : l1, 'l1e' : l1e, 'l2' : l2,
                                   'l2e' : l2e}

                        outFileName = func.create_outFile_name(
                            pl1, pl2, l1, l2, outputDir,
                            fileExtension=outFileExt)

                        arguments.append(tuple([data, pl1, pl2, Tstart, Tend,
                                                lambdas, params['NUM_SAMPLES'],
                                                outFileName, saveFile,
                                                overwrite, verbose]))

                else:
                    print("    !!!!! Nothing done for", pl1['pl_name'],
                          "and", pl2['pl_name'], "!!!!!")
                    continue

                if planetsInverted:
                    (pl1, pl2) = (pl2, pl1)

    if len(problematicSystems) != 0:
        print("\nThe folowwing systems are not usable :\n")
        [print(_) for _ in problematicSystems]

    ofNames = np.array([_[7] for _ in arguments])

    if len(ofNames) != np.unique(len(ofNames)):
        raise RuntimeError("There might be some duplicates in the file names.")

    # No parallelisation
    if params['NUM_CPU'] == 1:
        for arg in arguments:
            func.compute_transit_probability(*arg)

    # Parallelisation
    else:
        with Pool(processes=params['NUM_CPU']) as p:
            p.starmap(func.compute_transit_probability, arguments)
