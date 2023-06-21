import numpy as np
import matplotlib.pyplot as plt
import pickle
import os
import sys
from astropy import time as tm

if len(sys.argv) != 3:
	print("Usage :", sys.argv[0], "<inputDir> <outputDir>")
	quit(1)

inputDir  = sys.argv[1]
outputDir = sys.argv[2]
overwrite = False

def create_graph_path(datafName, outputDir):
	if outputDir[-1] == '/':
		return outputDir + datafName.replace(".pkl", "_graph.png")
	else:
		return outputDir + '/' + datafName.replace(".pkl", "_graph.png")


def collect_data_from_file(fName):
	infos = {'st_name'  : fName.split('/')[-1].split('_')[0],                      \
		 'pl1_name' : fName.split('/')[-1].split('_')[1],                      \
		 'pl2_name' : fName.split('/')[-1].split('_')[2],                      \
		 'l1'       : fName.split('/')[-1].split('_')[3],                      \
		 'l2'       : fName.split('/')[-1].split('_')[4].replace('.pkl', ''),  \
		 'fName'    : fName}
	    
	if infos['l1'][0] == 'm':
	    infos['l1'] = -1 * float(infos['l1'].replace('m', ''))
	    
	if infos['l2'][0] == 'm':
	    infos['l2'] = float(infos['l2'].replace('m', ''))
	
	with open(inputDir + "/" +infos['fName'], 'rb') as inFile:
	    infos['data'] = pickle.load(inFile)

	return infos


def create_graph(infos):
	
	# Creating the new x axis
	xaxis = [str(tm.Time(tm.Time(_, format='jd'), format='isot')).split('T')[0] for _ in infos['data']['Tmid']]

	fig, ax = plt.subplots(figsize=(8, 4), dpi=300)
	ax.set_title("Star {}  -  {}&{}  -  $\lambda_1={}$, $\lambda_2={}$".format(infos['st_name'], infos['pl1_name'], infos['pl2_name'], infos['l1'], infos['l2']))
	ax.set_xlabel("Mid transit time")
	ax.set_ylabel("Proba. (log)")
	
	ax.set_xticks(infos['data']['Tmid'])
	ax.set_xticklabels(xaxis, rotation = 30)
	
	ax.bar(infos['data']['Tmid'], infos['data']['prob'], width=1)
	ax.scatter(infos['data']['Tmid'], infos['data']['prob'], c='C0')
	ax.set_yscale('log')

	fig.savefig(infos['graph_path'], facecolor='w', bbox_inches='tight')
	plt.close(fig)


counter = 1
numgraphs = len(os.listdir(inputDir))


for fName in os.listdir(inputDir):

	if '.pkl' not in fName:
		continue
	
	if os.path.isfile(create_graph_path(fName, outputDir)) and not overwrite:
		print("Warning : File already exists and overwriting is not allowed :", fName)
		continue
	else:
		print("Creating graph file {:<60} from file {}".format(create_graph_path(fName, outputDir), fName))
		continue

#	print("{} / {}".format(counter, numgraphs))
#	infos = collect_data_from_file(fName)	
#	infos['graph_path'] = create_graph_path(infos['fName'], outputDir)
#	create_graph(infos)
#	counter += 1
    
