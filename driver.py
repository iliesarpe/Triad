import sys
import json
import os
import subprocess
from multiprocessing import Pool
import h5py
from concurrent.futures import ThreadPoolExecutor
import math

def dumpJson(data, filename):
	json_object = json.dumps(data)
	with open(filename, "w") as outfile:
		outfile.write(json_object)

hugedata = set(["graph500-scale23-ef16_adj","bn-human-BNU_1_0025919_S1_remapped", "bio-mouse-gene_remapped.hdf5"])
withLabs = set(["fb-CMU-Carnegie49.hdf5","proteins-all.hdf5", "soc-YouTube-ASU.hdf5"])
DBLPNAME = "dblp-2024-08-01.hdf5"
newproc = set(["fb-CMU-Carnegie49.hdf5", "graph500-scale23-ef16_adj.hdf5", "proteins-all.hdf5", "soc-YouTube-ASU.hdf5"])
dataQ = set(["proteins-all.hdf5", "fb-CMU-Carnegie49.hdf5"])
parallelData = set(["graph500-scale23-ef16_adj", "bn-human-BNU_1_0025919_S1_remapped"])
MODEEX = 0
MODEQVAL = 1
MODETOHDF5 = 2
MODEEXCOUNTS = 3
MODEPARTITION = 4
MODEGETSTATS = 5
MODETOTXT = 10
MODEFIXEDSS = 13
MODEPROCESSDBLPTEMPORAL = 15
MODERUNWITHADAPTEPS = 16
MODEPARALLEL = 17

folderDataCsv = "data/raw/"
folderDataProc = "data/proc/"
folderResultsQ = "results/qval/"
folderResultsOursFixed = "results/ouralg/fixedss/"
folderLabels = "data/raw/labelled/"
folderPlaintTxt = "data/raw/txt/"
folderResultsAdapt = "results/adapt/"
folderResultsParallel = "results/parallel/"

folderEx = "results/exactStats/"

folderDataLargCsv = "data/raw/"
folderDataLargeProc = "data/proc/"

dblpCSV = "dblp-2024-08-01.raw.csv"
outdblphdf5 = "dblp-2024-08-01.hdf5"
dblplabs = "dblp-2024-08-01.raw-authorsMap.json"

dblpTemporal = "dblp-2024-08-01-SnapshotsYears.json"
dblpTemporalOutFile = "dblpTemporal.hdf5"
dblpTemporalOutJSON = "dblpTemporalRes.json"

CPPEXE = "./clustcoeff"
folderTmp = "tmp/"
mode = sys.argv[1]


if mode == "tohdf5":
	tasks = []
	remap = False
	missingdata = int(sys.argv[2]) # 0 if no, 1 if yes
	if(missingdata):
		folderDataCsv = folderDataLargCsv 
		folderDataProc = folderDataLargeProc 
		remap = True
	for i, fr in enumerate(os.listdir(folderDataCsv)):
		if(os.path.isdir(folderDataCsv+fr)):
			continue
		filenoext = fr.split(".")[0]
		huge = False
		if filenoext in hugedata:
			huge = True
		outfile = filenoext + ".hdf5"
		dict_conf = {
		"dataset":folderDataCsv+fr, 
		"outFile":folderDataProc+outfile, 
		"remap": remap,
		"huge": huge,
		"mode":MODETOHDF5
		}
		jsonfilename = folderTmp + "conf_tohdf5_" + str(i)+".json"
		dumpJson(dict_conf, jsonfilename)
		if os.path.isfile(folderDataProc+outfile):
			print("File already exists, skipping", folderDataProc+outfile)
			continue
		task = [CPPEXE, jsonfilename]
		tasks.append(task)

	if(len(tasks) > 0):
		with ThreadPoolExecutor(max_workers=6) as executor:
			executor.map(subprocess.run, tasks)


if mode == "dblpTemporal":
	dict_conf = {
		"fileData":dblpTemporal, 
		"outFile":dblpTemporalOutFile, 
		"resfile":dblpTemporalOutJSON, 
		"mode":MODEPROCESSDBLPTEMPORAL 
		}
	jsonfilename = folderTmp + "conf_temporalDBLP.json"
	dumpJson(dict_conf, jsonfilename)
	exit(0)
	task = [CPPEXE, jsonfilename]
	subprocess.run(task)

elif mode == "exactTriCounts":
	tasks = []
	for i, fr in enumerate(os.listdir(folderDataProc)):
		filenoext = fr.split(".")[0]
		huge = False
		if filenoext in hugedata:
			huge = True
		oldproc = True
		if fr in newproc:
			oldproc = False
		f5 = h5py.File(folderDataProc+fr, 'r')
		if 'exact-counts' in f5.keys():
			print("skipping file:", fr, " as already has counts")
			f5.close()
			continue
		print("examinig", fr)
		f5.close()
		dict_conf = {
		"dataset":folderDataProc+fr, 
		"oldread":oldproc,
		"huge":huge,
		"mode":MODEEXCOUNTS
		}
		jsonfilename = folderTmp + "conf_exCnts" + str(i)+".json"
		dumpJson(dict_conf, jsonfilename)
		task = [CPPEXE, jsonfilename]
		tasks.append(task)
	print(tasks)

	tot_w = min(6, len(tasks))
	if(len(tasks) > 0):
		with ThreadPoolExecutor(max_workers=tot_w) as executor:
			executor.map(subprocess.run, tasks)

elif mode == "modePartition":
	tasks = []
	partBucketsDegs = 25
	partCores = 30
	files = []
	for i, fr in enumerate(os.listdir(folderDataProc)):
		filenoext = fr.split(".")[0]
		readpart = h5py.File(folderDataProc+fr, 'r')
		if "vertexSetPartitions" in readpart.keys():
			print("Skipping file as already has partitions enable recompute if you want to rewrite", fr)
			readpart.close()
			continue
		readpart.close()
		outfile = filenoext + ".hdf5"
		oldproc = True
		hasLabs = False
		locLab = ""
		if fr in withLabs:
			hasLabs = True
			locLab = folderLabels + "LABS-" + filenoext + ".txt"
		if fr in newproc:
			oldproc = False
		huge = False
		if filenoext in hugedata:
			huge = True
		if not huge:
			if not os.path.isfile(folderPlaintTxt + filenoext + ".txt"):
				print("generating",folderPlaintTxt + filenoext + ".txt")
				dict_conf_txt = {
				"dataset":folderDataProc+fr, 
				"oldread":oldproc,
				"mode":MODETOTXT,
				"outfile":folderPlaintTxt + filenoext + ".txt"
				}
				jsonftxt = folderTmp + "conf_toTxt" + str(i) +".json"
				dumpJson(dict_conf_txt, jsonftxt)
				subprocess.run([CPPEXE, jsonftxt])
		dict_conf = {
		"dataset":folderDataProc+fr, 
		"mode":MODEPARTITION,
		"bucketsDegs": partBucketsDegs,
		"oldread":oldproc,
		"labsPath":locLab,
		"huge":huge,
		"fileTxt":folderPlaintTxt + filenoext + ".txt",
		"hasLabels":hasLabs,
		"partsCores": partCores
		}
		jsonfilename = folderTmp + "conf_partition" + str(i)+".json"
		dumpJson(dict_conf, jsonfilename)
		files.append(jsonfilename)
		task = [CPPEXE, jsonfilename]
		tasks.append(task)

	print("dumping normal partitions")
	if(len(tasks) > 0):
		with ThreadPoolExecutor(max_workers=6) as executor:
			executor.map(subprocess.run, tasks)

elif mode == "modeTestQ":
	tasks = []
	itersVect = [0.0001, 0.001, 0.01]
	qvalVect = [0.001, 0.01, 0.1, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
	for i, fr in enumerate(dataQ):
		filenoext = fr.split(".")[0]
		oldproc = True
		if fr in newproc:
			oldproc = False
		paths = []
		f5 = h5py.File(folderDataProc+fr, 'r')
		for path in f5['/vertexSetPartitions'].keys():
			currpath = "/vertexSetPartitions/" + path
			paths.append(currpath)
		f5.close()
		print(paths)
		outfile = filenoext + ".hdf5"
		dict_conf = {
		"dataset":folderDataProc+fr, 
		"outJson": folderResultsQ+filenoext +".json",
		"mode":MODEQVAL,
		"oldread":oldproc,
		"itersPercVect": itersVect,
		"qvalVect": qvalVect,
		"paths": paths
		}
		jsonfilename = folderTmp + "conf_testQ" + str(i)+".json"
		dumpJson(dict_conf, jsonfilename)
		task = [CPPEXE, jsonfilename]
		tasks.append(task)
		print(task)

	if(len(tasks) > 0):
		with ThreadPoolExecutor(max_workers=6) as executor:
			executor.map(subprocess.run, tasks)

elif mode == "modeFixedSampleSizeOurs":
	tasks = []
	worksmall = 30
	sampleSize = [1000, 500, 200]
	for samp in sampleSize:
		for i, fr in enumerate(os.listdir(folderDataProc)):
			filenoext = fr.split(".")[0]
			oldproc = True
			sampleSize = 100
			pathEdges = "/edges"
			if fr in newproc:
				oldproc = False
				pathEdges = "/edges/batch-0"
			paths = []
			f5 = h5py.File(folderDataProc+fr, 'r')
			M = f5[pathEdges].attrs['edges']
			for path in f5['/vertexSetPartitions'].keys():
				currpath = "/vertexSetPartitions/" + path
				paths.append(currpath)
			f5.close()
			huge = False
			if filenoext in hugedata:
				huge = True
			if huge:
				sampleSize = int(M/(1.*samp))
				worksmall = 100
			else:
				sampleSize = int(M/(1.*samp))
				worksmall = 30
			print(paths)
			outfile = filenoext + ".hdf5"
			if not os.path.isdir(folderResultsOursFixed + "/" + str(samp)):
				os.mkdir(folderResultsOursFixed + "/" + str(samp))
			dict_conf = {
			"dataset":folderDataProc+fr, 
			"outJson": folderResultsOursFixed+"/" + str(samp) + "/" +filenoext +".json",
			"mode":MODEFIXEDSS,
			"oldread":oldproc,
			"worksmall":worksmall,
			"samples":sampleSize,
			"paths": paths
			}
			jsonfilename = folderTmp + "conf_fixedSSOurs" + str(i)+".json"
			dumpJson(dict_conf, jsonfilename)
			task = [CPPEXE, jsonfilename]
			tasks.append(task)

		if(len(tasks) > 0):
			with ThreadPoolExecutor(max_workers=8) as executor:
				executor.map(subprocess.run, tasks)

if mode == "exactCountsParts":
	tasks = []
	for i, fr in enumerate(os.listdir(folderDataProc)):
		filenoext = fr.split(".")[0]
		oldproc = True
		if fr in newproc:
			oldproc = False
		paths = []
		f5 = h5py.File(folderDataProc+fr, 'r')
		for path in f5['/vertexSetPartitions'].keys():
			currpath = "/vertexSetPartitions/" + path
			paths.append(currpath)
		f5.close()
		outfile = filenoext + ".hdf5"
		dict_conf = {
		"dataset":folderDataProc+fr, 
		"outJson": folderEx+filenoext +".json",
		"oldread":oldproc,
		"mode":MODEGETSTATS,
		"paths": paths
		}
		jsonfilename = folderTmp + "conf_getStats" + str(i)+".json"
		dumpJson(dict_conf, jsonfilename)
		task = [CPPEXE, jsonfilename]
		tasks.append(task)

	if(len(tasks) > 0):
		with ThreadPoolExecutor(max_workers=6) as executor:
			executor.map(subprocess.run, tasks)

elif mode == "modeCompAdapt":
	tasks = []
	epsilonV = [0.075]
	for i, fr in enumerate(os.listdir(folderDataProc)):
		tl = 600
		worksmall = 150
		filenoext = fr.split(".")[0]
		oldproc = True
		sampleSize = 100
		pathEdges = "/edges"
		if fr in newproc:
			oldproc = False
			pathEdges = "/edges/batch-0"
		paths = []
		f5 = h5py.File(folderDataProc+fr, 'r')
		M = f5[pathEdges].attrs['edges']
		for path in f5['/vertexSetPartitions'].keys():
			currpath = "/vertexSetPartitions/" + path
			paths.append(currpath)
		f5.close()

		huge = False
		if filenoext in hugedata:
			huge = True
			worksmall = 500
			tl = 3600
		outfile = filenoext + ".hdf5"
		dict_conf = {
		"dataset":folderDataProc+fr, 
		"outJson": folderResultsAdapt + filenoext +".json",
		"worksmall":worksmall,
		"epsilon":epsilonV,
		"timeLimit":tl,
		"mode":MODERUNWITHADAPTEPS,
		"oldread":oldproc,
		"paths": paths
		}
		jsonfilename = folderTmp + "conf_compAdapt" + str(i)+".json"
		dumpJson(dict_conf, jsonfilename)
		task = [CPPEXE, jsonfilename]
		tasks.append(task)

	if(len(tasks) > 0):
		with ThreadPoolExecutor(max_workers=10) as executor:
			executor.map(subprocess.run, tasks)

elif mode == "modeParallel":
	tasks = []
	worksmall = 30
	samp= 500
	for i, fr in enumerate(os.listdir(folderDataProc)):
		filenoext = fr.split(".")[0]
		oldproc = True
		sampleSize = 100
		pathEdges = "/edges"
		if fr in newproc:
			oldproc = False
			pathEdges = "/edges/batch-0"
		paths = []
		f5 = h5py.File(folderDataProc+fr, 'r')
		M = f5[pathEdges].attrs['edges']
		for path in f5['/vertexSetPartitions'].keys():
			currpath = "/vertexSetPartitions/" + path
			paths.append(currpath)
		f5.close()
		huge = False
		if filenoext not in parallelData:
			continue
		if filenoext in hugedata:
			huge = True
		if huge:
			sampleSize = int(M/(1.*samp))
			worksmall = 100
		else:
			sampleSize = int(M/(1.*samp))
			worksmall = 10
		outfile = filenoext + ".hdf5"
		dict_conf = {
		"dataset":folderDataProc+fr, 
		"outJson": folderResultsParallel + filenoext +".json",
		"mode":MODEPARALLEL,
		"oldread":oldproc,
		"worksmall":worksmall,
		"samples":sampleSize,
		"threads":5,
		"paths": paths
		}
		jsonfilename = folderTmp + "conf_parallel" + str(i)+".json"
		dumpJson(dict_conf, jsonfilename)
		task = [CPPEXE, jsonfilename]
		tasks.append(task)

	if(len(tasks) > 0):
		with ThreadPoolExecutor(max_workers=2) as executor:
			executor.map(subprocess.run, tasks)

