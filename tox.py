#! usr/bin/env/py

# Gaurav Sablok
# codeprog@icloud.com
# Toxpy: python version for the comparative genomics of the ToxDB species
# It provided you the complete machine leaning as well as comparative genomics.


from typing import Tuple

import seaborn as sns
from typing_extensions import List


def toxannotator(pathfile1, pathfile2):
	"""
	   This function takes the path of the Toxannotation file
				and compared for the sequences of the protein coding.
	"""
	id1:List = []
	id2:List = []
	proteincoordinates1:List[Tuple] = []
	proteincoordinates2: List[Tuple] = []
	strand1:List = []
	strand2:List = []
	with open(pathfile1, 'r') as pathopenfile1:
		for i in pathopenfile1.readlines():
			if not i.startswith("#"):
				line = i.strip().split("\t")
				start = int(line[3])
				stop = int(line[4])
				strand = str(line[6])
				annotation = str(line[8].split("\t")[0]).replace("ID=", "")
				id1.append(annotation)
				strand1.append(strand)
				proteincoordinates1.append((start,stop))
	with open(pathfile2, 'r') as pathopenfile2:
		for val in pathopenfile2.readlines():
			if not val.startswith("#"):
				line = val.strip().split("\t")
				start = int(line[3])
				stop = int(line[4])
				strand = str(line[6])
				annotation = str(line[8].split("\t")[0]).replace("ID=", "")
				strand2.append(strand)
				id2.append(annotation)
				proteincoordinates2.append((start,stop))
	common_proteins:List[List[Tuple]] = []
	for i in range(len(id1)):
		for val in range(len(id2)):
			if id1[i] == id2[val]:
				id1insert = id1[i]
				id1start = int(proteincoordinates1[i][0])
				id1end = int(proteincoordinates1[i][1])
				strand1 = strand1[i]
				id2insert = id2[val]
				id2start = int(proteincoordinates2[i][0])
				id2end = int(proteincoordinates2[i][1])
				strand2 = strand2[val]
				insertvalue1 = (id1insert, id1start, id1end, strand1)
				insertvalue2 = (id2insert, id2start, id2end, strand2)
				common_proteins.append([insertvalue1, insertvalue2])
	noncommon_proteins:List[list[Tuple]] = []
	for i in range(len(id1)):
		for val in range(len(id2)):
			if id1[i] != id2[val]:
				id1insert = id1[i]
				id1start = int(proteincoordinates1[i][0])
				id1end = int(proteincoordinates1[i][1])
				strand1 = strand1[i]
				id2insert = id2[val]
				id2start = int(proteincoordinates2[i][0])
				id2end = int(proteincoordinates2[i][1])
				strand2 = strand2[val]
				insertvalue1 = (id1insert, id1start, id1end, strand1)
				insertvalue2 = (id2insert, id2start, id2end, strand2)
				noncommon_proteins.append([insertvalue1, insertvalue2])
	proteinnames1:List = []
	proteinlength1:List = []
	proteinnames2:List= []
	proteinlength2:List = []
	for i in range(len(id1)):
		proteinnames1.append(id1[i])
		proteinlength1.append(int(proteincoordinates1[i][1]- proteincoordinates1[i][0]))
	for i in range(len(id2)):
		proteinnames2.append(id2[i])
		proteinlength2.append(int(proteincoordinates2[i][1]- proteincoordinates2[i][0]))
	barplot1 = sns.barplot(x = proteinnames1, y = proteinlength1, hue = proteinlength1)
	barplot2 = sns.barplot(x = proteinnames2, y = proteinlength2, hue = proteinlength2)
