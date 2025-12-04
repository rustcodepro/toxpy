import matplotlib.pyplot as plt
import seaborn as sns
from typing_extensions import List, Tuple

# Gaurav Sablok
# codeprog@icloud.com

def toxcompare_same_proteins_different_strand(pathfile1:str, pathfile2:str):
	"""
	   This function takes the path of the ToxannotationDB file
				and compared for the sequences of the protein coding.
	   This compares the annotation records of the two gff files and
				show the similar geneid and different strands their associated information.
	"""
	id1:List[str] = []
	id2:List[str] = []
	proteincoordinates1:List[Tuple[int, int]] = []
	proteincoordinates2: List[Tuple[int, int]] = []
	strand1:List[str] = []
	strand2:List[str] = []
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
			if id1[i] == id2[val] and strand1[i] != strand2[val]:
				id1insert = id1[i]
				id1start = int(proteincoordinates1[i][0])
				id1end = int(proteincoordinates1[i][1])
				strand1insert = strand1[i]
				id2insert = id2[val]
				id2start = int(proteincoordinates2[i][0])
				id2end = int(proteincoordinates2[i][1])
				strand2insert = strand2[val]
				insertvalue1 = (id1insert, id1start, id1end, strand1insert)
				insertvalue2 = (id2insert, id2start, id2end, strand2insert)
				common_proteins.append([insertvalue1, insertvalue2])
	proteinnames1:List[str] = []
	proteinlength1:List[int] = []
	proteinnames2:List[str]= []
	proteinlength2:List[int] = []
	for i in range(len(common_proteins)):
		proteinnames1.append(common_proteins[i][0][0])
		proteinlength1.append(int(common_proteins[i][0][2]- common_proteins[i][0][1]))
	for i in range(len(common_proteins)):
		proteinnames1.append(common_proteins[i][1][0])
		proteinlength1.append(int(common_proteins[i][1][2]- common_proteins[i][1][1]))
	barplot1 = sns.barplot(x = proteinnames1, y = proteinlength1, hue = proteinlength1)
	plt.savefig('barplot-1.png', dpi=300)
	barplot2 = sns.barplot(x = proteinnames2, y = proteinlength2, hue = proteinlength2)
	plt.savefig('barplot-2.png', dpi=300)
	with open("comparative-common.txt", 'w') as filewrite:
		filewrite.write("The comparative table for the different genes and protein")
		filewrite.write("id1"+ '\t' + "id1start" + '\t' + "id1end" + '\t' + "id1strand" + '\t' + "id2" + '\t' + "id2start" + '\t' + "id2end" + '\t' + "id2stand")
		for i in range(len(common_proteins)):
			filewrite.write(common_proteins[i][0][0] +'\t' + common_proteins[i][0][1] +'\t' + common_proteins[i][0][2] +'\t' + common_proteins[i][0][3] +'\t' + common_proteins[i][1][0] +'\t' + common_proteins[i][1][1] +'\t' +common_proteins[i][1][2]
				+'\t' + common_proteins[i][1][3])
		filewrite.close()
