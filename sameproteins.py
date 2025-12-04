
import matplotlib.pyplot as plt
import seaborn as sns
from typing_extensions import List, Tuple

# Gaurav Sablok
# codeprog@icloud.com

def toxcompare_same_proteins(pathfile1: str, pathfile2:str):
	"""
	   This function takes the path of the ToxannotationDB file
				and compared for the sequences of the protein coding.
	   This compares the annotation records of the two gff files and
				show the similar geneid and their associated information.
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
				if line[2] == "protein_coding_gene":
					start = int(line[3])
					stop = int(line[4])
					strand = str(line[6])
					annotation = str(line[8].split(";")[0]).replace("ID=", "")
					id1.append(annotation)
					strand1.append(strand)
					proteincoordinates1.append((start,stop))
	with open(pathfile2, 'r') as pathopenfile2:
		for val in pathopenfile2.readlines():
			if not val.startswith("#"):
				line = val.strip().split("\t")
				if line[2] == "protein_coding_gene":
					start = int(line[3])
					stop = int(line[4])
					strand = str(line[6])
					annotation = str(line[8].split(";")[0]).replace("ID=", "")
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
				strand1insert = strand1[i]
				id2insert = id2[val]
				id2start = int(proteincoordinates2[val][0])
				id2end = int(proteincoordinates2[val][1])
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
		proteinlength1.append(abs(int(common_proteins[i][0][2]- common_proteins[i][0][1])))
	for i in range(len(common_proteins)):
		proteinnames2.append(common_proteins[i][1][0])
		proteinlength2.append(abs(int(common_proteins[i][1][2]- common_proteins[i][1][1])))
	proteindifference:List[Tuple[str, int,int]] = []
	for i in range(len(common_proteins)):
		valuename:str = common_proteins[i][0][0]
		startdifference: int = common_proteins[i][0][1]-common_proteins[i][1][1]
		enddifference:int = common_proteins[i][0][2] - common_proteins[i][1][2]
		proteindifference.append((valuename, startdifference, enddifference))
	proteindifferencenames:list[str] = [proteindifference[i][0] for i in range(len(proteindifference)) ]
	proteinstartdifflength:list[int] = [proteindifference[i][1] for i in range(len(proteindifference)) ]
	proteinenddifferencelength:List[int] = [proteindifference[i][2] for i in range(len(proteindifference))]
	_ = sns.barplot(x = proteindifferencenames, y = proteinstartdifflength, dpi= 300)
	plt.savefig("startdifference.png")
	_ = sns.barplot(x = proteindifferencenames, y = proteinenddifferencelength, dpi = 300)
	plt.savefig("enddifference.png")
	_ = sns.barplot(x = proteinnames1, y = proteinlength1, hue = proteinlength1)
	plt.savefig('barplot-1.png', dpi=300)
	_ = sns.barplot(x = proteinnames2, y = proteinlength2, hue = proteinlength2)
	plt.savefig('barplot-2.png', dpi=300)
	with open("comparative-common.txt", 'w') as filewrite:
		filewrite.write("id1"+ '\t' + "id2" + '\t' + "id1start" + '\t' + "id2start" + '\t' + "id1end" + '\t' + "id2end" + '\t' + "id1strand1" + '\t' + "idstrand2" + "\n")
		for i in range(len(common_proteins)):
			filewrite.write(str(common_proteins[i][0][0]) +'\t' + str(common_proteins[i][0][1]) +'\t' + str(common_proteins[i][0][2]) +'\t' + str(common_proteins[i][0][3]) +'\t' + str(common_proteins[i][1][0]) +'\t'
				+ str(common_proteins[i][1][1]) +'\t' + str(common_proteins[i][1][2]) +'\t' + str(common_proteins[i][1][3]) + '\n')
		filewrite.close()
	with open("difference-common-protein.txt", 'w') as filecommon:
		filecommon.write("id" + '\t' + "start difference" + "end difference" + '\n')
		for i in range(len(proteindifference)):
				filecommon.write(str(proteindifference[i][0]) + '\t' + str(proteindifference[i][1]) + str(proteindifference[i][2]) + '\n')
		filecommon.close()
