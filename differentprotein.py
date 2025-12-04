from typing_extensions import List, Tuple

# Gaurav Sablok
# codeprog@icloud.com

def toxcompare_different_proteins(pathfile1:str, pathfile2:str):
	"""
	   This function takes the path of the ToxannotationDB file
				and compared for the sequences of the protein coding.
	   This compares the two annotation files and shows the missing
				proteins which are not present in them and their coordinates.

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
	noncommon_proteins_1:List[list[Tuple]] = []
	for i in range(len(id1)):
		for val in range(len(id2)):
			if id1[i] != id2[val]:
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
				noncommon_proteins_1.append([insertvalue1, insertvalue2])
	with open("comparative-common.txt", 'w') as filewrite_1:
		filewrite_1.write("The comparative table for the different genes and protein")
		filewrite_1.write("id1"+ '\t' + "id1start" + '\t' + "id1end" + '\t' + "id1strand" + '\t' + "id2" + '\t' + "id2start" + '\t' + "id2end" + '\t' + "id2stand" + '\n')
		for i in range(len(noncommon_proteins_1)):
			filewrite_1.write(str(noncommon_proteins_1[i][0][0]) +'\t' + str(noncommon_proteins_1[i][0][1]) +'\t' + str(noncommon_proteins_1[i][0][2]) +'\t' + str(noncommon_proteins_1[i][0][3]) +'\t' + str(noncommon_proteins_1[i][1][0]) +'\t' + str(noncommon_proteins_1[i][1][1]) +'\t' + str(noncommon_proteins_1[i][1][2])
				+'\t' + str(noncommon_proteins_1[i][1][3]) +'\n')
	filewrite_1.close()
