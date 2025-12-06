from typing_extensions import List, Tuple

# Gaurav Sablok
# codeprog@icloud.com

def mrna_comparison(pathfile1:str, pathfile2:str) -> str:
	"""
	This analyzes the mrna from the gff for both the gff.
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
				if line[2] == "mRNA":
					start = int(line[3])
					stop = int(line[4])
					strand = str(line[6])
					annotation = str(line[8].split(";")[1]).replace("Parent=", "")
					id1.append(annotation)
					strand1.append(strand)
					proteincoordinates1.append((start,stop))
	with open(pathfile2, 'r') as pathopenfile2:
		for val in pathopenfile2.readlines():
			if not val.startswith("#"):
				line = val.strip().split("\t")
				if line[2] == "mRNA":
					start = int(line[3])
					stop = int(line[4])
					strand = str(line[6])
					annotation = str(line[8].split(";")[1]).replace("Parent=", "")
					strand2.append(strand)
					id2.append(annotation)
					proteincoordinates2.append((start,stop))
	common_mrnas:List[List[Tuple]] = []
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
				common_mrnas.append([insertvalue1, insertvalue2])
# different start and same end for the binning.
	different_start_same_end = []
	for i in range(len(common_mrnas)):
		if common_mrnas[i][0][0] == common_mrnas[i][1][0] and common_mrnas[i][0][1] \
		!= common_mrnas[i][1][1] and common_mrnas[i][0][2] == common_mrnas[i][1][2]:
				different_start_same_end.append((common_mrnas[i][0][0], common_mrnas[i][0][2]-common_mrnas[i][1][2]))
	with open("different_start_same_end.txt", 'w') as filewrite:
			filewrite.write("mRNA" + '\t' + "shiftbp" + '\n')
			for i in range(len(different_start_same_end)):
				filewrite.write(str(different_start_same_end[i][0]) + '\t' + str(different_start_same_end[i][1]) + '\n')
			filewrite.close()
# same start and different end for the binning.
	same_start_different_end = []
	for i in range(len(common_mrnas)):
		if common_mrnas[i][0][0] == common_mrnas[i][1][0] and common_mrnas[i][0][1] \
		== common_mrnas[i][1][1] and common_mrnas[i][0][1] != common_mrnas[i][1][1]:
			same_start_different_end.append((common_mrnas[i][0][0], common_mrnas[i][0][2]-common_mrnas[i][1][2]))
	with open("same_start_different_end.txt", 'w') as filewrite:
			filewrite.write("mRNA" + '\t' + "shiftbp" + '\n')
			for i in range(len(same_start_different_end)):
					filewrite.write(str(same_start_different_end[i][0]) + '\t' + str(same_start_different_end[i][1]) + '\n')
			filewrite.close()
# ids same but different start and different end
	sameids_different_start_end = []
	for i in range(len(common_mrnas)):
		if common_mrnas[i][0][0] == common_mrnas[i][1][0] and common_mrnas[i][0][1] \
		!= common_mrnas[i][1][1] and common_mrnas[i][0][2] != common_mrnas[i][1][2]:
			sameids_different_start_end.append((common_mrnas[i][0][0], common_mrnas[i][0][2]-common_mrnas[i][1][2]))
	with open("sameids_different_start.txt", 'w') as filewrite:
			filewrite.write("mRNA" + '\t' + "shiftbp" + '\n')
			for i in range(len(sameids_different_start_end)):
					filewrite.write(str(sameids_different_start_end[i][0]) + '\t' + str(sameids_different_start_end[i][1]) + '\n')
			filewrite.close()
	return "The files have been written for the comparative analysis"
