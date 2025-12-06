from typing import List, Tuple

import matplotlib.pyplot as plt
import seaborn as sns

# Gaurav Sablok
# codeprog@icloud.com

def mrnasequencecompare(pathfile1: str, pathfile2:str, gff_file1: str, gff_file2:str) -> str:
	"""
	 This analyses the sequence as well as the mrna comparison.
	"""
	id1:List[str] = []
	id2:List[str] = []
	proteincoordinates1:List[Tuple[int, int]] = []
	proteincoordinates2: List[Tuple[int, int]] = []
	strand1:List[str] = []
	strand2:List[str] = []
	with open(gff_file1, 'r') as pathopenfile1:
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
	with open(gff_file2, 'r') as pathopenfile2:
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
	protein_1_names:List[str] = []
	protein_2_names:List[str] = []
	protein_1_seq:List[str] = []
	protein_2_seq:List[str] = []
	with open(pathfile1, 'r') as a:
		for i in a.readlines():
			if i.startswith(">"):
				protein_1_names.append(i.strip().split("|")[0].split("-")[0].replace(">", ""))
			if not i.startswith(">"):
				protein_1_seq.append(i.strip())
	with open(pathfile2, 'r') as b:
		for i in b.readlines():
			if i.startswith(">"):
				protein_2_names.append(i.strip().split("|")[0].split(".")[0].replace(">", ""))
			if not i.startswith(">"):
				protein_2_seq.append(i.strip())

	sameid_same_seq:List[Tuple[str, str]] = []
	for i in range(len(protein_1_names)):
		for val in range(len(protein_2_names)):
			if protein_1_names[i] == protein_2_names[val] and protein_1_seq[i] == protein_2_seq[val]:
				sameid_same_seq.append((protein_1_names[i],protein_2_names[val]))
	with open("sameid_same_seq", 'w') as filewrite:
		filewrite.write("id1" + '\t' + "id2" + '\n')
		for i in range(len(sameid_same_seq)):
			filewrite.write(str(sameid_same_seq[i][0]) + '\t' + str(sameid_same_seq[i][1] + '\n'))
	sameid_different_seq:List[Tuple[str,str]] = []
	for i in range(len(protein_1_names)):
		for val in range(len(protein_2_names)):
			if protein_1_names[i] == protein_2_names[val] and protein_1_seq[i] != protein_2_seq[val]:
				sameid_different_seq.append((protein_1_names[i],protein_2_names[val]))
	with open("sameid_different_seq", 'w') as filewrite1:
		filewrite1.write("id1" + '\t' + "id2" + '\n')
		for i in range(len(sameid_different_seq)):
			filewrite1.write(str(sameid_different_seq[i][0]) + '\t' + str(sameid_different_seq[i][1] + '\n'))
	sameidslength: int = len(sameid_same_seq)
	sameidsdifferentseqlength:int = len(sameid_different_seq)
	differentids = [i for i in protein_1_names if i not in protein_2_names]
	differentidslength = len(differentids)
	plotting_variable = [str(sameidslength), str(sameidsdifferentseqlength), str(differentidslength)]
	_ = sns.histplot(x = plotting_variable, y = ["sameids-length", "sameids-different-seqlength", "different-idslength"])
	plt.savefig("histogram of the difference.png", dpi = 300)
	return f"The number of the different ids are: {differentids}"
