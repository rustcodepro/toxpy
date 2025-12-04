#! usr/bin/env/py

# Gaurav Sablok
# codeprog@icloud.com

# Toxpy: python version for the comparative genomics of the ToxDB species
# It provided you the complete machine leaning as well as comparative genomics.
# # a complete package which performs the ToxDB comparison and also fits a LSTM
# and Convolutional Neural Network with L2 regularization
# why i implemented L2 regulrization as PTM have the end terminal
# modification and hence they will reach the vanishing gradient earlier.


import random

import click
import numpy as np
import seaborn as sns
from click.decorators import command, option
from pandas.core.generic import common
from sklearn.preprocessing import LabelEncoder
from tensorflow.keras.layers import LSTM, Dense, Embedding
from tensorflow.keras.models import Sequential
from tensorflow.keras.preprocessing.sequence import pad_sequences
from typing_extensions import Dict, List, Tuple


@click.group()
@click.command()
@click.option("-dnaseq", "-classificationfile", "-ptmsequences",
	default = "lstm training on the PTM proteins",
	help = "pass the dna sequences and the classification categories as a fasta file with the cateories in the header")
# Toxdb LSTM model
def tox_lstm(dnaseq,classificationfile, ptmsequences):
	"""
	  This function fits the LSTM model for the classification of the
			PTM(post translational modifications) using the Keras and Tensorflow
			LSTM layers. Optimized for the PTM.
	"""
	seqdict =  {}
	header:List[str] = []
	sequence:List[str] = []
	with open(dnaseq, 'r') as pathopenfile:
		for i in pathopenfile.readlines():
			if i.startswith(">"):
				header.append(i.strip().replace(">", ""))
			elif not i.startswith(">"):
				sequence.append(i.strip())
			else:
				pass
	for i in range(len(header)):
		seqdict[header[i]] = sequence[i]
	classification_id:List[int] = []
	with open(classificationfile, 'r') as clasificationids:
		for i in clasificationids.readlines():
			line = i.strip()
			classification_id.append(int(line))
	header:List[str] = []
	sequence:List[str] = []
	seqdictptm:Dict =  {}
	with open(ptmsequences, 'r') as pathopenfile:
		for i in pathopenfile.readlines():
			if i.startswith(">"):
				header.append(i.strip().replace(">", ""))
			elif not i.startswith(">"):
				sequence.append(i.strip())
			else:
				pass
	for i in range(len(header)):
		seqdictptm[header[i]] = sequence[i]
	def encode_dna_sequence(seq):
		mapping = {'A': 1, 'C': 2, 'G': 3, 'T': 4}
		return [mapping[base] for base in seq]
	finalsequences:List[str] = [seq for seq in seqdict.values()]
	encoded_sequences = [encode_dna_sequence(seq) for seq in finalsequences]
	max_length = max(len(seq) for seq in finalsequences)
	X = pad_sequences(encoded_sequences, maxlen=max_length, padding='post')
	y = np.array(classification_id)
	ptmsequence:List[str] = [sequences for sequences in ptmsequences.values()]
	model = Sequential([
    	Embedding(input_dim=5, output_dim=8, input_length=max_length),
     LSTM(16, return_sequences=False),
     Dense(8, activation='relu'),
     Dense(1, activation='sigmoid')
])
	model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
	model.fit(X, y, epochs=10, batch_size=2, verbose=1)
	testsequence = [encode_dna_sequence(seq) for seq in ptmsequence]
	max_length = max(len(seq) for seq in testsequence)
	X = pad_sequences(encoded_sequences, maxlen=max_length, padding='post')
	sequence = pad_sequences([testsequence], maxlen=max_length, padding='post')
	prediction = model.predict(sequence)
	print(f"Prediction for sequence: {'Positive' if prediction[0] > 0.5 else 'Negative'}")

#CNN model with L2 regularization
# Gaurav Sablok
# codeprog@icloud.com
@click.group()
@click.command()
@click.option("-dnaseq", "-classificationfile", "-ptmsequences",
	default = "lstm training on the PTM proteins",
	help = "pass the dna sequences and the classification categories as a fasta file with the cateories in the header")
def cnn_l2_regularization(dnaseq, classificationfile, ptmsequences):
	"""
	 This model fits the convolutional neural network(CNN) with L2 regularization to minimize the vanishing gradient
		and optimized for the Post Translational Modifications (PTM)
	"""
	seqdict:Dict[str, str] =  {}
	header:List[str] = []
	sequence:List[str] = []
	with open(dnaseq, 'r') as pathopenfile:
		for i in pathopenfile.readlines():
			if i.startswith(">"):
				header.append(i.strip().replace(">", ""))
			elif not i.startswith(">"):
				sequence.append(i.strip())
			else:
				pass
	for i in range(len(header)):
		seqdict[header[i]] = sequence[i]
	classification_id:List[int] = []
	with open(classificationfile, 'r') as clasificationids:
		for i in clasificationids.readlines():
			line = i.strip()
			classification_id.append(int(line))
	header:List[str] = []
	sequence:List[str] = []
	seqdictptm:Dict =  {}
	with open(ptmsequences, 'r') as pathopenfile:
		for i in pathopenfile.readlines():
			if i.startswith(">"):
				header.append(i.strip().replace(">", ""))
			elif not i.startswith(">"):
				sequence.append(i.strip())
			else:
				pass
	for i in range(len(header)):
		seqdictptm[header[i]] = sequence[i]
	def encode_dna_sequence(seq):
		mapping = {'A': 1, 'C': 2, 'G': 3, 'T': 4}
		return [mapping[base] for base in seq]
	finalsequences:List[str] = [seq for seq in seqdict.values()]
	encoded_sequences = [encode_dna_sequence(seq) for seq in finalsequences]
	max_length = max(len(seq) for seq in finalsequences)
	X = pad_sequences(encoded_sequences, maxlen=max_length, padding='post')
	y = np.array(classification_id)
	ptmsequence:List[str] = [sequences for sequences in ptmsequences.values()]
	X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
	model = Sequential()
	model.add(Conv1D(filters=32, kernel_size=8, activation='relu', input_shape=(X.shape[1], X.shape[2]),
                 kernel_regularizer=l2(0.01)))
	model.add(MaxPooling1D(pool_size=2))
	model.add(Conv1D(filters=64, kernel_size=8, activation='relu', kernel_regularizer=l2(0.01)))
	model.add(MaxPooling1D(pool_size=2))
	model.add(Flatten())
	model.add(Dense(64, activation='relu', kernel_regularizer=l2(0.01)))
	model.add(Dropout(0.2))
	model.add(Dense(1, activation='sigmoid'))
	model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
	early_stop = EarlyStopping(monitor='val_loss', patience=5, restore_best_weights=True)
	history = model.fit(X_train, y_train, epochs=50, batch_size=32, validation_split=0.2,
                    callbacks=[early_stop], verbose=1)
	y_pred = (model.predict(ptmsequence) > 0.5).astype(int)
	print(classification_report(y_test, y_pred))


@click.group()
@click.command()
@click.option("-pathfile1", "-pathfile2", default="pass a file", help=" pass the gff1 and gff2 for comparison")
def toxcompare_same_proteins(pathfile1, pathfile2):
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
	barplot2 = sns.barplot(x = proteinnames2, y = proteinlength2, hue = proteinlength2)
	barplot1.get_figure().savefig("barplot-1.png")
	barplot2.get_figure().savefig("barplot-2.png")
	with open("comparative-common.txt", 'w') as filewrite:
		filewrite.write("id1start"+ '\t' + "id2start" + '\t' + "id1start" + '\t' + "id2start" + '\t' + "id1end" + '\t' + "id2end" + '\t' + "id1strand1" + '\t' + "idstrand2")
		for i in range(len(common_proteins)):
			filewrite.write(common_proteins[i][0][0] +'\t' + common_proteins[i][0][1] +'\t' + common_proteins[i][0][2] +'\t' + common_proteins[i][0][3] +'\t' + common_proteins[i][1][0] +'\t' + common_proteins[i][1][1] +'\t' +common_proteins[i][1][2]
				+'\t' + common_proteins[i][1][3])
		filewrite.close()




@click.group()
@click.command()
@click.option("- pathfile1",
	"-pathfile2",
	default="pass a file",
	help = "This is for the comparison of the gff and writing down the comparative analysis")
def toxcompare_different_proteins(pathfile1, pathfile2):
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
	noncommon_proteins_1:List[list[Tuple]] = []
	for i in range(len(id1)):
		for val in range(len(id2)):
			if id1[i] != id2[val]:
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
				noncommon_proteins_1.append([insertvalue1, insertvalue2])
	with open("comparative-common.txt", 'w') as filewrite_1:
		filewrite_1.write("The comparative table for the different genes and protein")
		filewrite_1.write("id1"+ '\t' + "id1start" + '\t' + "id1end" + '\t' + "id1strand" + '\t' + "id2" + '\t' + "id2start" + '\t' + "id2end" + '\t' + "id2stand")
		for i in range(len(noncommon_proteins_1)):
			filewrite_1.write(noncommon_proteins_1[i][0][0] +'\t' + noncommon_proteins_1[i][0][1] +'\t' + noncommon_proteins_1[i][0][2] +'\t' + noncommon_proteins_1[i][0][3] +'\t' + noncommon_proteins_1[i][1][0] +'\t' + noncommon_proteins_1[i][1][1] +'\t' + noncommon_proteins_1[i][1][2]
							+'\t' + noncommon_proteins_1[i][1][3])
	filewrite_1.close()

@click.group()
@click.command()
@click.option("-pathfile1", "-pathfile2",
	default = "file for the comparison",
	help = "intersect the filesboth for the id and the sequences")
def toxcompare_same_proteins_with_same_strand(pathfile1, pathfile2):
	"""
	   This function takes the path of the ToxannotationDB file
				and compared for the sequences of the protein coding.
	   This compares the annotation records of the two gff files and
				show the similar geneid and similar strand and their associated information.
	"""
	id1:List[str] = []
	id2:List[str] = []
	proteincoordinates1:List[Tuple[int, int]] = []
	proteincoordinates2: List[Tuple[int,int]] = []
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
			if id1[i] == id2[val] and strand1[i] == strand2[val]:
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
	barplot2 = sns.barplot(x = proteinnames2, y = proteinlength2, hue = proteinlength2)
	barplot1.get_figure().savefig("barplot-1.png")
	barplot2.get_figure().savefig("barplot-2.png")
	with open("comparative-common.txt", 'w') as filewrite:
		filewrite.write("The comparative table for the different genes and protein")
		filewrite.write("id1"+ '\t' + "id1start" + '\t' + "id1end" + '\t' + "id1strand" + '\t' + "id2" + '\t' + "id2start" + '\t' + "id2end" + '\t' + "id2stand")
		for i in range(len(common_proteins)):
			filewrite.write(common_proteins[i][0][0] +'\t' + common_proteins[i][0][1] +'\t' + common_proteins[i][0][2] +'\t' + common_proteins[i][0][3] +'\t' + common_proteins[i][1][0] +'\t' + common_proteins[i][1][1] +'\t' +common_proteins[i][1][2]
				+'\t' + common_proteins[i][1][3])
		filewrite.close()


@click.group()
@click.command()
@click.option("-pathfile1", "-pathfile2",
	default = "file for the comparison",
	help = "intersect the filesboth for the id and the sequences")
def toxcompare_same_proteins_different_strand(pathfile1, pathfile2):
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
	barplot2 = sns.barplot(x = proteinnames2, y = proteinlength2, hue = proteinlength2)
	barplot1.get_figure().savefig("barplot-1.png")
	barplot2.get_figure().savefig("barplot-2.png")
	with open("comparative-common.txt", 'w') as filewrite:
		filewrite.write("The comparative table for the different genes and protein")
		filewrite.write("id1"+ '\t' + "id1start" + '\t' + "id1end" + '\t' + "id1strand" + '\t' + "id2" + '\t' + "id2start" + '\t' + "id2end" + '\t' + "id2stand")
		for i in range(len(common_proteins)):
			filewrite.write(common_proteins[i][0][0] +'\t' + common_proteins[i][0][1] +'\t' + common_proteins[i][0][2] +'\t' + common_proteins[i][0][3] +'\t' + common_proteins[i][1][0] +'\t' + common_proteins[i][1][1] +'\t' +common_proteins[i][1][2]
				+'\t' + common_proteins[i][1][3])
		filewrite.close()


@click.group
@click.command
@click.option("-gff_file1", "-gff_file2", "-fastafile1", "-fastafile2",
	default= "no default", help = "compare the gff based on the ids and the sequences")
def toxcompare_different_protein_same_sequences(gff_file1, gff_file2, fastafile1, fastafile2):
	"""
	   This function takes the path of the ToxannotationDB file
				and compared for the sequences of the protein coding.
	   This compares the annotation records of the two gff files and
				show the similar geneid and their associated information.
	   This is to check the annotation errors if by chance two ids have
				the same sequence.
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
				start = int(line[3])
				stop = int(line[4])
				strand = str(line[6])
				annotation = str(line[8].split("\t")[0]).replace("ID=", "")
				id1.append(annotation)
				strand1.append(strand)
				proteincoordinates1.append((start,stop))
	with open(gff_file2, 'r') as pathopenfile2:
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
	fastafile1_unpack = list(filter(None,[x.strip() for x in open(fastafile1).readlines()]))
	fastafile2_unpack = list(filter(None,[x.strip() for x in open(fastafile2).readlines()]))
	protein_1:Dict[str, str] = {}
	for i in fastafile1_unpack:
		if i.startswith(">"):
			transcript_path = i.strip().split("|")[2].replace("gene=", "")
			if i not in protein_1:
				protein_1[i] = ""
				continue
			protein_1[transcript_path] += i.strip()
	protein_2:Dict[str,str] = {}
	for i in fastafile2_unpack:
		if i.startswith(">"):
			transcript = i.strip().split("|")[2].replace("gene=", "")
			if i not in protein_2:
				protein_2[i] = ""
				continue
			protein_2[transcript] += i.strip()
	common_proteins:List[List[Tuple]] = []
	for i in range(len(id1)):
		for valrange in range(len(id2)):
			for (val, key) in protein_1.items():
				for (val2, key2) in protein_2.items():
					if id1[i] != id2[valrange] and  val != val2 and key == key2:
						id1insert = id1[i]
						id1start = int(proteincoordinates1[i][0])
						id1end = int(proteincoordinates1[i][1])
						strand1insert = strand1[i]
						id2insert = id2[valrange]
						id2start = int(proteincoordinates2[i][0])
						id2end = int(proteincoordinates2[i][1])
						strand2insert = strand2[valrange]
						sequence_1 = protein_1[val]
						sequence_2 = protein_2[val2]
						insertvalue1 = (id1insert, id1start, id1end, strand1insert, sequence_1)
						insertvalue2 = (id2insert, id2start, id2end, strand2insert, sequence_2)
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
	barplot2 = sns.barplot(x = proteinnames2, y = proteinlength2, hue = proteinlength2)
	barplot1.get_figure().savefig("barplot-1.png")
	barplot2.get_figure().savefig("barplot-2.png")
	with open("comparative-common.txt", 'w') as filewrite:
		filewrite.write("The comparative table for the different genes and protein")
		filewrite.write("id1"+ '\t' + "id1start" + '\t' + "id1end" + '\t' + "id1strand" + '\t' "id1sequence" + '\t' + "id2" + '\t' + "id2start" + '\t' + "id2end" + '\t' + "id2stand" + '\t' + "id2sequence")
		for i in range(len(common_proteins)):
			filewrite.write(common_proteins[i][0][0] +'\t' + common_proteins[i][0][1] +'\t' + common_proteins[i][0][2] +'\t' + common_proteins[i][0][3] +'\t' + common_proteins[i][1][4] + '\t' + common_proteins[i][1][0] +'\t' + common_proteins[i][1][1] +'\t' +common_proteins[i][1][2]
				+'\t' + common_proteins[i][1][3] + '\t' + common_proteins[i][1][4])
		filewrite.close()


@click.group
@click.command
@click.option("-gff_file1", "-gff_file2", "-fastafile1", "-fastafile2",
	default= "no default", help = "compare the gff based on the ids and the sequences")
def toxcompare_same_proteins_different_sequences(gff_file1, gff_file2, fastafile1, fastafile2):
	"""
	   This function takes the path of the ToxannotationDB file
				and compared for the sequences of the protein coding.
	   This compares the annotation records of the two gff files and
				show the similar geneid and their associated information.
	   This compares the gff and the fasta to have the same ids and the same sequences.
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
				start = int(line[3])
				stop = int(line[4])
				strand = str(line[6])
				annotation = str(line[8].split("\t")[0]).replace("ID=", "")
				id1.append(annotation)
				strand1.append(strand)
				proteincoordinates1.append((start,stop))
	with open(gff_file2, 'r') as pathopenfile2:
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
	fastafile1_unpack = list(filter(None,[x.strip() for x in open(fastafile1).readlines()]))
	fastafile2_unpack = list(filter(None,[x.strip() for x in open(fastafile2).readlines()]))
	protein_1 = {}
	for i in fastafile1_unpack:
		if i.startswith(">"):
			transcript_path = i.strip().split("|")[2].replace("gene=", "")
			if i not in protein_1:
				protein_1[i] = ""
				continue
			protein_1[transcript_path] += i.strip()
	protein_2 = {}
	for i in fastafile2_unpack:
		if i.startswith(">"):
			transcript = i.strip().split("|")[2].replace("gene=", "")
			if i not in protein_2:
				protein_2[i] = ""
				continue
			protein_2[transcript] += i.strip()
	common_proteins:List[List[Tuple]] = []
	for i in range(len(id1)):
		for val in range(len(id2)):
			for (val, key) in protein_1.items():
				for (val2,key2) in protein_2.items():
					if id1[i] == id2[val] == val == val2 and key !=key2:
						id1insert = id1[i]
						id1start = int(proteincoordinates1[i][0])
						id1end = int(proteincoordinates1[i][1])
						strand1insert = strand1[i]
						id2insert = id2[val]
						id2start = int(proteincoordinates2[i][0])
						id2end = int(proteincoordinates2[i][1])
						strand2insert = strand2[val]
						sequence_1 = protein_1[val]
						sequence_2 = protein_2[val2]
						insertvalue1 = (id1insert, id1start, id1end, strand1insert, sequence_1)
						insertvalue2 = (id2insert, id2start, id2end, strand2insert, sequence_2)
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
	barplot2 = sns.barplot(x = proteinnames2, y = proteinlength2, hue = proteinlength2)
	barplot1.get_figure().savefig("barplot-1.png")
	barplot2.get_figure().savefig("barplot-2.png")
	with open("comparative-common.txt", 'w') as filewrite:
		filewrite.write("The comparative table for the different genes and protein")
		filewrite.write("id1"+ '\t' + "id1start" + '\t' + "id1end" + '\t' + "id1strand" + '\t' "id1sequence" + '\t' + "id2" + '\t' + "id2start" + '\t' + "id2end" + '\t' + "id2stand" + '\t' + "id2sequence")
		for i in range(len(common_proteins)):
			filewrite.write(common_proteins[i][0][0] +'\t' + common_proteins[i][0][1] +'\t' + common_proteins[i][0][2] +'\t' + common_proteins[i][0][3] +'\t' + common_proteins[i][1][4] + '\t' + common_proteins[i][1][0] +'\t' + common_proteins[i][1][1] +'\t' +common_proteins[i][1][2]
				+'\t' + common_proteins[i][1][3] + '\t' + common_proteins[i][1][4])
		filewrite.close()
