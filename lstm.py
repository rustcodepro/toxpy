from typing import Dict, List

import numpy as np
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.layers import LSTM, Dense, Embedding
from tensorflow.keras.models import Sequential
from tensorflow.keras.preprocessing.sequence import pad_sequences
from tensorflow.keras.regularizers import l2

# Gaurav Sablok
# codeprog@icloud.com

def tox_lstm(dnaseq,classificationfile, ptmsequence):
	"""
	  This function fits the LSTM model for the classification of the
			PTM(post translational modifications) using the Keras and Tensorflow
			LSTM layers. Optimized for the PTM.
	"""
	seqdict =  {}
	headerdnaseq:List[str] = []
	sequencednaseq:List[str] = []
	with open(dnaseq, 'r') as pathopenfile:
		for i in pathopenfile.readlines():
			if i.startswith(">"):
				headerdnaseq.append(i.strip().replace(">", ""))
			elif not i.startswith(">"):
				sequencednaseq.append(i.strip())
			else:
				pass
	for i in range(len(headerdnaseq)):
		seqdict[headerdnaseq[i]] = sequencednaseq[i]
	classification_id:List[int] = []
	with open(classificationfile, 'r') as clasificationids:
		for i in clasificationids.readlines():
			line = i.strip()
			classification_id.append(int(line))
	header:List[str] = []
	sequence:List[str] = []
	seqdictptm:Dict =  {}
	with open(ptmsequence, 'r') as pathopenfile:
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
	ptmseq:List[str] = [sequences for sequences in ptmsequence.values()]
	model = Sequential([
    	Embedding(input_dim=5, output_dim=8, input_length=max_length),
     LSTM(16, return_sequences=False),
     Dense(8, activation='relu'),
     Dense(1, activation='sigmoid')
])
	model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
	model.fit(X, y, epochs=10, batch_size=2, verbose=1)
	testsequence = [encode_dna_sequence(seq) for seq in ptmseq]
	max_length = max(len(seq) for seq in testsequence)
	X = pad_sequences(encoded_sequences, maxlen=max_length, padding='post')
	sequence = pad_sequences([testsequence], maxlen=max_length, padding='post')
	prediction = model.predict(sequence)
	print(f"Prediction for sequence: {'Positive' if prediction[0] > 0.5 else 'Negative'}")
