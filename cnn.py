from typing import Dict, List

import numpy as np
import tensorflow as tf
from sklearn.metrics import classification_report
from sklearn.model_selection import train_test_split
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.layers import Conv1D, Dense, Dropout, Flatten, MaxPooling1D
from tensorflow.keras.models import Sequential
from tensorflow.keras.preprocessing.sequence import pad_sequences
from tensorflow.keras.regularizers import l2

# Gaurav Sablok
# codeprog@icloud.com

def toxcnn_l2_regularization(dnaseq, classificationfile, ptmsequence):
	"""
	 This model fits the convolutional neural network(CNN) with L2 regularization to minimize the vanishing gradient
		and optimized for the Post Translational Modifications (PTM)
	"""
	seqdict:Dict[str, str] =  {}
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
	X_train, _, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
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
	_ = model.fit(X_train, y_train, epochs=50, batch_size=32, validation_split=0.2,
                    callbacks=[early_stop], verbose=1)
	y_pred = (model.predict(ptmseq) > 0.5).astype(int)
	print(classification_report(y_test, y_pred))
