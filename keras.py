import random

import numpy as np
from sklearn.preprocessing import LabelEncoder
from tensorflow.keras.layers import LSTM, Dense, Embedding
from tensorflow.keras.models import Sequential
from tensorflow.keras.preprocessing.sequence import pad_sequences
from typing_extensions import Dict, List

# Gaurav Sablok
# codeprog@icloud.com

# LSTM model
# Gaurav Sablok
# codeprog@icloud.com
# a complete package which performs the ToxDB comparison and also fits a LSTM
# and Convolutional Neural Network with L2 regularization
# why i implemented L2 regulrization as PTM have the end terminal
# modification and hence they will reach the vanishing gradient earlier.
#

seqdict:Dict =  {}
classification_id:List[int] = []
seqdictptm:Dict =  {}
pathfile = ""
classification = ""
def fastaread(pathfile):
	header:List[str] = []
	sequence:List[str] = []
	with open(pathfile, 'r') as pathopenfile:
		for i in pathopenfile.readlines():
			if i.startswith(">"):
				header.append(i.strip().replace(">", ""))
			elif not i.startswith(">"):
				sequence.append(i.strip())
			else:
				pass
	for i in range(len(header)):
		seqdict[header[i]] = sequence[i]
		return seqdict

def classificationid(pathfileclassification):
	with open(pathfileclassification, 'r') as clasificationids:
		for i in clasificationids.readlines():
			line = i.strip()
			classification_id.append(int(line))
	return classification_id

def ptmclass(pathfile):
	header:List[str] = []
	sequence:List[str] = []
	with open(pathfile, 'r') as pathopenfile:
		for i in pathopenfile.readlines():
			if i.startswith(">"):
				header.append(i.strip().replace(">", ""))
			elif not i.startswith(">"):
				sequence.append(i.strip())
			else:
				pass
	for i in range(len(header)):
		seqdictptm[header[i]] = sequence[i]
		return seqdict


def encode_dna_sequence(seq):
    mapping = {'A': 1, 'C': 2, 'G': 3, 'T': 4}
    return [mapping[base] for base in seq]


# Toxdb LSTM model
def lstm():
	finalsequences:List[str] = [seq for seq in seqdict.values()]
	encoded_sequences = [encode_dna_sequence(seq) for seq in finalsequences]
	max_length = max(len(seq) for seq in finalsequences)
	X = pad_sequences(encoded_sequences, maxlen=max_length, padding='post')
	y = np.array(classification_id)
	ptmsequences:List[str] = [sequences for sequences in seqdictptm.values()]
	model = Sequential([
    Embedding(input_dim=5, output_dim=8, input_length=max_length),
    LSTM(16, return_sequences=False),
    Dense(8, activation='relu'),
    Dense(1, activation='sigmoid')  # Binary classification
])
	model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
	model.fit(X, y, epochs=10, batch_size=2, verbose=1)
	test_sequence = [encode_dna_sequence(seq) for seq in ptmsequences]
	max_length = max(len(seq) for seq in test_sequence)
	X = pad_sequences(encoded_sequences, maxlen=max_length, padding='post')
	sequence = pad_sequences([test_sequence], maxlen=max_length, padding='post')
	prediction = model.predict(sequence)
	print(f"Prediction for sequence: {'Positive' if prediction[0] > 0.5 else 'Negative'}")

#CNN model with L2 regularization
# Gaurav Sablok
# codeprog@icloud.com

def cnn_l2_regularization():
	finalsequences:List[str] = [seq for seq in seqdict.values()]
	encoded_sequences = [encode_dna_sequence(seq) for seq in finalsequences]
	max_length = max(len(seq) for seq in finalsequences)
	X = pad_sequences(encoded_sequences, maxlen=max_length, padding='post')
	y = np.array(classification_id)
	ptmsequences:List[str] = [sequences for sequences in seqdictptm.values()]
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
	y_pred = (model.predict(ptmsequences) > 0.5).astype(int)
	print(classification_report(y_test, y_pred))
