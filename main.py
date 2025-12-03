# Gaurav Sablok
# codeprog@icloud.com

from keras import cnn_l2_regularization, lstm
from tox import toxannotator


def main(pathfile1, pathfile2):

	toxannotator(pathfile1, pathfile2)
	lstm()
	cnn_l2_regularization()

	if __name__ == main:
		main(pathfile1, pathfile2)
