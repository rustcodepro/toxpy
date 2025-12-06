# toxpy

- python interface for ToxDB with annotation figures and tables.
- From bioinformatics data analysis to machine and deep learning model.
- a complete package which performs the ToxDB comparison and also fits a LSTM and Convolutional Neural Network with L2 regularization
```
Usage: tox.py [OPTIONS] COMMAND [ARGS]...

  Command line interface for comparing protein sequences from GFF and FASTA
  files for ToxDB. This also includes the machine and deep learning using the
  Keras and Tensorflow and integrated LSTM and Convolution Network with L2
  Regularization.

  Gaurav Sablok codeprog@icloud.com

Options:
  --help  Show this message and exit.

Commands:
  toxcnn-l2-regularization-option
                                  implement the CNN on the PTM
  toxcompare-a-option             same-protein
  toxcompare-b-option             different-protein
  toxcompare-c-option             samprotein-same-strand
  toxcompare-d-strand-option      sameprotein-different-strand
  toxcompare-e-option             same-protein-same-seq
  toxcompare-f-option             same-protein-different-seq
  toxcomplare-g-option            compare mrna and seq
  toxcomplare-h-option            compare mrna for the comparative analysis
  toxlstm-option                  implement the lstm on the PTM


```

```
id1	id2	id1start	id2start	id1end	id2end	id1strand1	idstrand2
TGME49_302060	1	454	-	TGME49_302060	1	454	-
TGME49_302057	513	956	-	TGME49_302057	513	956	-
TGME49_302055	977	1342	-	TGME49_302055	977	1342	-
TGME49_300601	1349	1756	-	TGME49_300601	1349	1756	-
```
![](https://github.com/rustcodepro/toxpy/blob/main/barplot-1.png)

![](https://github.com/rustcodepro/toxpy/blob/main/barplot-2.png)


Gaurav Sablok \
codeprog@icloud.com
