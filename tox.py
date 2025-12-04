#! usr/bin/env/ python3

# Gaurav Sablok
# codeprog@icloud.com

# Toxpy: python version for the comparative genomics of the ToxDB species
# It provided you the complete machine leaning as well as comparative genomics.
# a complete package which performs the ToxDB comparison and also fits a LSTM
# and Convolutional Neural Network with L2 regularization
# why i implemented L2 regulrization as PTM have the end terminal
# modification and hence they will reach the vanishing gradient earlier.



import click

from cnn import toxcnn_l2_regularization
from differentprotein import toxcompare_different_proteins
from differentproteinsameseq import toxcompare_different_protein_same_sequences
from lstm import tox_lstm
from sameproteindifferentseq import toxcompare_same_proteins_different_sequences
from sameproteindiffstrand import toxcompare_same_proteins_different_strand
from sameproteins import toxcompare_same_proteins
from sameproteinsamstrand import toxcompare_same_proteins_with_same_strand


@click.group()
def cli():
    """Command line interface for comparing protein sequences from GFF and FASTA files for ToxDB.
       This also includes the machine and deep learning using the Keras and Tensorflow and
       integrated LSTM and Convolution Network with L2 Regularization.

       Gaurav Sablok
       codeprog@icloud.com
    """
    pass

@cli.command()
@click.argument('pathfile1', type=click.Path(exists=True))
@click.argument('pathfile2', type=click.Path(exists=True))
def toxcompare_a_option(pathfile1, pathfile2):
	"""
	same-protein


	This function takes the path of the ToxannotationDB file
	and compared for the sequences of the protein coding.
	This compares the annotation records of the two gff files
	and show the similar geneid and their associated information."""
	toxcompare_same_proteins(pathfile1, pathfile2)

@cli.command()
@click.argument('pathfile1', type=click.Path(exists=True))
@click.argument('pathfile2', type=click.Path(exists=True))
def toxcompare_b_option(pathfile1, pathfile2):
	"""
	different-protein


 	This function takes the path of the ToxannotationDB file
	and compared for the sequences of the protein coding.
 	This compares the two annotation files and shows the missing
	proteins which are not present in them and their coordinates.

	"""
	toxcompare_different_proteins(pathfile1, pathfile2)


@cli.command()
@click.argument('pathfile1', type=click.Path(exists=True))
@click.argument('pathfile2', type=click.Path(exists=True))
def toxcompare_c_option(pathfile1, pathfile2):
	"""
	samprotein-same-strand


 	This function takes the path of the ToxannotationDB file
	and compared for the sequences of the protein coding.
 	This compares the annotation records of the two gff files and
	show the similar geneid and similar strand and their associated information.
	"""
	toxcompare_same_proteins_with_same_strand(pathfile1, pathfile2)


@cli.command()
@click.argument('pathfile1', type=click.Path(exists=True))
@click.argument('pathfile2', type=click.Path(exists=True))
def toxcompare_d_strand_option(pathfile1, pathfile2):
	"""
	sameprotein-different-strand


	This function takes the path of the ToxannotationDB file
	and compared for the sequences of the protein coding.
 	This compares the annotation records of the two gff files and
	show the similar geneid and different strands their associated information.
	"""
	toxcompare_same_proteins_different_strand(pathfile1, pathfile2)


@cli.command()
@click.argument('gff_file1', type=click.Path(exists=True))
@click.argument('gff_file2', type=click.Path(exists=True))
@click.argument('fastafile1', type=click.Path(exists=True))
@click.argument('fastafile2', type=click.Path(exists=True))
def toxcompare_e_option(gff_file1, gff_file2, fastafile1, fastafile2):
    """
    same-protein-same-seq


    This function takes the path of the ToxannotationDB file and compared for the sequences of the protein coding.
	  This compares the annotation records of the two gff files and
		show the similar geneid and their associated information.
	  This is to check the annotation errors if by chance two ids have the same sequence."""
    toxcompare_different_protein_same_sequences(gff_file1, gff_file2, fastafile1, fastafile2)

@cli.command()
@click.argument('gff_file1', type=click.Path(exists=True))
@click.argument('gff_file2', type=click.Path(exists=True))
@click.argument('fastafile1', type=click.Path(exists=True))
@click.argument('fastafile2', type=click.Path(exists=True))
def toxcompare_f_option(gff_file1, gff_file2, fastafile1, fastafile2):
    """
    same-protein-different-seq


    This function takes the path of the ToxannotationDB file and compared for the sequences of the protein coding.
	  This compares the annotation records of the two gff files and
		show the similar geneid and their associated information.
	  This compares the gff and the fasta to have the same ids and the same sequences."""
    toxcompare_same_proteins_different_sequences(gff_file1, gff_file2, fastafile1, fastafile2)

@cli.command()
@click.argument('dnaseq', type=click.Path(exists=True))
@click.argument('classificationfile', type=click.Path(exists=True))
@click.argument('ptmsequences', type=click.Path(exists=True))
def toxlstm_option(dnaseq, classificationfile, ptmsequences):
	"implement the lstm on the PTM"
	tox_lstm(dnaseq, classificationfile, ptmsequences)

@cli.command()
@click.argument('dnaseq', type=click.Path(exists=True))
@click.argument('classificationfile', type=click.Path(exists=True))
@click.argument('ptmsequences', type=click.Path(exists=True))
def toxcnn_l2_regularization_option(dnaseq, classificationfile, ptmsequences):
	"implement the CNN on the PTM"
	toxcnn_l2_regularization(dnaseq, classificationfile, ptmsequences)


if __name__ == '__main__':
    cli()
