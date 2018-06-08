#!/usr/bin/env python3

import time
import argparse
import os
import re
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, generic_dna
from Bio.SeqFeature import SeqFeature, FeatureLocation
import samtools_lookup
import pandas as pd

# VARIABLES
genes = "/home/sfrenk/Documents/Resources/Seq/WS251/genes.gtf"
transposons = "/home/sfrenk/Documents/Resources/Seq/transposons/ce11_repbase/ce11_transposons.gff3"

###############################################################################

def add_feature(input_file, gff_file):

	record = SeqIO.read(input_file, "gb")
	record.seq.alphabet= generic_dna

	chrom, start, end = get_coords_from_record(record)

	# Get features
	gff = open(gff_file, "r")
	for line in gff:

		feature_chrom = line.strip().split("\t")[0]
		feature_group = line.strip().split("\t")[1]
		feature_start = int(line.strip().split("\t")[3])
		feature_end = int(line.strip().split("\t")[4])
		feature_strand = line.strip().split("\t")[6]
		feature_name = re.search("(name|gene_id)[= ]?([^;]+)", line.strip().split("\t")[8]).group(2)

		if feature_chrom == chrom and feature_start < start and feature_end > end:

			annot_start = (feature_start-start)
			annot_end = (feature_end-start)

			if feature_strand == "-":
				annot_strand = -1
			else:
				annot_strand = 1

			feature = SeqFeature(FeatureLocation(start = (annot_start), end = (annot_end) + 1), type = 'misc_feature', qualifiers = {"label" : feature_name, "ugene_name" : feature_name, "ugene_group" : feature_group}, strand = annot_strand)
			record.features.append(feature)

	gff.close()

	return(record)

def get_coords_from_record(record):

	# Get chromosomal coords of Ugene file from the ID

	chrom = re.search("^([^:]+)", record.id).group(1)
	start = int(re.search(":([0-9]+),", record.id).group(1))
	end = int(re.search(",([0-9]+)$", record.id).group(1))

	return(chrom, start, end)

def get_features(record, feature = None):

	# Get coordinates of record
	chrom, start, end = get_coords_from_record(record)

	feature_files = [genes, transposons]

	if feature == None:
		feature_files.append(feature)

		feature_regex = re.compile("(gene_name |Name=)([^;]*)")

		for i in feature_files:
			table = pd.read_csv(i, sep = "\t", header = None, usecols = [0,3,4,6,8], names = ["chrom", "start", "end", "strand", "feature"])

			# Find features within the specified genomic range
			table = table.loc[(table.chrom == chrom) & (table.start > start) & (table.end < end),]

			for k, row in table.iterrows():
				
				name = str(re.search(feature_regex, row.feature).group(2))
				
				if row.strand == "-":
					strand = -1
				else:
					strand = 1

				feature = SeqFeature(FeatureLocation(start = (row.start - args.start), end = (row.end - args.start) + 1), type = 'misc_feature', qualifiers = {"label" : name, "ugene_name" : name}, strand = strand)
				record.features.append(feature)

	return(record)


def get_sequence(chrom, start, end, name = None):

	seq_object = samtools_lookup.get_seq(chrom, start, end).seq

	# Create a record
	coords = chrom + ":" + str(start) + "," + str(end)
	record = SeqRecord(seq_object,
	                   id = coords,
	                   name = name,
	                   description='Created using ugene_fetch.py')

	return(record)


def save_output(record, output_filename):

	output_file = open(output_filename, 'w')
	SeqIO.write(record, output_file, 'genbank')
	output_file.close()


def main():

	parser = argparse.ArgumentParser(description = "Creates Ugene file with sequence and features using coordinates provided")

	parser.add_argument("-c", "--chromosome")
	parser.add_argument("-s", "--start", type = int)
	parser.add_argument("-e", "--end", type = int)
	parser.add_argument("-n", "--name", help = "Name of sequence (by default, the name is constructed from the coordinates)", default = None)
	parser.add_argument("-d", "--directory", help = "Output directory (default: current directory)", default = ".")
	parser.add_argument("-f", "--feature", help = "Optional gff file containing genomic coordinates of extra features", default = None)
	parser.add_argument("-a", "--add_feature", help = "Add features from gff file (using -f/--feature argument) to existing Ugene file", default = None)

	args = parser.parse_args()

	if args.add_feature != None:

		ugene_record = add_feature(args.add_feature, args.feature)
		output_name = args.add_feature

	else:

		if args.name == None:
			name = chrom + "_" + str(start) + "_" + str(end)
		else:
			name = args.name

		seq = get_sequence(args.chromosome, args.start, args.end, name)
		ugene_record = get_features(seq, args.feature)
		output_name = args.directory + "/" + name + ".gb"
	
	save_output(ugene_record, output_name)


if __name__ == '__main__':

	main()
