#!/usr/bin/env python3

import time
import argparse
import os
import re
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation
import samtools_lookup
import pandas as pd

# VARIABLES
genes = "/home/sfrenk/Documents/Resources/Seq/WS251/genes.gtf"
transposons = "/home/sfrenk/Documents/Resources/Seq/transposons/ce11_repbase/ce11_transposons.gff3"

###############################################################################

parser = argparse.ArgumentParser(description="Sequence and annotation from entrez and UCSC and create Ugene project folder")

parser.add_argument("-c", "--chromosome")
parser.add_argument("-s", "--start", type = int)
parser.add_argument("-e", "--end", type = int)
parser.add_argument("-n", "--name", help = "Name of sequence (by default, the name is constructed from the coordinates)", default = None)
parser.add_argument("-d", "--directory", help = "Output directory (default: current directory)", default = ".")
parser.add_argument("-f", "--feature", help = "optional gff file containing genomic coordinates of extra features", default = None)

args = parser.parse_args()


# Set the master directory where project directories are kept
#ugene_dir = "/home/sfrenk/Documents/Ugene/"
#seq_dir = ugene_dir + args.name + "/"

#if os.path.exists(ugene_dir + seq_dir):
#	print("ERROR: " + ugene_dir + args.name + " already exists. Choose a different name.")
#	sys.exit()
#else:
#	os.mkdir(ugene_dir + args.name)

# GET SEQUENCE
seq_object = samtools_lookup.get_seq(args.chromosome, args.start, args.end).seq

# Create a record
coords = args.chromosome + ":" + str(args.start) + "," + str(args.end)
record = SeqRecord(seq_object,
                   id = coords,
                   name = args.name,
                   description='Created using ugene_fetch.py')

#GET FEATURES

feature_files = [genes, transposons]

if args.feature:
	feature_files.append(args.feature)

feature_regex = re.compile("(gene_name |Name=)([^;]*);")

for i in feature_files:
	table = pd.read_csv(i, sep = "\t", header = None, usecols = [0,3,4,6,8], names = ["chrom", "start", "end", "strand", "feature"])

	# Find features within the specified genomic range
	table = table.loc[(table.chrom == args.chromosome) & (table.start > args.start) & (table.end < args.end),]

	for k, row in table.iterrows():
		
		name = str(re.search(feature_regex, row.feature).group(2))
		
		if row.strand == "-":
			strand = -1
		else:
			strand = 1

		feature = SeqFeature(FeatureLocation(start = (row.start - args.start), end = (row.end - args.start) + 1), type = 'misc_feature', qualifiers = {"label" : name, "ugene_name" : name}, strand = strand)
		record.features.append(feature)

# Save output
output_filename = args.directory + "/" + args.name + ".gb"
output_file = open(output_filename, 'w')
SeqIO.write(record, output_file, 'genbank')
output_file.close()
