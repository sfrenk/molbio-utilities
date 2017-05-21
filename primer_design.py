#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import primer3
import argparse

class primerSet(object):

    """ A class to store multiple sets of primers. This class has the following attributes:

        fwd: F1 primer sequence
        rev: R1 primer sequence
        size: product size
        fwd_tm: F1 Tm
        rev_tm: R1_Tm
        fwd_gc: F1 GC content (%)
        rev_gc: R1 GC content (%)
    """

    def __init__(self, fwd, rev, size, fwd_tm, rev_tm, fwd_gc, rev_gc):
        
        self.fwd = fwd
        self.rev = rev
        self.product_size = size
        self.fwd_tm = fwd_tm
        self.rev_tm = rev_tm
        self.fwd_gc = fwd_gc
        self.rev_gc = rev_gc


    def summary(self):

        """ Return primer set info as a string
        """

        return("Fwd: " + self.fwd + "\n" + "Rev: " + self.rev + "\n" + "Product size: " + str(self.size) + "\n" + "Fwd Tm: " + str(self.fwd_tm) + "\n" + "Rev Tm: " + str(self.rev_tm) + "\n" + "Fwd GC: " + str(self.fwd_gc) + "\n" "Rev GC: " + str(self.fwd_gc) + "\n")


def primerSet_from_primer3(p3, p_index):

    """ Extract a primer pair from a primer3 output object and return a primerSet object for this pair
    """
        
    fwd_seq = p3['PRIMER_LEFT_' + str(p_index) + '_SEQUENCE']
    rev_seq = p3['PRIMER_RIGHT_' + str(p_index) + '_SEQUENCE']
    product_size = p3['PRIMER_PAIR_' + str(p_index) + '_PRODUCT_SIZE']
    fwd_TM = round(p3['PRIMER_LEFT_' + str(p_index) + '_TM'], 1)
    rev_TM = round(p3['PRIMER_RIGHT_' + str(p_index) + '_TM'], 1)
    fwd_GC = round(p3['PRIMER_LEFT_' + str(p_index) + '_GC_PERCENT'], 1)
    rev_GC = round(p3['PRIMER_RIGHT_' + str(p_index) + '_GC_PERCENT'], 1)


    return(primerSet(fwd_seq, rev_seq, product_size, fwd_TM, rev_TM, fwd_GC, rev_GC))


def design_primers(seq, size = None, target = None):

    """ Run pimer3 on input sequence, specifying target amplicon size. Returns a list of three primer sets.

        Note: the optional target parameter is in the form of a list [x,y] where x = target start coordinate and y = target length
    """

    seq_params = {
                    'SEQUENCE_ID': "input_sequence",
                    'SEQUENCE_TEMPLATE': seq 
                }
    
    primer_params = {
                    'PRIMER_TASK': 'generic',
                    'PRIMER_PICK_LEFT_PRIMER': 1,
                    'PRIMER_PICK_INTERNAL_OLIGO': 0,
                    'PRIMER_PICK_RIGHT_PRIMER': 1,
                    'PRIMER_NUM_RETURN': 3,
                    'PRIMER_OPT_SIZE': 20,
                    'PRIMER_MIN_SIZE': 18,
                    'PRIMER_MAX_SIZE': 25,
                    'PRIMER_OPT_TM': 60.0,
                    'PRIMER_MIN_TM': 57.0,
                    'PRIMER_MAX_TM': 63.0,
                    'PRIMER_MIN_GC': 20.0,
                    'PRIMER_MAX_GC': 80.0,
                    'PRIMER_MAX_POLY_X': 5,
                    'PRIMER_SALT_MONOVALENT': 50.0,
                    'PRIMER_DNA_CONC': 50.0,
                    'PRIMER_MAX_NS_ACCEPTED': 0,
                    'PRIMER_MAX_SELF_ANY': 12,
                    'PRIMER_MAX_SELF_END': 8,
                    'PRIMER_PAIR_MAX_COMPL_ANY': 12,
                    'PRIMER_PAIR_MAX_COMPL_END': 8
                    }

    if target != None:
        seq_params['SEQUENCE_TARGET'] = target
        if size == None:

            # If target is specified but there is no size limit, set max product size to target length + 200bp

            primer_params['PRIMER_PRODUCT_SIZE_RANGE'] = [[target[1], target[1] + 200]]


    if size != None:
        primer_params['PRIMER_PRODUCT_SIZE_RANGE'] = size


    primer = (primer3.bindings.designPrimers(seq_params, primer_params))
    
    primer_object = [primerSet_from_primer3(primer, x) for x in range(3)]
    
    return(primer_object)


def main():

    parser = argparse.ArgumentParser(description="Get primers to amplify from a given template")

    parser.add_argument("file", help = "Input file")
    parser.add_argument("-f", "--format", help = "Input file format (eg fasta, genbank). Default: fasta", default = "fasta")
    parser.add_argument("-s", "--size", help = "Size of PCR product (default: 500bp)", default = 500, type = int)

    args = parser.parse_args()

    record = SeqIO.read(args.file, args.format)

    input_seq = str(record.seq)

    primer_seq = design_primers(input_seq, args.size)
    for i in range(3):

        print("\nprimer set " + str(i) + "\n" + primer_seq[i].summary())

if __name__ == '__main__':
    main()
