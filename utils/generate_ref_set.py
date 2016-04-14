#!/usr/bin/env
#
# Lina Faller, April 2016
#
# Usage: python generate_ref_set.py -s V_REGION_START -e V_REGION_END -i REF_FASTA -o OUTFILE -a AMBIGUOUS_SEQUENCES
#
# This script reads in a fasta file and generates a set of unique probes, one for each input
# fasta sequence.
#

import sys
import os
import argparse
import logging
from Bio import SeqIO

########################################################################
# set up command line argument parsing
def parse_args():

	parser = argparse.ArgumentParser(description='This script reads in a fasta file and generates a set of unique probes, one for each input fasta sequence. If the start and/or end values are outside the range of the sequence, the script defaults to using the whole sequence for identifying appropriate probes.')

	parser.add_argument('-i', '--input', action="store", help='File containing the input.', type = str)
	parser.add_argument('-o', '--output', action="store", help='File containing the output.', type = str)
	parser.add_argument('-a', '--ambiguous', action="store", help='File containing the ambiguous sequences where no unique kmer could be found.', type = str)
	parser.add_argument('-s', '--start', action="store", help='Starting value of v-region to consider. Default value: 1 (i.e. start identifying probes at first base in sequence)', type=int, required = False, default = 1)
	parser.add_argument('-e', '--end', action="store", help='Ending value of v-region to consider. Default value: 0 (i.e. stop identifying probes at last base in sequence)', type=int, required = False, default = 0)
	parser.add_argument('-l', '--log', action="store", help='File containing logging information.', type = str )
	
	return parser.parse_args()

###############################################################################
# return all substrings of a given length k
def get_all_substrings( k, input_string ):
	length = len(input_string)
	list_of_strings = list()
	
	for i in xrange(length): # loop over the starting position
		current_substring = input_string[i:i+k]
		if len(current_substring) == k:
			list_of_strings.append( current_substring )
	return list_of_strings
	
########################################################################
# annotate_kmers
#
# Return a dict of kmers that exist exactly once in the fasta sequence
def annotate_kmers( header, seq, kmin, kmax ):

	current_dict_of_kmers = dict()
	
	for k in range(kmin, kmax):
	
		# get all substrings for this k
		current_kmers = get_all_substrings(k,seq)
		logging.debug( "k: {0}".format(k) )
		logging.debug( "Current kmers: {0}".format(current_kmers) )
		
		# add this set of kmers to the master dict
		for kmer in current_kmers:
			current_dict_of_kmers[ kmer ] = header
	
	return current_dict_of_kmers
	
########################################################################
# Main
def main():
	
	# parse command line arguments
	arguments = parse_args()

	# variables used
	list_of_kmers = dict()
	kmin = 15
	kmax = 40
	seq_start = int(arguments.start)
	seq_end = int(arguments.end)
	master_kmer_dict = dict()
	seq_to_kmer_dict = dict()
	seqs_in_input = dict()
	count_records = 0
	count_records_with_unique_kmer = 0
	
	# do you want to debug things?
	logging.basicConfig( filename = arguments.log, level = logging.INFO, format = '%(asctime)s %(message)s' )
		
	logging.info( "####################################" )
	logging.info( "Program started" )
	logging.info( "Call: " + str(sys.argv) )	
	logging.info( "k min: " + str(kmin) )
	logging.info( "k max: " + str(kmax) )
	
	########################################################################
	# process fasta file	
	for record in SeqIO.parse(arguments.input, "fasta"):
		
		header = record.description
		seq = str(record.seq)
		seq_len = len(seq)
		
		# remember that we saw this sequence
		seqs_in_input[ header ] = seq
		
		logging.info( "Annotating fasta sequence: {0}\tLength: {1}\t({2})".format(header, seq_len, seq) )
			
		# check if appropriate truncating values
		if seq_start < 0:
			logging.info( "Invalid start parameter ({0}). Use default of 0 instead.".format( seq_start ) )
			seq_start = 0
		if seq_end > seq_len:
			logging.info( "Invalid end parameter ({0}). Use default of length of sequence ({1}) instead.".format( seq_end, seq_len ) )
			seq_end = seq_len
			
		truncated_sequence = seq[seq_start:seq_end]
		logging.info( "Truncate the sequence to [{0},{1}]. Length now {2}.".format(seq_start, seq_end, len(truncated_sequence) ) )

		# annotate kmers for this sequence
		current_dict_of_kmers = annotate_kmers( header, truncated_sequence, kmin, kmax )
		
		# append this header to the master dict
		for kmer in current_dict_of_kmers:
			if kmer in master_kmer_dict:
				master_kmer_dict[ kmer ].append(header)
			else:
				master_kmer_dict[ kmer ] = [ header ]
		
		count_records += 1
	
	########################################################################
	# determine kmers unique for a single sequence
	logging.debug( "Master list:" )
	for kmer in master_kmer_dict:
		list_of_sequences = master_kmer_dict[kmer]
		num_seqs = len(list_of_sequences)
		
		logging.debug( "{0}\t{1}\t{2}".format( kmer, num_seqs, list_of_sequences ) )
	
		# if this kmer uniquely identifies one sequence, save it
		# note, currently if there are two unique kmers for a given sequence, it 
		# just saves the last one it found
		if num_seqs == 1:
			seq_to_kmer_dict[ list_of_sequences[0] ] = kmer
			
	logging.info( "Processed {0} fasta sequences.".format( count_records ) )

	########################################################################
	# determine kmers unique for a single sequence
	logging.info( "Write output file: {0}".format( arguments.output ) )
	
	OUT = open( arguments.output, 'w' ) or die( "Cannot open {0} for output!\n\n".format( arguments.output ) )
	
	for seq_id in seq_to_kmer_dict:
		# write output
		OUT.write( "{0}\t{1}\n".format( seq_id, seq_to_kmer_dict[ seq_id ] ) )
		
		# remove from list of known fasta headers
		del seqs_in_input[seq_id]
		
		# keep count
		count_records_with_unique_kmer += 1
	
	OUT.close()
	
	logging.info( "Found unique kmers for {0} sequences (out of {1} sequences).".format( count_records_with_unique_kmer, count_records ) )
	
	########################################################################
	# output the remaining sequences that did not have a unique kmer to identify them
	logging.info( "Write ambiguous fasta sequences to: {0}".format(arguments.ambiguous) )
	AMBIGUOUS = open( arguments.ambiguous, 'w' ) or die( "Cannot open {0} for output!\n\n".format( arguments.ambiguous ) )
	
	for seq_id in seqs_in_input:
		AMBIGUOUS.write( ">{0}\n{1}\n".format( seq_id, seqs_in_input[ seq_id ] ) )
		
	AMBIGUOUS.close()
	
if __name__ == '__main__':
	main()
