#!/usr/bin/env
#
# Lina Faller, April 2016
#
# Usage: python annotate.py -r REF_FILE -i FASTA -o OUTFILE -l LOG
#
# This script reads in a reference file with two columns. The first column contains an ID
# and the second column contains a short fragment of nucleotide sequence that uniquely
# identifies this specific reference sequence.
# 
# The script outputs a list of all the reference IDs and the counts indicating how many
# sequences in the input matched this reference. There will be one additional ID, labeled
# "unmatched" which indicates the number of reads that could not be matched.
#

import sys
import os
import argparse
import logging
from Bio import SeqIO

########################################################################
# set up command line argument parsing
def parse_args():

	parser = argparse.ArgumentParser(description="This script reads in a reference file with two columns. The first column contains an ID and the second column contains a short fragment of nucleotide sequence that uniquely identifies this specific reference sequence.\nThe script outputs a list of all the reference IDs and the counts indicating how many sequences in the input matched this reference. There will be one additional ID, labeled 'unmatched' which indicates the number of reads that could not be matched.")

	parser.add_argument('-i', '--input', action="store", help='File containing the input.', type = str)
	parser.add_argument('-o', '--output', action="store", help='File containing the output.', type = str)
	parser.add_argument('-r', '--reference', action="store", help='File containing the pairs of reference id and unique probe.', type = str)
	parser.add_argument('-l', '--log', action="store", help='File containing logging information.', type = str )
	
	return parser.parse_args()

########################################################################
# find_probe
#
# Return a string indicating the reference_id for this sequence. 
# Return "unmatched" if no reference_id was found.
def find_probe( sequence, reference_ids ):
	
	logging.info( "Looking at sequence: {0}".format( sequence ) )
	
	for ref_id in reference_ids:
		
		current_probe = reference_ids[ ref_id ]
				
		if current_probe in sequence:
			logging.info( "Identified reference id ({0}) with probe ({1})".format( ref_id, current_probe ) )
			return ref_id
	
	logging.info( "No reference id was identified. Probe is considered unmatched." )		
	return 'unmatched'
	
########################################################################
# process_reference_file
#
# Return a dict where:
#     key = reference_id
#     value = unique probe sequence
def process_reference_file( REFERENCE ):

	reference_ids = dict()
	
	with open( REFERENCE ) as IN:
		for line in IN:
			line = line.strip()
			items = line.split( "\t" )
			reference_ids[ items[0] ] = items[1]

	return reference_ids
	
########################################################################
# Main
def main():
	
	# parse command line arguments
	arguments = parse_args()

	# variables used
	all_probes = dict()
	
	# do you want to debug things?
	logging.basicConfig( filename = arguments.log, level = logging.INFO, format = '%(asctime)s %(message)s' )
		
	logging.info( "####################################" )
	logging.info( "Program started" )
	logging.info( "Call: " + str(sys.argv) )	

	########################################################################
	# process reference file	
	ref_ids = process_reference_file( arguments.reference )
	
	########################################################################
	# process fasta file	
	for record in SeqIO.parse(arguments.input, "fasta"):
		
		header = record.description
		seq = str(record.seq)
		seq = seq.strip()
		
		current_ref_id = find_probe( seq, ref_ids )
		
		if current_ref_id in all_probes:
			all_probes[ current_ref_id ] += 1
		else:
			all_probes[ current_ref_id ] = 1
		
	########################################################################
	# write output
	logging.info( "Write output file: {0}".format( arguments.output ) )
	
	OUT = open( arguments.output, 'w' ) or die( "Cannot open {0} for output!\n\n".format( arguments.output ) )
	
	for probe in all_probes:
		# write output
		OUT.write( "{0}\t{1}\n".format( probe, all_probes[ probe ] ) )
	
	OUT.close()
	
	logging.info( "Program finished running." )
	
if __name__ == '__main__':
	main()
