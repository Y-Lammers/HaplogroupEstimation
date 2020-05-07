#!/usr/bin/env python

# Usage: OutputVariants.py --sam [input SAM file] --vcf [SNP list] 
# --window [windows size]> [output table]

# author: Youri Lammers
# contact: youri.lammers@uit.no / youri.lammers@gmail.com

# import the modules used by the script
import os, argparse, itertools, sys, re
from collections import defaultdict

# Retrieve the commandline arguments
parser = argparse.ArgumentParser(description=
	"""Read a SAM and vcf SNP file and output the haplotypes counts
	across the SNP windows of the provided size""")

parser.add_argument('--sam', type=argparse.FileType('r'), nargs='?',
	default=sys.stdin, help='SAM file, use - in order to read from stdin')
parser.add_argument('--vcf', metavar='', type=str,
	help='SNP file in vcf format')
parser.add_argument('--window', metavar='', type=int, default=35,
	help='Window size in which to look for SNPs. Default=35')
args = parser.parse_args()


def parse_SNP():

	# Extract the SNP positions from the vcf file.

	# Create the lists for storing the raw SNP positions and the
	# SNP windows
	SNPlist, SNPwindows = defaultdict(list), defaultdict(list)

	# parse through the snp file
	for line in open(args.vcf):

		# skip if it is the header
		if line[0] == '#': continue

		# split the tsv format
		line=line.strip().split('\t')

		# append the SNP to the list
		SNPlist[line[0]].append(int(line[1]))


	# for each reference detected in the SNP list	
	for ref in SNPlist:

		# Parse through the list of SNPs 
		for pos in range(0,len(SNPlist[ref])):

			# store the SNPs for any window in a temporary list
			# and keep track of the window length
			tempwindow, inc = [SNPlist[ref][pos]], 1

			# keep adding SNPs to the temporary window so long as
			# the maximum window size is not exceeded
			try:
				while (SNPlist[ref][pos+inc] - 
					SNPlist[ref][pos]) <= args.window:

					tempwindow.append(SNPlist[ref][pos+inc])
					inc+=1
			except:
				pass

			# if at least two SNPs are in a temporary window
			# continue to saving the window information
			if len(tempwindow) >= 2:

				# if the reference already has windows recorded,
				# check the amount of overlap between the old
				# and new window prior to saving
				if len(SNPwindows[ref]) > 0:

					# if two or more SNP positions are
					# shared with the last window added,
					# skip the window as it would inflate
					# the number of windows
					if len(list(set(tempwindow) & 
						set(SNPwindows[ref][-1]))) >= 2:
						continue

				# save the window
				SNPwindows[ref].append(tempwindow)

	# return the window information
	return SNPwindows


def translate_cigar(cigar, seq):

	# use the sequence and cigar code to reconstruct the read alignment
	
	# create an empty string and extract the individual cigar components
	cseq, match = '', re.findall(r'\d+[A-Z]+',cigar)

	# for each cigar component, build up the read
	for i in match:

		# get the length and cigar info for each component
		clength, ctype = int(i[:-1]), i[-1]

		# build up the read according to the cigar type
		if ctype == "M":
			cseq += seq[:clength]
			seq = seq[clength:]
		elif ctype == "I":
			seq = seq[clength:]
		elif ctype == "D":
			cseq += "*"*clength
		else:
			pass

	# return the build up sequence
	return cseq


def parse_sam(SNPwindows):

	# parse through the SAM file and see which reads are overlapping
	# with known SNP windows. If a read is overlapping, translate the
	# aligned read and extract the haplotype present.

	# store the usable reads in a new dictionary
	haplo_reads = defaultdict(list)

	# open the sam file
	for line in args.sam:

		# split the tsv formated sam file		
		line = line.strip().split('\t')

		# get the reference info
		ref = line[2]

		# if the reference contains no SNP windows, skip the read
		if ref not in SNPwindows: continue

		# get the start and stop info, the read and obtain the 
		# translated cigar sequence
		start, cigar, seq = int(line[3]),line[5],line[9]
		cseq = translate_cigar(line[5],line[9])

		# calculate the length and end of the read
		length = len(cseq)
		end = start + length

		# go through the windows recorded for the reference
		for snp in SNPwindows[ref]:

			# if the read overlaps a snp window, store it.
			if start < snp[0] and end > snp[-1]:

				# store it with a key that is build up out of
				# the reference and the SNP positions
				window = '{0}:{1}'.format(ref,
					'-'.join([str(p) for p in snp]))

				haplo_reads[window].append([start,length,
					end,cigar,seq,cseq])

	# keep track of the total number of haplotypes and the variants present
	tot_haplo_counts, tot_haplo_vars = defaultdict(int), {}

	# Parse through the SNP windows with reads
	for sset in haplo_reads:

		# get the reference, total coverage and SNP info from the window
		ref = sset.split(':')[0]
		totcov = len(haplo_reads[sset])
		snp = [int(p) for p in sset.split(':')[1].split('-')]

		# keep track of how often each haplotype is present
		haplo_count = defaultdict(int)

		# parse through the reads
		for read in haplo_reads[sset]:

			# build up the present haplotype by appending the
			# recorded nucleotides for each SNP in the window
			var = ''

			for s in snp:

				# append the variant present
				var += read[5][(s-read[0])]

			# store the haplotype detected and keep track of
			# the abundance
			haplo_count[var] += 1

		# remove low abundance variants (less than 3 occurences or
		# a frequency lower than 0.15)
		#varlist = list(haplo_count.keys())
		for var in list(haplo_count.keys()):
			if haplo_count[var] < 3 or float(haplo_count[var])/totcov <= 0.15:
				del haplo_count[var]

		# if no haplotypes survived the above filtering, continue with
		# the next haplotype window
		if len(haplo_count) == 0: continue

		# Add the remaining haplotypes and count information 
		# to the total haplotype dictionaries
		tot_haplo_counts[len(haplo_count)] += 1
		tot_haplo_vars[sset] = haplo_count


	# Get the maximum number of haplotypes recorded
	maxcount = max(tot_haplo_counts.keys())
	
	# Get all the detected haplotype positions and sort them
	# print the haplotypes and count for each snp window to the stdout

	# get the windows and sort them
	hapwinkeys = list(tot_haplo_vars.keys())
	hapwinkeys.sort()

	# for each window, set up the output info
	for hap in hapwinkeys:

		# create a list with the haplotype info and counts		
		var = ['\t'.join([i,str(tot_haplo_vars[hap][i])])
				for i in tot_haplo_vars[hap]]

		# append blanks if the number of variants is below the max
		# for the entire dataset
		while len(var) < maxcount: var.append('')

		# format and print the results to the stdout
		print('{0}\t{1}'.format(hap.replace(':','\t'),'\t'.join(var)))

def main ():

	# get the SNP window information
	SNPwindows = parse_SNP()

	# extract the haplotypes from the SAM file
	parse_sam(SNPwindows)

if __name__ == '__main__':
	main()
