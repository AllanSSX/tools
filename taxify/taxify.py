#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import sys 
import argparse
import os
import TaxCollection

class InputObject():
	def __init__(self, args):
		self.assembly_type, self.assembly_file = self.getAssembly(args)
		self.blast_files, self.blast_order = self.getDictOfBLASTFiles(args)
		self.taxdb = self.getDictOfTaxDB(args.taxdb)
		self.outfiles = self.getDictOfOutfiles(args.o) 
		self.taxrule = args.taxrule 
		self.rank = self.checkRank(args.rank)

	def checkRank(self, rank):
		if rank in RANKS:
			return rank
		else:
			sys.exit("[ERROR] : \'{}\' is not a valid taxonomic rank. Please select one of the following taxonomic ranks: {}.".format(rank, ", ".join(RANKS)))

	def printParameters(self):
		print "\n\tAssembly file :\n\t\t- {} (type = '{}')".format(self.assembly_file, self.assembly_type)
		print "\tBLAST file(s) :"
		for blast_name, blast_file in self.blast_files.items():
			print "\t\t- {} : {}".format(blast_name, blast_file)
		print "\tTaxDB files :"
		for taxdb_name, taxdb_file in self.taxdb.items():
			print "\t\t- {} : {}".format(taxdb_name, taxdb_file)
		if len(self.blast_order) > 1:
			print "\tTaxonomification rule : \"{}\"".format(self.taxrule)
			if self.taxrule == 'A':
				print "\t\t - Taxid with maximal sum of scores is chosen"
			elif self.taxrule == 'B':
				print "\t\t - Confidence in taxid decreases with order in BLAST libs :\n\t\t {}".format(str(self.blast_order))
			else:
				pass
		print "\tTaxonomy will be inferred for the following taxonomic rank :\n\t\t - \"{}\"".format(self.rank.capitalize())
		print "\tOutfiles :\n\t\t - Tax-file : {}\n\t\t - Stats-file : {}\n".format(self.outfiles['taxs'], self.outfiles['stats']) 

	def getAssembly(self, args):
		assembly = {'unknown' : args.a}
		assembly = {k: v for k, v in assembly.items() if v}	
		try:
			# get assembly type, assembly file and check whether file exists
			[(assembly_type, assembly_file)] = assembly.items()
			if os.path.exists(assembly_file):
				pass
			else:
				sys.exit("[ERROR] : Assembly file {} does not exist.".format(assembly_file))		
		except ValueError, e:
			# throw ValueError if there are too many (>1) elements to unpack from assembly   
			sys.exit("[ERROR] : Please specify an assembly file.")	
		if assembly_type == 'unknown':
			self.exclude_assembly_cov = True
		return assembly_type, assembly_file

	def getDictOfBLASTFiles(self, args):
		blasts = {}
		order = []  
		blast_files = args.blast
		set_of_files = set()
		blast_count = 0
		for blast_file in blast_files:
			if os.path.exists(blast_file):
				if blast_file in set_of_files:
					sys.exit("[ERROR] : BLAST file %s listed more than once." %blast_file)		
				set_of_files.add(blast_file)
				blast_count += 1
				blast_lib = "BLAST_" + str(blast_count)
				blasts[blast_lib] = blast_file
				order.append(blast_lib)
			else:
				sys.exit("[ERROR] : BLAST file %s does not exist." %blast_file)	
		if not blasts:
			sys.exit("[ERROR] : Please specify at least one BLAST file.")	
		return blasts, order

	def getDictOfTaxDB(self, taxdb_dir):
		taxdb = {}
		if not taxdb_dir:
			sys.exit("[ERROR] : Please specify the path to NCBI TaxDB.")
		if os.path.exists(taxdb_dir):
			if not os.path.exists(taxdb_dir + "nodes.dmp"):
				sys.exit("[ERROR] : 'nodes.dmp' from NCBI TaxDB could not be found in %s." %taxdb_dir)
			if not os.path.exists(taxdb_dir + "names.dmp"):
				sys.exit("[ERROR] : 'names.dmp' from NCBI TaxDB could not be found in %s." %taxdb_dir)	
		else:
			sys.exit("[ERROR] : NCBI TaxDB could not be found in %s." %taxdb_dir)
		taxdb['nodes'] = taxdb_dir + "nodes.dmp"
		taxdb['names'] =  taxdb_dir + "names.dmp"
		return taxdb

	def getDictOfOutfiles(self, outprefix):
		outfiles = {}
		if outprefix:
			outfiles['taxs'] = outprefix + ".taxs.txt"
			outfiles['stats'] = outprefix + ".taxs.stats.txt"
		else:
			outfiles['taxs'] = "taxs.txt"
			outfiles['stats'] = "taxs.stats.txt"
		return outfiles

def getInput():
	parser = argparse.ArgumentParser(
		prog='makeblobs.py',
		usage = '%(prog)s -a <ASSEMBLY> -cas <CAS> -blast <BLAST> -taxdb <PATH_TO_TAXDB> -o <OUTPUT> [-h]',
		add_help=True)
	# only ONE assembly file
	parser.add_argument('-a', metavar = 'ASSEMBLY_FASTA', default='', help='Assembly file')
	# multiple BLAST files
	parser.add_argument('-blast', metavar = 'BLAST_FILE', default=[], nargs='+', help='BLAST file') 
	parser.add_argument('-rank', metavar = 'TAX_RANK', default='phylum', help='Select target taxonomic rank (species, genus, order, phylum, superkingdom). Default = phylum') 
	parser.add_argument('-taxrule', metavar = 'A or B', default='A', help='Tax-rule on how to deal with multiple BLAST libs. A : "higher bitscore wins", B : "Decreasing trust in BLAST libs"') 
	parser.add_argument('-taxdb', metavar = 'TAX_DUMP', default='', help='Path to NCBI taxdb (nodes.dmp, names.dmp)') 
	parser.add_argument('-o', metavar = 'OUTPUT_PREFIX', default='', help='Output prefix') 
	# Version number
	parser.add_argument('-v', action='version', version='%(prog)s version 0.1')
	args = parser.parse_args()
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	return InputObject(args)

def parseInput(parameters):
	'''
	1. Create a TaxCollection Object using RANKS, rank_level 
	... technically only rank_level is needed if we agree that we 
	do not need the "whole" taxonomy (kingdom, phylum, order, genus, species)
	I think it is enough to only have one taxonomic rank ... 
	'''
	data = TaxCollection.TaxCollection(RANKS, parameters.rank)
	'''
	2. Save outfile-dict as part of TaxCollection object
	'''
	data.outfiles = parameters.outfiles
	'''
	3. Parse contigs of assembly into Tax objects in TaxCollection object
	'''
	data.getTaxsFromAssembly(parameters.assembly_file, parameters.assembly_type, parameters.exclude_assembly_cov)
	'''
	4. Parse BLAST files into Tax objects in TaxCollection object and group into taxonomic bins by translating taxids to taxonomic group
	'''
	data.getTaxForTaxs(parameters.blast_files, parameters.taxdb)
	'''
	5. Infer a "taxonomy" from the taxonomic bins based on taxrule for each of the blobs
	'''
	data.getConsensusTaxForTaxs(parameters.taxrule, parameters.blast_order)
	'''
	6. Return TaxCollection object
	'''
	return data

if __name__ == "__main__":
	
	__version__ = 0.1

	# ranks are being parsed from nodes/names ... could be only rank_level
	RANKS = ['species', 'genus', 'family', 'order', 'phylum', 'superkingdom']

	# Put all parameters of the run into InputObject
	parameters = getInput()
	
	# Print InputObject so that user gets feedback about whats happening
	parameters.printParameters()
	
	# # Parse files according to InputObject into TaxCollection object 
	data = parseInput(parameters)
	
	# # Do stats on TaxCollection object so that the stats file can be printed
	data.getStats()
	# 
	# # Write Output to the two files *blobplot.txt and *stats.txt
	data.writeOutput(__version__)

