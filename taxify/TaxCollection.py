#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from MiscFunctions import n50, keyWithMaxVal
import numpy
import time
import subprocess
import commands
import re
import sys

class Tax():
	'''
	Tax contain the information of a contig in an assembly. 
	'''

	def __init__(self, header, seq):
		self.name = header
		self.seq = seq.upper()
		self.length = len(seq)
		self.corrected_length = len(self.seq) - self.seq.count('N') 
		self.gc = float((self.seq.count('G') + self.seq.count('C') ) / self.corrected_length ) if self.corrected_length > 0 else 0.0
		self.covs = dict()
		self.tax = dict() 

class TaxCollection():
	'''
	The TaxCollection class contains all methods needed 
	for parsing the input files and doing the necessary computations 
	for arriving at the final results. 
	'''
	def __init__(self, target_ranks, rank_level):
		self.contigs = dict()
		self.index = list()
		self.outfile = str()
		self.stats = dict()
		self.cov_libs = list()
		self.blast_libs = list()
		self.target_ranks = target_ranks
		self.rank_level = rank_level
		self.blast_order = list()

	def addTax(self, tax):
		if not tax.name in self.contigs: 
			self.contigs[tax.name] = tax
		else: 
			sys.exit("[ERROR] - Sequence header {} occurs more than once".format(tax.name))
		self.index.append(tax.name)

	def getTaxsFromAssembly(self, assembly_file, assembly_type, exclude_assembly_cov):
		print "[STATUS] - Parsing assembly %s" % (assembly_file)
		header, seq = '', ''
		with open(assembly_file) as fh:
			for line in fh:
				line_data = line.rstrip("\n")
				if line.startswith('>'):
					if (seq): 
						tax = Tax(header, seq) 
						self.addTax(tax)
						seq = ''
						if assembly_type == 'unknown' or exclude_assembly_cov:
							pass
						else:
							cov = float(self.parseCovFromHeader(header, assembly_type))
							self.addTaxCov(header, assembly_type, cov)
					header = line.rstrip("\n").lstrip(">")
				else:
					seq += line.rstrip("\n").upper() 
			tax = Tax(header, seq) 
			self.addTax(tax)
			if assembly_type == 'unknown' or exclude_assembly_cov:
				pass
			else:
				cov = float(self.parseCovFromHeader(header, assembly_type))
				self.addTaxCov(header, assembly_type, cov)
		if not assembly_type == 'unknown' and not exclude_assembly_cov:
			self.cov_libs.append(assembly_type)

	def getTaxForTaxs(self, blast_files, taxdb):
		'''
		- Parses NCBI nodes, names and all BLAST files
		- gets taxonomy for each BLAST hit taxid (getTaxonomy)
		- for each BLAST file it sums up the bitscore by rank_level (e.g. phylum) of hit taxonomy and stores them in tax.tax[blast_lib]
		- then it goes through all taxs again and sets those without hits to 'no-hit'
		'''
		nodes_dict = self.parse_taxdb_nodes(taxdb['nodes'])
		names_dict = self.parse_taxdb_names(taxdb['names'])

		for blast_lib, blast_file in blast_files.items():
			self.blast_libs.append(blast_lib)
			with open(blast_file) as fh:
				for line in fh:
					line_data = line.rstrip("\n").split("\t")
					blast_line_re = re.compile(r"^(\S+)\t(\S+)\t\s*(\S+)") # turns out that blastn output is not tab delimited but tab/(+space) delimited
					match = blast_line_re.search(line)
					if match:
						try:
							qseqid, taxid, bitscore = match.group(1), int(match.group(2).split(";")[0]), float(match.group(3)) # if more than one taxid is found ... first one will be used
						except ValueError, e:
							sys.exit("[ERROR] : BLAST output does not seem to be in the right format ('6 qseqid staxids bitscore ... ')\n" + line)
						
						if qseqid in self.contigs: 
							taxonomy = {}
							taxonomy = self.getTaxonomy(taxid, taxonomy, nodes_dict, names_dict) # infers taxonomy based on 
							if not blast_lib in self.contigs[qseqid].tax:
								self.contigs[qseqid].tax[blast_lib] = {}
							taxonomic_group = taxonomy[self.rank_level]
							self.contigs[qseqid].tax[blast_lib][taxonomic_group] = self.contigs[qseqid].tax[blast_lib].get(taxonomic_group, 0) + bitscore							
						else:
							print ("[WARN] - {} in {} does not seem to be part of the assembly.".format(qseqid, blast_file))
					else:
						sys.exit("[ERROR] - BLAST results in {} do not seem to be in the right format. Please run BLAST with the following option : -outfmt '6 qseqid staxids bitscore'".format(blast_file))

			for contig_name in self.contigs:
				if not blast_lib in self.contigs[contig_name].tax:
					self.contigs[contig_name].tax[blast_lib] = {} 
					self.contigs[contig_name].tax[blast_lib]['no-hit'] = 0

	def getTaxonomy(self, taxid, taxonomy, nodes_dict, names_dict):
		'''
		gets target_ranks from nodes/names based on taxid and returns taxonomy
		'''
		if taxid in nodes_dict: 
			rank = nodes_dict[int(taxid)]['rank']
			parent = int(nodes_dict[int(taxid)]['parent'])
			if taxid == 1 or rank == 'superkingdom':
				# finish if taxid == 1 (root) or rank == superkingdom
				taxonomy[rank] = names_dict[int(taxid)] 
				for rank in self.target_ranks:
					# set rank to undef if ranks in target ranks was not populated 
					taxonomy[rank] = taxonomy.get(rank, "undef")
				return taxonomy
			else:
				if rank in self.target_ranks:
					taxonomy[rank] = names_dict[int(taxid)] 
				self.getTaxonomy(parent, taxonomy, nodes_dict, names_dict)
		else:
			# if taxid is not in nodes_dict then print warning and set all ranks to undef
			print "[WARN] - Taxid %s not found in TAXDB. This is probably because you have BLASTed against a newer version of NCBI nt. You should update your TAXDB." % taxid
			taxonomy = {rank: "undef" for rank in self.target_ranks}
		return taxonomy

	def parse_taxdb_names(self, infile):
		print "[STATUS] - Parsing names.dmp from NCBI taxdb" 
		names_dict = {}
		with open(infile) as fh:
			for line in fh:
				fields = line.split("\t")
				if fields[6] == "scientific name":
					names_dict[int(fields[0])] = fields[2]
		return names_dict

	def parse_taxdb_nodes(self, infile):
		print "[STATUS] - Parsing nodes.dmp from NCBI taxdb" 
		nodes_dict = {}
		with open(infile) as fh:
			for line in fh:
				fields = line.split("\t")
				nodes_dict[int(fields[0])]={ 'parent' : fields[2] , 'rank' : fields[4] }
		return nodes_dict

	def getConsensusTaxForTaxs(self, taxrule, blast_order):
		''' 
		- Based on taxrule ("A" or "B") and the blast_order (list in order in which blast files where specified) 
		it calculates the consensus taxonomy for each tax 
		- if taxrule == A:
			- it puts all taxonomic groups in a dict with their summed scores as values
			- if a taxonomic group occurs in hits of more than one BLAST file, the highest score is used
		- if taxrule == B:
			- taxonomic groups are put in the dict with their summed scores as values IF they come from the first BLAST file
			- If there was no hit then take the taxonomic groups from the next one  	
		- The highest scoring taxonomic group is selected as consensus taxonomy for each tax
		'''
		for contig_name in self.contigs:
			dict_for_tax_merging = {}
			for blast_lib in blast_order:
				for tax, score in sorted(self.contigs[contig_name].tax[blast_lib].items(), key=lambda x: x[1], reverse=True):
					# loops through tax/score with decreasing score
					if taxrule == 'A':
						if not tax in dict_for_tax_merging:
							dict_for_tax_merging[tax] = score
						else:
							if score > dict_for_tax_merging[tax]:
								dict_for_tax_merging[tax] = score 
					elif taxrule == 'B':
						if blast_lib == blast_order[0]:
							# First blast_lib
							dict_for_tax_merging[tax] = score
						else:
							if len(dict_for_tax_merging) <= 1 and ('no-hit' in dict_for_tax_merging):
								dict_for_tax_merging[tax] = score	

			tax = keyWithMaxVal(dict_for_tax_merging)
			self.contigs[contig_name].tax['tax'] = {}
			self.contigs[contig_name].tax['tax'][tax]=dict_for_tax_merging[tax]
		self.blast_libs.append('tax')

	def getStats(self):
		'''
		Calculates all different kind of stats for each taxonomic group based on each BLAST file and the consensus... 
		'''
		self.stats['count'] = {}
		self.stats['span']= {}
		self.stats['n50']= {}
		self.stats['lengths']= {}
		self.stats['gc']= {}
		self.stats['total_count'] = 0
		self.stats['total_span'] = 0
		self.stats['total_n50'] = 0
		self.stats['total_lengths'] = []
		self.stats['total_gc'] = {'raw' : [], 'mean' : 0.0, 'stdev' : 0.0}

		for contig_name in self.contigs:
			tax = self.contigs[contig_name]
			self.stats['total_count'] += 1
			self.stats['total_span'] += tax.length
			self.stats['total_lengths'].append(tax.length)
			self.stats['total_gc']['raw'].append(tax.gc)

			for blast_lib in self.blast_libs:
				
				bestTax = keyWithMaxVal(tax.tax[blast_lib])
				if not blast_lib in self.stats['count']:
					self.stats['count'][blast_lib] = {}
					self.stats['span'][blast_lib] = {}
					self.stats['lengths'][blast_lib] = {}
					self.stats['gc'][blast_lib] = {}
				self.stats['count'][blast_lib][bestTax] = self.stats['count'][blast_lib].get(bestTax, 0) + 1	
				self.stats['span'][blast_lib][bestTax] = self.stats['span'][blast_lib].get(bestTax, 0) + tax.length

				if not bestTax in self.stats['gc'][blast_lib]:
					self.stats['gc'][blast_lib][bestTax] = {'raw' : [], 'mean' : 0.0, 'stdev' : 0.0}
					self.stats['lengths'][blast_lib][bestTax] = []	
				self.stats['gc'][blast_lib][bestTax]['raw'].append(tax.gc) 
				self.stats['lengths'][blast_lib][bestTax].append(tax.length)

		for blast_lib in self.blast_libs:
			# calculate N50
			for tax, list_of_lengths in self.stats['lengths'][blast_lib].items():
				if not blast_lib in self.stats['n50']:
					self.stats['n50'][blast_lib] = {}
				self.stats['n50'][blast_lib][tax] = n50(list_of_lengths)
			self.stats['total_n50'] = n50(self.stats['total_lengths'])

			# calculate total gc mean/stdev
			for tax in self.stats['gc'][blast_lib]:
				self.stats['gc'][blast_lib][tax]['mean'] = "{0:.2f}".format(numpy.mean(self.stats['gc'][blast_lib][tax]['raw']))
				self.stats['gc'][blast_lib][tax]['stdev'] = "{0:.2f}".format(numpy.std(self.stats['gc'][blast_lib][tax]['raw']))
		
		self.stats['total_gc']['mean'] = "{0:.2f}".format(numpy.mean(self.stats['total_gc']['raw']))
		self.stats['total_gc']['stdev'] = "{0:.2f}".format(numpy.std(self.stats['total_gc']['raw']))

	def writeOutput(self, version):
		
		'''
		Writes outputfiles:
		- stats.txt which contains tables with the stats calculated through getStats()
		- taxs.txt which contains the taxs
		'''

		header = '# maketaxs.py v{}\n# {} {}\n# {}\n'.format(version, time.strftime("%Y-%m-%d"), time.strftime("%H:%M:%S"), " ".join(sys.argv))
		
		# Writing of stats.txt
		stats_fh = open(self.outfiles['stats'], 'w')
		stats_fh.write(header)
		stats_string = ''
		for blast_lib in self.blast_libs:
			stats_string += "\tTAX:{:>10}\t{:>10}\t{:>10}\t{:>10}\t{:^10}".format(blast_lib, "contigs", "span", "N50", "GC")
			stats_string += "\n\t" + (("-" * 10) + "\t") * (5 + len(self.cov_libs)) + "\n"
			for tax, score in sorted(self.stats['span'][blast_lib].items(), key=lambda x: x[1], reverse = True):
				stats_string += "\t{:>10}\t{:>10}\t{:>10}\t{:>10}\t{:>10}".format(tax, self.stats['count'][blast_lib][tax], self.stats['span'][blast_lib][tax], self.stats['n50'][blast_lib][tax],  str(self.stats['gc'][blast_lib][tax]['mean']) + " SD:" + str(self.stats['gc'][blast_lib][tax]['stdev']))
				stats_string += "\n"
			stats_string += "\t{:>10}\t{:>10}\t{:>10}\t{:>10}\t{:>10}".format("Total", self.stats['total_count'], self.stats['total_span'], self.stats['total_n50'],  str(self.stats['total_gc']['mean']) + " SD:" + str(self.stats['total_gc']['stdev']))
			stats_string += "\n\n"
		stats_fh.write(stats_string)
		stats_fh.close()

		# Writing of taxs.txt
		taxs_fh = open(self.outfiles['taxs'], 'w')
		taxs_fh.write(header)
		taxs_fh.write("# contig_id\tlength\tgc\ttaxonomy\n")
		taxs_string = ''
		for contig_name in self.index:
			tax = self.contigs[contig_name]
			taxs_string += "{}\t{}\t{:.3f}".format(tax.name, tax.length, tax.gc)
			tax_string = '\t'
			for blast_lib, tax in tax.tax.items():
				tax_string += "{}=".format(blast_lib)
				for phylum, score in sorted(tax.items(), key=lambda x: x[1], reverse = True):
					tax_string += "{}:{},".format(phylum, score)
				tax_string = tax_string[0:-1] + ";"
			taxs_string += tax_string[0:-1]
			taxs_string += "\n"
		taxs_fh.write(taxs_string)
		taxs_fh.close()