#!/usr/bin/env python

#_______________________________________________________#

import os
import shutil
import glob
import sys
import argparse
from optparse import OptionParser
from argparse import RawTextHelpFormatter

from sequana import bedtools
from sequana.reporting import report_mapping
from sequana.reporting import report_chromosome
from sequana.reporting import report_main
from sequana import logger
# from sequana import sequana_debug_level
from sequana import GenomeCov

from easydev import shellcmd, mkdirs
from easydev.console import purple

from pylab import show, figure, savefig

#_______________________________________________________#

def getArgs():
	parser = argparse.ArgumentParser(description="")
	parser.add_argument('-i',dest="input",type=str,required=True,help="Input file (BED)")
	parser.add_argument('-f',dest="window_size",type=int,default=4000,help="")
	parser.add_argument('-y',dest="ylim",type=int,default=100,help="ylim")
	parser.add_argument('-l',dest="length",type=int,default=20,help="length (inch)")
	parser.add_argument('-w',dest="width",type=int,default=10,help="width (inch)")
	parser.add_argument('-o',dest="outname",type=str,help="")
	
	args = parser.parse_args()
	
	return args

def main(options):
	
	if options.input.endswith(".bam"):
		bedfile = options.input.replace(".bam", ".bed")
		shellcmd("bedtools genomecov -d -ibam %s > %s" % (options.input, bedfile))
	
	elif options.input.endswith(".bed"):
		bedfile = options.input
	else:
		raise ValueError("Input file must be a BAM or BED file")
	
	###
	
	gc = GenomeCov(bedfile, -4, 4, 0.5, 0.5)
	if len(gc.chr_list) == 1:
		chrom = gc.chr_list[0]
		run_analysis(gc, chrom, 1, options)
	else:
		raise ValueError("Error, more than 1 chr in the BED file")

def run_analysis(gc, chrom, chrom_index, options):
	# chrom.running_median(n=options.window_size, circular=options.circular)
	chrom.running_median(n=options.window_size)
	
	stats = chrom.get_stats(output="dataframe")
	stats.set_index("name", inplace=True)
	
	DOC = stats.ix['DOC'].Value
	if DOC < 8:
		k = 1
	else:
		k = 2
	
	chrom.compute_zscore(k=k)
	
	# Let us save the thresholds first and then change it to compute centralness
	thresholds = chrom.thresholds.copy()
	
	chrom.thresholds.low = -3
	chrom.thresholds.high = 3
	c3 = chrom.get_centralness()
	
	chrom.thresholds.low = -4
	chrom.thresholds.high = 4
	c4 = chrom.get_centralness()
	chrom.thresholds = thresholds.copy()   # Get back to the original values
	
	figure(1)
	# chrom.plot_coverage()
	plot_coverage_perso(chrom, options)

def plot_coverage_perso(self, options, fontsize=16, 
		rm_lw=1, rm_color="#0099cc", rm_label="Running median", 
		th_lw=1, th_color="r", th_ls="--", main_color="k", main_lw=1, 
		main_kwargs={}):

	""" Plot coverage as a function of base position.

	:param filename:

	In addition, the running median and coverage confidence corresponding to
	the lower and upper  zscore thresholds

	.. note:: uses the thresholds attribute.
	"""
	
	import pylab
	# z = (X/rm - \mu ) / sigma

	high_zcov = (self.thresholds.high * self.best_gaussian["sigma"] +
			self.best_gaussian["mu"]) * self.df["rm"]
	low_zcov = (self.thresholds.low * self.best_gaussian["sigma"] +
			self.best_gaussian["mu"]) * self.df["rm"]

	pylab.clf()
	
	ax = pylab.gca()
	ax.set_facecolor('#eeeeee')
	
	pylab.figure(figsize=(options.length, options.width))
	
	pylab.xlim(0,self.df["pos"].iloc[-1])
	axes = []
	labels = []
	# p1, = pylab.plot(self.df["cov"], color=main_color, label="Coverage",
	# 		linewidth=main_lw, **main_kwargs)
	# axes.append(p1)
	# labels.append("Coverage")
	if rm_lw>0:
		p2, = pylab.plot(self.df["rm"],
				color=rm_color,
				linewidth=rm_lw,
				label=rm_label)
		axes.append(p2)
		labels.append(rm_label)

	pylab.legend(axes, labels, loc="best")
	pylab.xlabel("Position", fontsize=fontsize)
	pylab.ylabel("Per-base coverage", fontsize=fontsize)
	pylab.grid(True)
	pylab.ylim(0, options.ylim)
	try:
		pylab.tight_layout()
	except:
		pass
	
	pylab.savefig(options.outname+'.pdf', dpi=400, orientation='landscape', format="pdf")
		
if __name__ == "__main__":
   args = getArgs()
   main(args)
   