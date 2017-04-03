#!/usr/bin/env python
# -*- coding: utf-8 -*-

def keyWithMaxVal(d):
	v=list(d.values())
	k=list(d.keys())
	return k[v.index(max(v))]

def n50(list_of_lengths):
	total_span = 0
	sorted_list_of_lengths=sorted(list_of_lengths, reverse=True)
	for contig_length in sorted_list_of_lengths:
		total_span += contig_length
	teoN50 = total_span/2.0
	running_sum = 0
	N50 = 0
	for contig_length in sorted_list_of_lengths:
		running_sum += contig_length
		if teoN50 <= running_sum:
			N50 = contig_length
			break
	return N50
