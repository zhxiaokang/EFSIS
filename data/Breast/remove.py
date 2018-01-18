#!/usr/bin/env python

# Copyright (C) 2015-2016 Nicola Lazzarini
# School of Computing Science, Newcastle University

# This file is part of the RGIFE heuristic.

# RGIFE is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

# RGIFE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details:
# http://www.gnu.org/licenses/gpl.html

from __future__ import division
import numpy as np
import subprocess
import datasetParser

def set_seed(seed):
	np.random.seed(seed)
	return

def filter_attributes_fold(path_file, path_total, indexes, atts_to_remove, iteration):
	#remove the the attributes (being tested) from each fold
	data = False
	matrix = []
	with open (path_file + "tmp.arff","wt") as output:
		output.write("@relation it"+str(iteration)+"\n")
		for line in open (path_total):
			if line in ['\n', '\r\n', ' ']:
				continue
			if line.startswith("@att"):
				att_name = line.split()[1]
				if att_name in atts_to_remove:
					continue
				output.write(line)
			if line.startswith("@data"):
				output.write("@data\n")
				data = True
				continue
			if data:
				matrix.append(np.array (line.split(",")))
		filtered = np.delete (matrix,indexes,1)

		for row in filtered:
			output.write(",".join(row))
	
	subprocess.call(["mv",path_file + "tmp.arff", path_total])
	return

def filter_whole_dataset(path_file, new_name, iteration, binarised_attribute_mapping, dataset_name):
	#remove the the attributes (being tested) from the whole dataset
	atts_to_keep = []
	for line in open (path_file):
		if line in ['\n', '\r\n', ' ']:
			continue
		if line.lower().startswith("@att"):
			attribute = line.split()[1]
			if attribute in binarised_attribute_mapping:
				att_to_keep = binarised_attribute_mapping[attribute]
			else:
				att_to_keep = attribute
			if att_to_keep not in atts_to_keep:
				atts_to_keep.append(att_to_keep)
	data = False
	matrix = []
	indexes = []
	index = 0
	with open (new_name,"wt") as output:
		output.write("@relation it"+str(iteration)+"\n")
		for line in open (dataset_name):
			if line in ['\n', '\r\n', ' ']:
				continue
			if line.lower().startswith("@att"):
				att_name = line.split()[1]
				if att_name not in atts_to_keep:
					indexes.append(index)
					index += 1
					continue
				index += 1
				output.write(line)
			if line.lower().startswith("@data"):
				output.write("@data\n")
				data = True
				continue
			if data:
				matrix.append(np.array (line.split(",")))
		filtered = np.delete (matrix,indexes,1)
		for row in filtered:
			output.write(",".join(row))
	return
	
def remove_attributes(block_size, actual_iteration, reference_iteration, starting_index, white_list, path_results, name, n_folds ,dataset_name, binarised_attribute_mapping, pre_processing):
	#identify the attributes that needs to be removed based on the current rank
	#create the folder for the new iteration and copy the dataset from the reference iteration folder
	args = ["mkdir","-p", path_results + "it"+str(actual_iteration)]
	subprocess.call(args)
	args = ["cp", path_results + "it"+str(reference_iteration)+"/"+name+"_it"+str(reference_iteration)+".arff", path_results + "it"+str(actual_iteration)+"/"+name+"_it"+str(actual_iteration)+".arff"]
	subprocess.call(args)
	#copy the Training and the Test set from the previous iteration
	for fold in range(0,n_folds):
		subprocess.call("cp " + path_results + "it"+str(reference_iteration)+"/TrainFold"+str(fold) + " " + path_results + "it"+str(actual_iteration)+"/", shell = True)
		subprocess.call("cp " + path_results + "it"+str(reference_iteration)+"/TestFold"+str(fold) + " " + path_results + "it"+str(actual_iteration)+"/", shell = True)
		
	#avoid to copy the original training file TO REMOVE 
	#subprocess.call("rm " + path_results + "it"+str(actual_iteration)+"/"+name+"_it"+str(reference_iteration)+".arff", shell = True)
	ranked = []
	for line in open(path_results+ "it"+str(reference_iteration)+"/attribute_score.txt"):
		attribute = line.split()[0]
		if attribute not in white_list:
			ranked.append(attribute)
	
	print "Atts of reference dataset {0}".format(len(ranked)+len(white_list))
	num_to_remove = block_size
	if num_to_remove >= len(ranked):
		return -1
	#if there are not enough elements for the last block
	if starting_index + num_to_remove >= len (ranked):
		to_remove = ranked[starting_index:]
	else:
		to_remove = ranked[starting_index:(starting_index + num_to_remove)]
	
	#Get the indexes of the attributes to remove
	indexes = datasetParser.get_index_atts (path_results + "it"+str(actual_iteration)+"/TrainFold0")
	indexes_to_remove = []
	atts_to_remove = {}
	for att in to_remove:
		indexes_to_remove.append(indexes[att])
		atts_to_remove[att] = True
	#Remove for each fold
	for fold in range(0,n_folds):
		filter_attributes_fold(path_results+"it"+str(actual_iteration)+"/", path_results+"it"+str(actual_iteration)+"/TrainFold"+str(fold), indexes_to_remove, atts_to_remove, actual_iteration)
		filter_attributes_fold(path_results+"it"+str(actual_iteration)+"/", path_results+"it"+str(actual_iteration)+"/TestFold"+str(fold), indexes_to_remove, atts_to_remove, actual_iteration)
	#remove files from the whole dataset
	if pre_processing == 0:
		filter_attributes_fold(path_results+"it"+str(actual_iteration)+"/", path_results + "it"+str(actual_iteration)+"/"+name + "_it"+str(actual_iteration)+".arff", indexes_to_remove, atts_to_remove, actual_iteration)
	else:
		#check which binarised atts are kept and filter the original (no binarise atts) dataset according to them
		filter_whole_dataset(path_results+"it"+str(actual_iteration)+"/TrainFold0", path_results + "it"+str(actual_iteration)+"/"+name + "_it"+str(actual_iteration)+".arff", actual_iteration, binarised_attribute_mapping, dataset_name)
	return starting_index + num_to_remove
