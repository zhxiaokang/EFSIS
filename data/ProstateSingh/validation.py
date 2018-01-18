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
import preprocessing
import datasetParser
import remove
import numpy as np
import scipy.stats as ss
import random
from sklearn.metrics.pairwise import manhattan_distances


def set_seed(seed):
	random.seed(seed)
	np.random.seed(seed)
	return
	
def cvs_folds(name, path_file, categorical_attributes, binarised_attribute_mapping, n_folds, db_scv):
	#generate the train/test folds using either: SCV or DB_SCV
	if categorical_attributes.lower() == "yes":
		samples_for_distance = []
		preprocessing.binarise_attributes(name, binarised_attribute_mapping)
		binarised_name = name.split(".")[0] + "_Binarised"
		data = False
		for line in open (binarised_name):
			if line.lower().startswith("@data"):
				data = True
				continue
			if data:
				sample = [float(i) if i != '?' else 0 for i in line.split(",")[:-1]]
				samples_for_distance.append(sample)
		samples_for_distance = np.array(samples_for_distance)
	class_mapping = {}
	class_samples = {}
	attributes = []
	samples = []
	labels = []
	data_matrix = []

	previous_line = ""
	header = ""
	data = False
	index = 0
	#read the original file and get attributes info (header to be copied in each fold)
	for line in open (name):
		if line in ['\n', '\r\n', ' ']:
			continue
		if line.startswith("%"):
			continue
		if line.lower().startswith("@rel"):
			header = header+line
			previous_line = line
			continue
		if line.lower().startswith("@att"):
			attribute = line.split()[1]
			attributes.append(attribute)
			header = header+line
			previous_line = line
			continue
		if line.lower().startswith("@data"):
			classes = [l.strip(" {}'\n\r") for l in previous_line.split("{")[1].strip("}").split(",")]
			for i,label in enumerate(classes):
				class_mapping[label] = i
				class_samples[i] = []
			data = True
			continue
		if data:
			label = line.split(",")[-1].strip(" ' \n\r")
			labels.append(class_mapping[label])
			class_samples[class_mapping[label]].append(index)
			index += 1                   
			sample = [i if i != '?' else 0 for i in line.split(",")[:-1]]
			samples.append(sample)
			data_matrix.append(line)
	
	samples = np.array(samples)
	labels = np.array(samples)

	folds = {}
	
	for fold in range (0, n_folds):
		folds[fold] = []
		
	if (db_scv):
		if categorical_attributes.lower() == "no":
			samples_for_distance = samples 
		neighbours = {}
		#rows that belong to a specific class
		n_samples = samples.shape[0]
		dist = manhattan_distances(samples_for_distance,samples_for_distance)
		for i,sample in enumerate(samples_for_distance):
			neighbours[i] = [0]*(n_samples)
			ranked =  ss.rankdata(dist[i,:],method ='ordinal')
			for j,element in enumerate(ranked):
			#position = ranked[element]
				neighbours[i][int(element)-1] = j

		#db_scv
		for class_label in class_samples:
			sample_subset = class_samples[class_label]
			n = len(class_samples[class_label])//n_folds
			residuals = len(class_samples[class_label]) - (n * n_folds)
			e = random.choice(sample_subset)
			i = 0
			while (len(sample_subset) > residuals):
				folds[i].append(e)
				sample_subset.remove(e)
				#if there are no more instances
				if len(sample_subset) == 0:
					break
				i = (i+1) % n_folds
				#get the closest samples of e
				index = 1
				closest = neighbours[e][index]
				while (closest not in sample_subset):
					index += 1
					closest = neighbours[e][index]
				e = closest
			
			#folds that will have the extra sample
			fold_extra_samples = random.sample(range(0,n_folds),len(sample_subset))
			for fold in fold_extra_samples:
				e = random.choice(sample_subset)
				folds[fold].append(e)
				sample_subset.remove(e)
	else:
		#scv
		for class_label in class_samples:
			#divide the sample of the class within each fold
			sample_subset = class_samples[class_label]
			n = len(class_samples[class_label])//n_folds
			for fold in range (0,n_folds):
				for i in range(0,n):
					e = random.choice(sample_subset)
					folds[fold].append(e)
					sample_subset.remove(e)

			#folds that will have the extra sample
			fold_extra_samples = random.sample(range(0,n_folds),len(sample_subset))
			for fold in fold_extra_samples:
				e = random.choice(sample_subset)
				folds[fold].append(e)
				sample_subset.remove(e)
	
	#write the train/test fold files
	for fold in range (0,n_folds):
		test_indexes = folds[fold]
		with open (path_file+"/TrainFold"+str(fold),"wt") as output_train:
			with open (path_file+"/TestFold"+str(fold),"wt") as output_test:

				output_train.write(header)
				output_test.write(header)
				
				output_train.write("@data\n")
				output_test.write("@data\n")
				
				for sample in range(0,len(data_matrix)):
					if sample in test_indexes:
						output_test.write(data_matrix[sample])
					else:
						output_train.write(data_matrix[sample])
	

def loocv_folds(name, path_file):
	#generate the train/test folds using LOOCV
	header = ""
	data_matrix = []
	data = False

	for line in open (name):
		if line.startswith("@"):
			header = header+line
		if line.lower().startswith("@data"):
			data = True
			continue
		if data:
			data_matrix.append(line)

	n_samples = len(data_matrix)

	for fold in range (0,n_samples):
		with open (path_file + "/TrainFold"+str(fold),"wt") as output:
			output.write(header)
			for i,sample in enumerate(data_matrix):
				if i != fold:
					output.write(sample)
		output.close()
		with open (path_file + "/TestFold"+str(fold),"wt") as output:
			output.write(header)
			output.write(data_matrix[fold])
		output.close()

def fix_new_folds (iteration, current_attributes, path_results, n_folds):
	#when having multiple REPETITONS the new folds, after being generate contains all the binarised variables from the current set of variables
	#current_attribute -> current set of attributes
	#remove all the binarised variables that are not in current_attributes
	indexes = datasetParser.get_index_atts (path_results + "it"+str(iteration)+"/TrainFold0")
	atts_to_remove = {}
	indx_to_remove = []
	for line in open(path_results + "it"+str(iteration)+"/TrainFold0"):
		if line.startswith("@att"):
			att = line.split()[1]
			if att not in current_attributes:
				atts_to_remove[att] = True
				indx_to_remove.append(indexes[att])
	#Remove for each fold
	for fold in range(0,n_folds):
		remove.filter_attributes_fold(path_results+"it"+str(iteration)+"/", path_results+"it"+str(iteration)+"/TrainFold"+str(fold), indx_to_remove, atts_to_remove, iteration)
		remove.filter_attributes_fold(path_results+"it"+str(iteration)+"/", path_results+"it"+str(iteration)+"/TestFold"+str(fold), indx_to_remove, atts_to_remove, iteration)
	return	
