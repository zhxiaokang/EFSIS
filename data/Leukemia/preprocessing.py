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
from collections import Counter
import subprocess
import numpy as np
import random


def set_seed(seed):
	random.seed(seed)
	np.random.seed(seed)
	return
	
def binarise_attributes(data_set, binarised_attribute_mapping):
	#transforms every categorical attribute with N values in N binary (0,1) variables
	new_name = data_set.split(".")[0] + "_Binarised"
	data = False
	att_index = 0
	mapping = {}
	isCategorical = {}
	#read the data set
	with open (new_name, "wt") as output:
		with open (data_set) as f:
			lines = f.readlines()
			for line_index, line in enumerate(lines):
				#skip empty lines
				if line in ['\n', '\r\n', ' ']:
					continue
				if line.startswith("%"):
					continue
				if line.lower().startswith("@attribute"):
					#get the next line and check if it's the class attribute
					j = line_index + 1
					while(True):
						next_line = lines[j]
						if next_line not in ['\n', '\r\n', ' ']:
							break
						j += 1	
					#do not binarise the class attribute
					if next_line.lower().startswith("@data"):
						isCategorical[att_index] = False
						output.write(line)
						continue
					attribute_name = line.split()[1]
					attribute_type = line.lower().split()[2]
					if attribute_type == "numeric" or attribute_type == "real"  or attribute_type == "string" or attribute_name.lower() == "id":
						output.write(line)
						isCategorical[att_index] = False
					else:
						isCategorical[att_index] = True
						#categorical attribute
						values = line.split("{")[1].strip("\n\r }").split(",")
						mapping[att_index] = {}
						if len(values) > 2:
							new_values = [0]*len(values)
							mapping[att_index]["?"] = new_values
							for i,value in enumerate(values):
								value = value.strip(" ")
								new_values = [0]*len(values)
								new_values[i] = 1
								mapping[att_index][value] = new_values
								output.write("@attribute is_{0}_{1} {{0,1}}\n".format(attribute_name, value))
								binarised_attribute_mapping["is_{0}_{1}".format(attribute_name, value)] = attribute_name
						else:
							output.write("@attribute {0} {{0,1}}\n".format(attribute_name))	
							mapping[att_index]["?"] = [0]
							for i,value in enumerate(values):
								value = value.strip(" ")
								mapping[att_index][value] = [i]
					att_index += 1
					continue				
				if line.lower().startswith("@"):
					output.write(line)
					if line.lower().startswith("@data"):
						data = True
					continue
				if data:
					new_line = ""
					values = line.split(",")
					for i, value in enumerate(values):
						value = value.strip("\n")
						if isCategorical[i]:
							new_line += ",".join(str(x) for x in mapping[i][value]) + ","
						else:
							new_line += value + ","
					output.write(new_line[:-1] + "\n")
	return


def preprocess_dataset(pre_processing, path_results, actual_iteration, dataset_name, n_folds, binarised_attribute_mapping, isCategorical, relation_name):
	#pre-process the dataset by applying: mv imputation and bianarisation
	
	#copy the original files
	subprocess.call(["mkdir","-p", path_results + "it"+str(actual_iteration) + "/OriginalFiles"])
	subprocess.call(["cp", dataset_name, path_results+"it"+str(actual_iteration)+"/OriginalFiles/"])
	
	for fold in range (0,n_folds):
		subprocess.call(["cp", path_results+"it"+str(actual_iteration)+"/TrainFold"+str(fold), path_results+"it"+str(actual_iteration)+"/OriginalFiles/"])
		subprocess.call(["cp", path_results+"it"+str(actual_iteration)+"/TestFold"+str(fold), path_results+"it"+str(actual_iteration)+"/OriginalFiles/"])
		
	if pre_processing == 1 or pre_processing == 2:
		for fold in range (0,n_folds):
			compute_missing_values_mean(path_results+"it"+str(actual_iteration)+"/TrainFold"+str(fold),path_results+"it"+str(actual_iteration)+"/TestFold"+str(fold), isCategorical, relation_name)
		if pre_processing == 1:
			#missing values and cat.attributes
			for fold in range (0,n_folds):
				binarise_attributes(path_results+"it"+str(actual_iteration)+"/TrainFold"+str(fold)+"_NoMissing", binarised_attribute_mapping)
				binarise_attributes(path_results+"it"+str(actual_iteration)+"/TestFold"+str(fold)+"_NoMissing", binarised_attribute_mapping)
				subprocess.call(["mv",  path_results+"it"+str(actual_iteration)+"/TrainFold"+str(fold)+"_NoMissing_Binarised", path_results+"it"+str(actual_iteration)+"/TrainFold"+str(fold)])
				subprocess.call(["mv",  path_results+"it"+str(actual_iteration)+"/TestFold"+str(fold)+"_NoMissing_Binarised", path_results+"it"+str(actual_iteration)+"/TestFold"+str(fold)])
				subprocess.call(["rm",  path_results+"it"+str(actual_iteration)+"/TrainFold"+str(fold)+"_NoMissing"])
				subprocess.call(["rm",  path_results+"it"+str(actual_iteration)+"/TestFold"+str(fold)+"_NoMissing"])
		else:
			#only missing values
			for fold in range (0,n_folds):
				subprocess.call(["mv",  path_results+"it"+str(actual_iteration)+"/TrainFold"+str(fold)+"_NoMissing", path_results+"it"+str(actual_iteration)+"/TrainFold"+str(fold)])
				subprocess.call(["mv",  path_results+"it"+str(actual_iteration)+"/TestFold"+str(fold)+"_NoMissing", path_results+"it"+str(actual_iteration)+"/TestFold"+str(fold)])
	else:
		for fold in range (0,n_folds):
			binarise_attributes(path_results+"it"+str(actual_iteration)+"/TrainFold"+str(fold), binarised_attribute_mapping)
			binarise_attributes(path_results+"it"+str(actual_iteration)+"/TestFold"+str(fold), binarised_attribute_mapping)
			subprocess.call(["mv",  path_results+"it"+str(actual_iteration)+"/TrainFold"+str(fold)+"_Binarised", path_results+"it"+str(actual_iteration)+"/TrainFold"+str(fold)])
			subprocess.call(["mv",  path_results+"it"+str(actual_iteration)+"/TestFold"+str(fold)+"_Binarised", path_results+"it"+str(actual_iteration)+"/TestFold"+str(fold)])

	return

def compute_missing_values_mean (training_set, test_set, isCategorical, relation_name):
	#impute the missing values by using the mean values across samples
	
	#read training file
	header_training = ""
	training_matrix = []
	training_labels = []
	for line in open (training_set):
		if line in ['\n', '\r\n', ' '] or line.lower().startswith("@relation"):
			continue
		if line.startswith("@"):
			header_training += line
		else:
			line = line.strip("\n\r ")
			line.replace(" ", "")
			sample = line.split(",")[:-1]
			sample_label = line.split(",")[-1]
			training_matrix.append(sample)
			training_labels.append(sample_label)
	
	#read test file
	header_test = ""
	test_matrix = []
	test_labels = []
	for line in open (test_set):
		if line in ['\n', '\r\n', ' ']:
			continue
		if line.startswith("%"):
			continue
		if line.startswith("@"):
			header_test += line
		else:
			line = line.strip("\n\r ")
			line.replace(" ", "")
			sample = line.split(",")[:-1]
			sample_label = line.split(",")[-1]
			test_matrix.append(sample)
			test_labels.append(sample_label)
	
	training_matrix = np.array(training_matrix, dtype = "|S15")
	test_matrix = np.array(test_matrix, dtype = "|S15")
	
	imputed_values = {}
	#data_matrix is a numpy array
	for column in range(0,training_matrix.shape[1]):
		#compute the mean for every attributes (MV can be only in test set)
		no_missing = [x for x in training_matrix[:,column] if x != '?']
		if isCategorical[column]:
			lst = Counter(no_missing)
			highest_count = max([lst[i] for i in lst])
			# if there are more values with the highest frequncy then pick a random one
			values = [i for i in lst if lst[i] == highest_count]
			random.shuffle(values)
			imputed = values[0]
			imputed_values[column] = imputed
		else:
			mean = np.mean(np.array(no_missing,dtype='f'))
			imputed_values[column] = mean
	
	#replace missing values 
	for column in range(0,training_matrix.shape[1]):
		if "?" in training_matrix[:,column]:
			imputed = str(imputed_values[column])
			training_matrix[training_matrix[:,column] == "?",column] = imputed
		if "?" in test_matrix[:,column]:
			imputed = str(imputed_values[column])
			test_matrix[test_matrix[:,column] == "?",column] = imputed
	
	#write the new arff files
	new_training_file = training_set.split(".")[0] #remove .arff
	with open (new_training_file + "_NoMissing", "wt") as output:
		output.write("@relation " + relation_name + "_NoMissingValues\n")
		output.write(header_training) 
		for index, row in enumerate(training_matrix):
			output.write(",".join(row))
			output.write("," + training_labels[index]+"\n")
	
	new_test_file = test_set.split(".")[0] #remove .arff
	with open (new_test_file + "_NoMissing", "wt") as output:
		output.write("@relation " + relation_name + "_NoMissingValues\n")
		output.write(header_test) 
		for index, row in enumerate(test_matrix):
			output.write(",".join(row))
			output.write("," + test_labels[index]+"\n")
	
	return
