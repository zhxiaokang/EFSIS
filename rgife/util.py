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
import sys
import subprocess
import os
import numpy as np
import ConfigParser
import tarfile

def set_seed(seed):
	np.random.seed(seed)
	return

def get_selected_attributes(path_results):
	atts = []
	for line in open ( path_results + "BestIteration/selected_best_data.arff"):
		if line.startswith("@attribute"):
			atts.append(line.split()[1])
	return atts[:-1]
	
def write_summary(path_results, metric, best_performance):
	selected_attributes = get_selected_attributes(path_results)
	with open (path_results + "summary.txt", "wt") as output:
		output.write("RGIFE selected {0} attributes that provide {1} = {2}\n".format(len(selected_attributes), metric, best_performance))
		output.write("== Selected attributes ==\n")
		for att in selected_attributes:
			output.write(att + "\n")

def copy_final_results(path_results, best_iteration, reference_iteration, pre_processing, name, binarised_attribute_mapping, dataset_name):
	
	args = ["mkdir","-p", path_results + "BestIteration"]
	subprocess.call(args)
	args = ["mkdir","-p", path_results + "ReferenceIteration"]
	subprocess.call(args)
	
	if pre_processing == 0:
		args = ["cp", path_results + "it"+str(best_iteration)+"/"+name+"_it"+str(best_iteration)+".arff", path_results + "BestIteration/selected_best_data.arff"]
		subprocess.call(args)
		args = ["cp", path_results + "it"+str(reference_iteration)+"/"+name+"_it"+str(reference_iteration)+".arff", path_results + "ReferenceIteration/selected_reference_data..arff"]
		subprocess.call(args)
		return
		
	#if there are categorical attributes
	best_attributes = []
	reference_attributes = []
	for line in open (path_results + "it"+str(best_iteration)+"/"+name+"_it"+str(best_iteration)+".arff"):
		if line in ['\n', '\r\n', ' ']:
			continue
		if line.lower().startswith("@att"):
			attribute = line.split()[1]
			if attribute in binarised_attribute_mapping:
				attribute_to_save = binarised_attribute_mapping[attribute]
			else:
				attribute_to_save = attribute
			if attribute_to_save not in best_attributes:
				best_attributes.append(attribute_to_save)

	for line in open (path_results + "it"+str(reference_iteration)+"/"+name+"_it"+str(reference_iteration)+".arff"):
		if line in ['\n', '\r\n', ' ']:
			continue
		if line.lower().startswith("@att"):
			attribute = line.split()[1]
			if attribute in binarised_attribute_mapping:
				attribute_to_save = binarised_attribute_mapping[attribute]
			else:
				attribute_to_save = attribute
			if attribute_to_save not in reference_attributes:
				reference_attributes.append(attribute_to_save)
	data = False
	matrix = []
	indexes = []
	index = 0
	with open ( path_results + "BestIteration/selected_best_data.arff","wt") as output:
		output.write("@relation BestIt\n")
		for line in open (dataset_name):
			if line in ['\n', '\r\n', ' ']:
				continue
			if line.lower().startswith("@att"):
				att_name = line.split()[1]
				if att_name not in best_attributes:
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
	
	data = False
	matrix = []
	indexes = []
	index = 0
	with open ( path_results + "ReferenceIteration/selected_reference_data.arff","wt") as output:
		output.write("@relation ReferenceIt\n")
		for line in open (dataset_name):
			if line in ['\n', '\r\n', ' ']:
				continue
			if line.lower().startswith("@att"):
				att_name = line.split()[1]
				if att_name not in reference_attributes:
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
	
	
def get_ordinal_attributes(ordinal_attributes_list):
	ordinal_attributes = {}
	
	for line in open(ordinal_attributes_list):
		att_name = line.strip("\t \n")
		ordinal_attributes[att_name] = True
	
	return ordinal_attributes
	
def parse_configuration (configuration_file):
	block_type = "RBS"
	validation_schema = "10CV"
	trees = 3000
	missing_values = "no"
	categorical_attributes = "no"
	cs_rf = "no"
	cost = []
	metric = "robust_accuracy"
	repetitions = 1
	depth = None
	different_folds = "no"
	seed = -1
	tolerance_samples = 1
	white_list = []
	ordinal_attributes = {}
	cv_schema = "db_scv"
	config = ConfigParser.ConfigParser()
	config.read(configuration_file)
	
	if config.has_option("parameters", 'block_type'):
		block_type = config.get("parameters", "block_type").upper()
		good_values = ["RBS", "ABS"]
		if block_type not in good_values:
			 print "Invalid block_type: {0} ({1})".format(block_type, good_values)
			 sys.exit(1)
	
	if config.has_option("parameters", 'validation'):
		validation_schema = config.get("parameters", "validation").upper()
		good_values = ["10CV", "LOOCV"]
		if validation_schema not in good_values:
			 print "Invalid validation: {0} ({1})".format(validation_schema, good_values)
			 sys.exit(1)
	
	if config.has_option("parameters", 'trees'):
		try:
			trees = int(config.get("parameters", "trees"))
		except ValueError:
			 print "Invalid trees parameter:, requires int values"
			 sys.exit(1)
	 
	if config.has_option("parameters", 'missing_values'):
		missing_values = config.get ("parameters", "missing_values").lower()
		good_values = ["yes", "no"]
		if missing_values not in good_values:
			 print "Invalid missing_values: {0} ({1})".format(missing_values, good_values)
			 sys.exit(1)
					
	if config.has_option("parameters", 'categorical_attributes'):
		categorical_attributes = config.get ("parameters", "categorical_attributes")
		good_values = ["yes", "no"]
		if categorical_attributes not in good_values:
			 print "Invalid categorical_attributes: {0} ({1})".format(categorical_attributes, good_values)
			 sys.exit(1)
	
	if config.has_option("parameters", 'cs_rf'):
		cs_rf = config.get ("parameters", "cs_rf").lower()
		good_values = ["yes", "no"]
		if cs_rf not in good_values:
			 print "Invalid cs_rf: {0} ({1})".format(cs_rf, good_values)
			 sys.exit(1)
		else:
			 if cs_rf == "yes":
				cost = config.get("parameters", "misclassification_cost").split(",")
				for i,value in enumerate(cost):
					try:
						value = int(value)
					except:
						print "Invalid cost parameter:, requires int values"
						sys.exit(1)
	
	if config.has_option("parameters", 'metric'):
		metric = config.get ("parameters", "metric").lower()
		good_values = ["accuracy", "robust_accuracy", "fscore", "gmean", "auc", "overall_auc"]
		if metric not in good_values:
			 print "Invalid metric: {0} ({1})".format(metric, good_values)
			 sys.exit(1)
	
	if config.has_option("parameters", 'repetitions'):
		try:
			repetitions = int(config.get("parameters", "repetitions"))
		except ValueError:
			print "Invalid repetitions parameter:, requires int values"
			sys.exit(1)

	if config.has_option("parameters", 'max_depth'):
		try:
			depth = int(config.get("parameters", "max_depth"))
		except ValueError:
			print "Invalid max_depth parameter:, requires int values"
			sys.exit(1)
	
	if config.has_option("parameters", 'different_folds'):
		different_folds = config.get ("parameters", "different_folds").lower()
		good_values = ["yes", "no"]
		if different_folds not in good_values:
			 print "Invalid different_folds: {0} ({1})".format(different_folds, good_values)
			 sys.exit(1)
	
	if config.has_option("parameters", 'seed'):
		try:
			seed = int(config.get("parameters", "seed"))
		except ValueError:
			 print "Invalid seed parameter:, requires int values"
			 sys.exit(1)
	
	
	if config.has_option("parameters", 'tolerance_samples'):
		tolerance_samples = config.get ("parameters", "tolerance_samples")
		try:
			tolerance_samples = float(config.get("parameters", "tolerance_samples"))
		except ValueError:
			 print "Invalid tolerance_samples parameter: requires int/float values"
			 sys.exit(1)
		tolerance_samples = config.get ("parameters", "tolerance_samples")
		if "." in tolerance_samples:
			tolerance_samples = float(config.get ("parameters", "tolerance_samples"))
		else:
			tolerance_samples = int(config.get ("parameters", "tolerance_samples"))
	
	if config.has_option("parameters", 'cv_schema'):
		cv_schema = config.get ("parameters", "cv_schema").upper()
		good_values = ["DB_SCV", "SCV"]
		if cv_schema not in good_values:
			 print "Invalid cv_schema: {0} ({1})".format(cv_schema, good_values)
			 sys.exit(1)
		
	if config.has_option("parameters", 'white_list'):
		white_list = config.get ("parameters", "white_list").split(",")		

	if config.has_option("parameters", 'ordinal_attributes'):
		ordinal_attributes_list = config.get ("parameters", "ordinal_attributes")
		ordinal_attributes = get_ordinal_attributes(ordinal_attributes_list)
	
	
	return block_type, validation_schema, trees, depth, missing_values, categorical_attributes, white_list, cs_rf, cost, metric, repetitions, ordinal_attributes, different_folds, tolerance_samples, seed, cv_schema

def make_tarfile():
    
    with tarfile.open(os.getcwd()+"/iterations.tar.gz", "w:gz") as tar:
        files = os.listdir(".")
        for folder in files:
			if folder.startswith("it"):
				tar.add(folder)
