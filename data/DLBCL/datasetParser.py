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

def read_arff_file (name):
	class_mapping = {}
	attributes = []
	samples = []
	labels = []
	previous_line = ""
	data = False
	ids = []
	for line in open (name):
		if line in ['\n', '\r\n', ' ']:
			continue
		line = line.strip("\n\r")
		if line.lower().startswith("@att"):
			attribute = line.split()[1]
			attributes.append(attribute)
			previous_line = line
			continue
		if line.lower().startswith("@data"):
			classes = previous_line.split("{")[1].strip("}").split(",")
			for i,label in enumerate(classes):
				class_mapping[label.strip(" ").strip("\n").strip("}")] = i
			data = True
			continue
		if data:
			label = line.strip("\n\r").split(",")[-1]
			labels.append(class_mapping[label])
			sample = [float(i) for i in line.split(",")[:-1]]
			samples.append(sample)
	return samples, labels, attributes, class_mapping


def get_attributes_info(name, attribute_definitions, attributes_type, isCategorical, attribute_indexes):
	index = 0
	previous_line = ""
	for line in open (name):
		#skip empty lines
		if line in ['\n', '\r\n', ' ']:
			continue
		if line.startswith("@"):
			field = line.split()[0].lower().strip("@")
			if field == "relation":
				relation_name = line.split()[1]
			if field  == "attribute":
				previous_line = line
				attribute_name = line.split()[1]
				attribute_type = line.split()[2]
				attribute_definitions.append(line)
				if attribute_type.lower() == "numeric" or attribute_type.lower() == "real"  or attribute_type.lower() == "string" :
					attributes_type[attribute_name] = attribute_type.lower()
					isCategorical[index] = False
				elif attribute_type.startswith("{"):
					cat_values = line.split("{")[1].strip("\n\r }")
					values = cat_values.replace(" ","").split(",")
					attributes_type[attribute_name] = values
					isCategorical[index] = True
				attribute_indexes[attribute_name] = index
				index += 1
				continue
			if field == "data":
				class_name = previous_line.split()[1]
				continue
	return relation_name, class_name


def get_attributes_name(dataset_name):
	attributes_name = []
	for line in open (dataset_name):
		if line.lower().startswith ("@att"):
			attributes_name.append (line.split()[1])
		continue
	return attributes_name[:-1] #remove class attribute	

def get_num_samples(dataset_name):
	count = 0
	for line in open (dataset_name):
		line = line.strip(' \n\r')
		if line.startswith("@") or line.startswith("%") or len(line)==0:
			continue
		else:
			count += 1
	return count

def get_num_classes(dataset_name):
	previous_line = ""
	for line in open (dataset_name):
		if line in ['\n', '\r\n', ' ']:
			continue
		if line.lower().startswith("@att"):
			previous_line = line
			continue
		if line.lower().startswith("@data"):
			classes = previous_line.split("{")[1].strip("}").split(",")
			return len(classes)

def get_nattributes (reference_iteration, path_results, name):
	attributes_name = []
	for line in open (path_results + "it"+str(reference_iteration)+"/"+name+"_it"+str(reference_iteration)+".arff"):
		if line.lower().startswith ("@att"):
			attributes_name.append (line.split()[1])
		continue
	return len(attributes_name[:-1]) 

def get_index_atts(path_total):
	indexes = {}
	index = 0
	for line in open (path_total):
		if line.lower().startswith("@att"):
			line = line.strip("\n\r")
			att_name = line.split()[1]
			indexes[att_name] = index
			index += 1
	return indexes
