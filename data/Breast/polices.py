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

import sys

def get_attributes(path):
	attributes = []
	for line in open (path):
		if line.startswith("RGIFE"):
			selected_attributes = int(line.split()[2])
			best_performance = float(line.split()[-1])
			continue
		if line.startswith("=="):
			continue
		attributes.append(line.split()[0])
	return best_performance, selected_attributes, attributes

def get_signatures():
	union = set()

selected_attributes = []
performances = []

path_results = str(sys.argv[1])
runs = int(sys.argv[2])
max_run = 0
min_run = 0
min_best = 0
max_best = 0
min_atts = 10000000000
max_atts = 0
selected_attributes_runs = []
union = set()

for run in range(1,runs + 1):
	path = "{0}/run{1}/summary.txt".format(path_results, run)
	performance, selected_attributes, attributes = get_attributes(path)
	
	for attribute in attributes:
		union.add(attribute)
	selected_attributes_runs.append(attributes)
	
	if selected_attributes < min_atts or (selected_attributes == min_atts and performance > min_best ):
		min_atts = selected_attributes
		min_run = run
		min_best = performance
	
	if selected_attributes > max_atts or (selected_attributes == max_atts and performance > max_best ):
		max_atts = selected_attributes
		max_run = run
		max_best = performance

with open (path_results + "max_model.txt", "wt") as output:
	for attribute in selected_attributes_runs[max_run - 1]:
		output.write("{0}\n".format(attribute))
output.close()
with open (path_results + "min_model.txt", "wt") as output:
	for attribute in selected_attributes_runs[min_run - 1]:
		output.write("{0}\n".format(attribute))
output.close()
with open (path_results + "union_model.txt", "wt") as output:
	for attribute in union:
		output.write("{0}\n".format(attribute))
output.close()
