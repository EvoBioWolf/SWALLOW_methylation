#!/bin/bash

awk '$6 != "NA"' /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/1_Population_Structure/4_pixy/HR.168.invariant/HR.168.invariant.eu_dxy.filter.txt > /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/1_Population_Structure/4_pixy/2_DXY/HR.168.invariant.eu_dxy.noNA.txt