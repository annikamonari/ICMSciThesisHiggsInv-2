#!/bin/bash

declare -a BGs=("bg_qcd" "bg_zll" "bg_wjets_ev" "bg_wjets_muv" "bg_wjets_tauv" "bg_top" "bg_vv" "bg_zjets_vv")
declare -a params=("HiddenLayers" "NCycles" "LearningRate")

# Create directory "analysis" if it doesn't exist
if [ ! -d "$analysis" ]; then
	echo "Creating directory analysis..."
    mkdir analysis
fi

# Loop through params, create directories in analysis
for i in "${params[@]}"
do
   if [ ! -d "analysis/MLP_varying_$i" ]; then
   	 echo "Creating directory analysis/MLP_varying_$i..."
     mkdir analysis/MLP_varying_$i
   else
   	 echo "Directory $i already exists."
   fi
   
   for j in "${BGs[@]}"
   do
   	  if [ ! -d "analysis/MLP_varying_$i/$j" ]; then
   	 	echo "Creating directory analysis/MLP_varying_$i/$j..."
     	   mkdir analysis/MLP_varying_$i/$j
   	  else
   	 	echo "Directory analysis/MLP_varying_$i/$j already exists."
   	  fi

   done   
done