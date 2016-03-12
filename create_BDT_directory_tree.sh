#!/bin/bash

declare -a NTrees=("NTrees=10" "NTrees=20" "NTrees=40" "NTrees=60" "NTrees=80" "NTrees=100")
declare -a BGs=("bg_qcd" "bg_zll" "bg_wjets_ev" "bg_wjets_muv" "bg_wjets_tauv" "bg_top" "bg_vv" "bg_zjets_vv")
declare -a SeparationType=("SeparationType=GiniIndex" "SeparationType=CrossEntropy" 
							"SeparationType=MisClassificationError" "SeparationType=SDivSqrtSPlusB")
declare -a AdaBoostBeta=("AdaBoostBeta=0.1" "AdaBoostBeta=0.2" "AdaBoostBeta=0.4" 
						 "AdaBoostBeta=0.5" "AdaBoostBeta=0.6" "AdaBoostBeta=0.8")
declare -a nCuts=("nCuts=5" "nCuts=7" "nCuts=10" "nCuts=12" "nCuts=15" "nCuts=20")
declare -a params=("NTrees" "SeparationType" "AdaBoostBeta" "nCuts")

# Create directory "analysis" if it doesn't exist
if [ ! -d "$analysis" ]; then
	echo "Creating directory analysis..."
    mkdir analysis
fi

# Loop through params, create directories in analysis
for i in "${params[@]}"
do
   if [ ! -d "analysis/BDT_varying_$i" ]; then
   	 echo "Creating directory analysis/BDT_varying_$i..."
     mkdir analysis/BDT_varying_$i
   else
   	 echo "Directory $i already exists."
   fi
   
   for j in "${BGs[@]}"
   do
   	  if [ ! -d "analysis/BDT_varying_$i/$j" ]; then
   	 	echo "Creating directory analysis/BDT_varying_$i/$j..."
     	   mkdir analysis/BDT_varying_$i/$j
   	  else
   	 	echo "Directory analysis/BDT_varying_$i/$j already exists."
   	  fi

   done   
done