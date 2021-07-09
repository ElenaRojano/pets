#! /usr/bin/env bash

current=`pwd`
hpo_file=$current/../external_data/hp.obo

#Launch Cohort Analyzer
mkdir -p $current/cohort_analyzer_results
cd $current/cohort_analyzer_results
	coPatReporter.rb -i $current/hummu_congenital_full_dataset.txt -o $current/cohort_analyzer_results -p phenotypes -c chr -d patient_id -s start -e stop -m lin
cd..
# Launch Reg2Phen
mkdir -p $current/reg2phen_results
cd $current/reg2phen_results
	reg2phen.rb -t $current/associations_file.txt -p $current/genes.txt -b $hpo_file -P -g -H -o $current/results/patient1Genes.txt -F $current/reg2phen_results/patient1Genes.html
cd ..
# Launch Phen2Reg
mkdir -p $current/phen2reg_results
cd $current/phen2reg_results
	phen2reg.rb -t $current/associations_file.txt -M 50 -p $current/example_patient_hpos.txt -k -y 0 -d prednum -i $current/hpo2ci.txt -r fisher -f $hpo_file -P 0.05 -b 1 -m -T -Q > $current/phen2reg_results/single_phens.txt
cd ..
