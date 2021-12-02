#! /usr/bin/env ruby
# Rojano E. & Seoane P., March 2019
# Code to prepare data to get the associations between pathological phenotypes (HPO) and genomic regions (SOR)

ROOT_PATH = File.dirname(__FILE__)
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'pets'))

##############################
#LIBRARIES
##############################
require 'optparse'
require 'constants.rb'
require 'pets'

###############################
#METHODS
###############################

def build_tripartite_network(patient_data, patients_by_cluster, add_parents)
	tripartite_network = []
	patients_by_cluster.each do |patient, node_ids|
		node_ids.each do |node_id|
			tripartite_network << [node_id, patient] 
		end
	end
	hpo = Cohort.get_ontology(Cohort.act_ont)
	patient_data.each_profile do |id, profile|
		profile = profile.map{|term| hpo.get_ancestors(term)}.flatten.uniq if add_parents == 'root'
	  profile.each do |term|
			tripartite_network << [term, id] 
	  end
	end
	return tripartite_network
end

##############################
#OPTPARSE
##############################

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"

  options[:chromosome_col] = nil
  opts.on("-c", "--chromosome_col INTEGER/STRING", "Column name if header is true, otherwise 0-based position of the column with the chromosome") do |data|
    options[:chromosome_col] = data
  end

  options[:id_col] = nil
  opts.on("-d", "--pat_id_col INTEGER/STRING", "Column name if header is true, otherwise 0-based position of the column with the patient id") do |data|
    options[:id_col] = data
  end

  options[:end_col] = nil
  opts.on("-e", "--end_col INTEGER/STRING", "Column name if header is true, otherwise 0-based position of the column with the end mutation coordinate") do |data|
    options[:end_col] = data
  end
  
 	options[:ont_col] = nil
  opts.on("-p", "--hpo_term_col INTEGER/STRING", "Column name if header true or 0-based position of the column with the HPO terms") do |data|
    options[:ont_col] = data
  end

  options[:start_col] = nil
  opts.on("-s", "--start_col INTEGER/STRING", "Column name if header is true, otherwise 0-based position of the column with the start mutation coordinate") do |data|
    options[:start_col] = data
  end

  options[:separator] = '|'
  opts.on("-S", "--hpo_separator STRING", "Set which character must be used to split the HPO profile. Default '|'") do |data|
    options[:separator] = data
  end

  options[:names] = false
  opts.on("-n", "--hpo_names", "Define if the input HPO are human readable names. Default false") do
    options[:names] = true
  end

  options[:header] = true
  opts.on("-H", "--header", "File has a line header. Default true") do 
    options[:header] = false
  end

  #===================================================================

  options[:input_file] = nil
  opts.on("-i", "--input_file PATH", "Input file with patients for parsing phenotypes to HPO codes") do |value|
    options[:input_file] = value
  end

  options[:output_file] = 'tripartite_network.txt'
  opts.on("-o", "--output_file PATH", "Output file for the tripartite network") do |value|
    options[:output_file] = value
  end

  options[:cluster_file] = 'cluster_coords.txt'
  opts.on("-u", "--cluster_file PATH", "Cluster coords output file that will be used to translate SOR nodes") do |value|
    options[:cluster_file] = File.basename(value)
  end 

  options[:excluded_hpo] = nil
  opts.on("-x", "--excluded_hpo PATH", "List of HPO phenotypes to exclude (low informative)") do |excluded_hpo|
    options[:excluded_hpo] = excluded_hpo
  end

  options[:tag] = 'A'
  opts.on("-m", "--mutation_type STRING", "Type of patient mutation, either it is a deletion (d) or duplication (D)") do |type|
    options[:tag] = type
  end

  options[:hpo_file] = nil
  opts.on("-O", "--hpo_file PATH", "Input HPO file for extracting HPO codes") do |value|
    options[:hpo_file] = value
  end

  options[:add_parents] = nil
  opts.on("-r", "--parents STRING", "'root' to add all parents until the ontology root. 'coh' to add parents until the most general term in the cohort(TODO).") do |value|
    options[:add_parents] = value
  end 

  opts.on_tail("-h", "--help", "Show this message") do
    puts opts
    exit
  end

end.parse!

###############################
#MAIN
###############################
output_folder = File.dirname(File.expand_path(options[:output_file]))
Dir.mkdir(output_folder) if !File.exists?(output_folder)

hpo_file = options[:hpo_file]
hpo_file = !ENV['hpo_file'].nil? ? ENV['hpo_file'] : HPO_FILE if hpo_file.nil?
Cohort.load_ontology(:hpo, hpo_file, options[:excluded_hpo])
Cohort.act_ont = :hpo

patient_data, rejected_hpos_L, rejected_patients_L = Cohort_Parser.load(options)
rejected_hpos_C, rejected_patients_C = patient_data.check
rejected_hpos = rejected_hpos_L | rejected_hpos_C
rejected_patients = rejected_patients_L + rejected_patients_C
patient_data.remove_incomplete_records
patient_data.index_vars
patients_by_cluster, sors = patient_data.generate_cluster_regions(:reg_overlap, options[:tag], 1)
tripartite_network = build_tripartite_network(patient_data, patients_by_cluster, options[:add_parents])

write_array(rejected_hpos, File.join(output_folder, 'missing_hpo_names'))
write_array(sors, File.join(output_folder, options[:cluster_file]))
write_array(tripartite_network, options[:output_file])