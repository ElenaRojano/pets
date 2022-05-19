#! /usr/bin/env ruby

ROOT_PATH = File.dirname(__FILE__)
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'pets'))

require 'optparse'
require 'pets'

###############
#OPTIONS
###############

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
  
  options[:header] = true
  opts.on("-H", "--header", "File has a line header. Default true") do 
    options[:header] = false
  end

  options[:output_file] = 'paco_file_with_hpo_names.txt'
  opts.on("-o", "--output_file PATH", "Output paco file with HPO names") do |data|
    options[:output_file] = data
  end

  options[:input_file] = nil
  opts.on("-P", "--input_file PATH", "Input file with PACO extension") do |value|
    options[:input_file] = value
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

  options[:n_phens] = nil
  opts.on("--n_phens INTEGER", "Remove records with N or less phenotypes") do |data|
    options[:n_phens] = data.to_i
  end

  options[:save_mode] = :default
  opts.on("-m", "--save_mode STRING", "Set output data mode") do |data|
    options[:save_mode] = data.to_sym
  end

  options[:names] = false
  opts.on("-n", "--hpo_names", "Define if the input HPO are human readable names. Default false") do
    options[:names] = true
  end

  options[:translate] = false
  opts.on("-t", "--translate", "Set to translate from hpo codes to names. By default, ther is not translation") do 
    options[:translate] = true
  end
end.parse!

###############
#MAIN
###############
hpo_file = !ENV['hpo_file'].nil? ? ENV['hpo_file'] : HPO_FILE
Cohort.load_ontology(:hpo, hpo_file)
Cohort.act_ont = :hpo

patient_data, rejected_hpos, rejected_patients = Cohort_Parser.load(options)
rejected_patients_by_phen = patient_data.filter_by_term_number(options[:n_phens]) if !options[:n_phens].nil?
patient_data.save(options[:output_file], options[:save_mode], options[:translate])
