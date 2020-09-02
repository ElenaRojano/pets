#! /usr/bin/env ruby

ROOT_PATH = File.dirname(__FILE__)
EXTERNAL_DATA = File.expand_path(File.join(ROOT_PATH, '..', 'external_data'))
HPO_FILE = File.join(EXTERNAL_DATA, 'hp.obo')
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'pets'))

require 'generalMethods.rb'
require 'optparse'
require 'semtools'

###############
#METHODS
###############

def translate_hpo(patient_data, hpo, translate)
	patients_with_hpo_names = {}
	patient_data.each do |patientID, patient_record|
		hpos, chr, start, stop = patient_record
    if translate == 'names'
      # hpos, rejected = hpo.translate_codes2names(hpos)
      hpos, rejected = hpo.translate_ids(hpos)
    elsif translate =='codes'
      # hpos, rejected = hpo.translate_names2codes(hpos)
      hpos, rejected = hpo.translate_names(hpos)
      STDERR.puts(" The ontology names '#{rejected.join(',')}' were not found") if !rejected.empty?
    end
    patient_record[0] = hpos
  end
end

def save_translated_file(patients_with_hpo_names, output_file)
	File.open(output_file, 'w') do |f|
  	patients_with_hpo_names.each do |id, patient_record|
      hpos, chr, start, stop = patient_record
  		id = id.gsub(/_i[0-9]+$/,'')
      f.puts "#{id}\t#{hpos.join('|')}\t#{[chr, start, stop].join("\t")}"
  	end
	end
end

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

  options[:pat_id_col] = nil
  opts.on("-d", "--pat_id_col INTEGER/STRING", "Column name if header is true, otherwise 0-based position of the column with the patient id") do |data|
    options[:pat_id_col] = data
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

  options[:hpo_col] = nil
  opts.on("-p", "--hpo_term_col INTEGER/STRING", "Column name if header true or 0-based position of the column with the HPO terms") do |data|
  	options[:hpo_col] = data
  end

  options[:start_col] = nil
  opts.on("-s", "--start_col INTEGER/STRING", "Column name if header is true, otherwise 0-based position of the column with the start mutation coordinate") do |data|
  	options[:start_col] = data
  end

  options[:hpo_separator] = '|'
  opts.on("-S", "--hpo_separator STRING", "Set which character must be used to split the HPO profile. Default '|'") do |data|
    options[:hpo_separator] = data
  end

  options[:translate] = nil
  opts.on("-t", "--translate STRING", "Set 'names' to translate from hpo codes to names or set 'codes' to translate from hpo names to codes. By default, ther is not translation") do |data|
    options[:translate] = data
  end
end.parse!

###############
#MAIN
###############
hpo_file = ENV['hpo_file']
hpo_file = HPO_FILE if hpo_file.nil?

patient_data = load_patient_cohort(options)
if !options[:translate].nil?
  # hpo = Ontology.new
  # hpo.load_data(hpo_file)
  hpo = OBO_Handler.new(file: hpo_file, load_file: true)
  translate_hpo(patient_data, hpo, options[:translate])
end
save_translated_file(patient_data, options[:output_file])
