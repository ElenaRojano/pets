#! /usr/bin/env ruby

ROOT_PATH = File.dirname(__FILE__)
EXTERNAL_DATA = File.expand_path(File.join(ROOT_PATH, '..', 'external_data'))
HPO_FILE = File.join(EXTERNAL_DATA, 'hp.obo')
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'gephepred'))

require 'generalMethods.rb'
require 'optparse'

###############
#METHODS
###############

def translate_codes_to_terms(patient_data, hpo_storage)
	patients_with_hpo_names = {}
	hpo_names = []
	patient_data.each do |patientID, hpos_and_cnvs|
		hpos = hpos_and_cnvs.shift.split('|')
		hpos.each do |hpo|
			hpo_names << hpo_storage[hpo][1]
		end
		hpos_and_cnvs << hpo_names.join('|')
		patients_with_hpo_names[patientID] = hpos_and_cnvs
		hpo_names = []
	end
	return patients_with_hpo_names
end

def save_translated_file(patients_with_hpo_names, output_file)
	handler = File.open(output_file, 'w')
	patients_with_hpo_names.each do |id, data|
		patientID = id.gsub(/_i[0-9]/,'')
		handler.puts "#{patientID}\t#{data.join("\t")}"
	end
	handler.close
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
  opts.on("-H", "--header", "Set if the file has a line header. Default true") do 
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


end.parse!


###############
#MAIN
###############

hpo_storage = load_hpo_file(HPO_FILE)
patient_data, $patient_number = load_patient_cohort(options)
patients_with_hpo_names = translate_codes_to_terms(patient_data, hpo_storage)

save_translated_file(patients_with_hpo_names, options[:output_file])


Process.exit	
