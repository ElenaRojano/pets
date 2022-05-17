#! /usr/bin/env ruby

ROOT_PATH = File.dirname(__FILE__)
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'pets'))

require 'optparse'
require 'pets'

#############################
## METHODS
#############################
def load_index(path_index)
  vcf_index = {}
  File.open(path_index).each do |line|
    id, path = line.chomp.split("\t")
    vcf_index[id] = path
  end
  return vcf_index
end


##########################
#OPT-PARSER
##########################

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

  options[:genome_assembly] = 'hg38'
  opts.on("-G", "--genome_assembly STRING", "Genome assembly version. Please choose between hg18, hg19 and hg38. Default hg38") do |data|
    options[:genome_assembly] = data
  end

  options[:header] = true
  #chr\tstart\tstop
  opts.on("-H", "--header", "Set if the file has a line header. Default true") do 
    options[:header] = false
  end

  options[:input_file] = nil
  opts.on("-i", "--input_file PATH", "Input file with patient data") do |data|
    options[:input_file] = data
  end

  options[:vcf_index] = nil
  opts.on("-I", "--vcf_index PATH", "VCF file with patient id pointing to vcf path") do |data|
    options[:vcf_index] = data
  end

  options[:names] = false
  opts.on("-n", "--hpo_names", "Define if the input HPO are human readable names. Default false") do
    options[:names] = true
  end

  options[:output_folder] = nil
  opts.on("-o", "--output_file PATH", "Output folder") do |data|
    options[:output_folder] = data
  end

  options[:ont_col] = nil
  opts.on("-p", "--hpo_term_col INTEGER/STRING", "Column name if header true or 0-based position of the column with the HPO terms") do |data|
  	options[:ont_col] = data
  end

  options[:separator] = '|'
  opts.on("-S", "--hpo_separator STRING", "Set which character must be used to split the HPO profile. Default '|'") do |data|
  	options[:separator] = data
  end

  options[:start_col] = nil
  opts.on("-s", "--start_col INTEGER/STRING", "Column name if header is true, otherwise 0-based position of the column with the start mutation coordinate") do |data|
  	options[:start_col] = data
  end

  opts.on_tail("-h", "--help", "Show this message") do
    puts opts
    exit
  end

end.parse!

#############################################################
## MAIN
#############################################################
hpo_file = !ENV['hpo_file'].nil? ? ENV['hpo_file'] : HPO_FILE
Cohort.load_ontology(:hpo, hpo_file, options[:excluded_hpo])
Cohort.act_ont = :hpo

patient_data, rejected_hpos_L, rejected_patients_L = Cohort_Parser.load(options)
rejected_hpos_C, rejected_patients_C = patient_data.check(hard=true)
patient_data.link2ont(Cohort.act_ont)

vcf_index = load_index(options[:vcf_index]) if !options[:vcf_index].nil?
patient_data.export_phenopackets(options[:output_folder], options[:genome_assembly], vcf_index: vcf_index)
