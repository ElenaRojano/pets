#! /usr/bin/env ruby

# Script to compare a patient cohort against MONDO ontology defined diseases
# @author Fernando Moreno Jabato <jabato(at)uma(dot)es>

ROOT_PATH = File.dirname(__FILE__)
EXTERNAL_DATA = File.expand_path(File.join(ROOT_PATH, '..', 'external_data'))
MONDO_FILE = File.join(EXTERNAL_DATA, 'mondo.obo')
HPO_FILE = File.join(EXTERNAL_DATA, 'hp.obo')
EXTERNAL_CODE = File.expand_path(File.join(ROOT_PATH, '..', 'external_code'))
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'pets'))

require 'optparse'
require 'semtools'
require 'coPatReporterMethods.rb'

##########################
# FUNCTIONS
##########################
def read_patient_profiles(file, header: true, col_sep: "\t", profile_sep: "|")
  patients = {}
  File.open(file).each do |line|
    if header
      header = false
    else
    	line.chomp!
    	id, hpos = line.split(col_sep)
    	hpos = hpos.split(profile_sep)
    	patients[id.to_sym] = hpos
    end
  end
  return patients
end


##########################
# OPT-PARSER
##########################

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"

  options[:input_file] = nil
  opts.on("-i", "--input_file PATH", "Input file with patient data") do |data|
    options[:input_file] = data
  end

  options[:output_file] = nil
  opts.on("-o", "--output_file PATH", "Output file with comparisson data") do |data|
    options[:output_file] = data
  end

  options[:verbose] = false
  opts.on("-v", "--verbose", "Activate verbose mode") do
    options[:verbose] = true
  end


  opts.on_tail("-h", "--help", "Show this message") do
    puts opts
    exit
  end

end.parse!


##########################
# MAIN
##########################
# Load ontologies
puts("Loading MONDO ontology ...") if options[:verbose] ### Verbose point
mondo = Ontology.new(file: MONDO_FILE, load_file: true)
puts("MONDO ontology loaded\nLoading Human Phenotype Ontology ...") if options[:verbose] ### Verbose point
hpo = Ontology.new(file: HPO_FILE, load_file: true)
puts("HP ontology loaded") if options[:verbose] ### Verbose point
# Load cohort 
puts("Loading patients :: " + options[:input_file])
patients = read_patient_profiles(options[:input_file])
puts("Patients loaded (" + patients.length.to_s + ")") if options[:verbose] ### Verbose point
patients.each do |id, hpos|
	patients[id], rejected_codes = hpo.check_ids(hpos.map{|term| term.to_sym}, substitute: false)
  warn("Patient (" + id.to_s + ") contains not allowed terms.")if !rejected_codes.empty?
end
hpo.load_profiles(patients)
hpo.clean_profiles(store: true, remove_alternatives: false)

# Take MONDO profiles
mondo.calc_dictionary(:xref, select_regex: /(HP:[0-9]*)/, store_tag: :HP, multiterm: true, substitute_alternatives: false)
mondo_profiles = mondo.dicts[:HP][:byTerm].clone
mondo_profiles = mondo_profiles.each{|mondo,hpos| mondo_profiles[mondo] = hpos.map{|hp| hp.to_sym}}
puts("Obtaining MONDO profiles with HPOs ("+ mondo_profiles.length.to_s + ")") if options[:verbose] ### Verbose point

# Compare
puts("Comparing patients against MONDO profiles. It can take a while ...") if options[:verbose] ### Verbose point
sims = hpo.compare_profiles(external_profiles: mondo_profiles, sim_type: :resnick, ic_type: :resnick, bidirectional: false, against_external: true) # Compare patients agains MONDO

# Export
puts("Exporting results") if options[:verbose] ### Verbose point
sim_pairs_file = File.join(options[:output_file],"sim_pairs")
sims_pairs = write_profile_pairs(sims, sim_pairs_file)
# File.open(options[:output_file], "w") { |f| f.write sims_matrix.to_json }

# Plot
puts("Rendering plots ...") if options[:verbose] ### Verbose point
system("#{File.join(EXTERNAL_CODE, 'plot_heatmap.R')} -d #{sim_pairs_file} -o #{File.join(options[:output_file],"sim_heatmap")} -m max -p -s -c 'MONDO terms' -r 'Patients'")    



puts("Program finish") if options[:verbose] ### Verbose point
