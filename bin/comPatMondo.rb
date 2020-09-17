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
require 'csv'
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


def read_expand_file(file)
  expand_table = {}
  header = true
  CSV.open(file, "r", { :col_sep => "\t" }).each do |row|
      # Row Format
      #  0 - Disease DB
      #  1 - Disease identifier
      #  2 - Disease name
      #  3 - Negation (qualifier)
      #  4 - HPO ID
      #  5 - Reference
      #  6 - Evidence code
      #  7 - Onset
      #  8 - Frequency HPO
      #  9 - Modifier
      #  10 - Sub-ontology
      #  11 - Alternative names
      #  12 - Curators
      #  13- Frequency Raw
      #  14 - Sex
      # row = row[0].split("\t")
    if header
      header = false
    else
      dis_id = row[5]
      dis_id = row[0] + ':' + row[1] if dis_id.nil?
      dis_id.gsub!('ORPHA','Orphanet')
      dis_id.gsub!('OMIM','OMIMPS')
      if !expand_table.include?(dis_id.to_sym)
        expand_table[dis_id] = [row[4].chomp.to_sym]
      else
        expand_table[dis_id] << row[4].chomp.to_sym
      end
    end
  end
  return expand_table
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
  opts.on("-o", "--output_file PATH", "Output path") do |data|
    options[:output_file] = data
  end
  
  options[:expand] = nil
  opts.on("-m", "--expand FILE", "Expand file with HPO annotation for OMIM/Orphanet diseases") do |data|
    options[:expand] = data
  end

  options[:verbose] = false
  opts.on("-v", "--verbose", "Activate verbose mode") do
    options[:verbose] = true
  end

  options[:export] = false
  opts.on("-e", "--export", "Flag to export ontology items as JSON files") do 
    options[:export] = true
  end

  options[:timer] = false
  opts.on("-t", "--time", "Calculates time used for main steps") do 
    options[:timer] = true
  end

  opts.on_tail("-h", "--help", "Show this message") do
    puts opts
    exit
  end

end.parse!


##########################
# MAIN
##########################
if options[:verbose]
  puts "Launch configuration:"
  puts "\tInput file :: " + options[:input_file].inspect
  puts "\tOutput path :: " + options[:output_file].inspect
  puts "\tExpand file :: " + options[:expand].inspect
  puts "\tVerbose mode :: " + options[:verbose].inspect
  puts "\tExport flag :: " + options[:export].inspect
  puts "\tTimer flag :: " + options[:timer].inspect
end

timeflag_init = Process.clock_gettime(Process::CLOCK_MONOTONIC) if options[:timer] # TIME FLAG

# Load ontologies
puts("Loading MONDO ontology ...") if options[:verbose] ### Verbose point
mondo = Ontology.new(file: MONDO_FILE, load_file: true)
puts("MONDO ontology loaded\nLoading Human Phenotype Ontology ...") if options[:verbose] ### Verbose point
hpo = Ontology.new(file: HPO_FILE, load_file: true)
puts("HP ontology loaded") if options[:verbose] ### Verbose point
# Load cohort 
timeflag_ontologyLoaded = Process.clock_gettime(Process::CLOCK_MONOTONIC) if options[:timer] # TIME FLAG
puts("Loading patients :: " + options[:input_file])
patients = read_patient_profiles(options[:input_file])
puts("Patients loaded (" + patients.length.to_s + ")") if options[:verbose] ### Verbose point
patients.each do |id, hpos|
	patients[id], rejected_codes = hpo.check_ids(hpos.map{|term| term.to_sym}, substitute: false)
  warn("Patient (" + id.to_s + ") contains not allowed terms.")if !rejected_codes.empty?
end
hpo.load_profiles(patients, reset_stored: false)
hpo.clean_profiles(store: true, remove_alternatives: false)
timeflag_patientsLoaded = Process.clock_gettime(Process::CLOCK_MONOTONIC) if options[:timer] # TIME FLAG

# Take MONDO profiles
mondo.calc_dictionary(:xref, select_regex: /(HP:[0-9]*)/, store_tag: :HP, multiterm: true, substitute_alternatives: false)
mondo_profiles = mondo.dicts[:HP][:byTerm].clone
mondo_profiles = mondo_profiles.each{|mondo,hpos| mondo_profiles[mondo] = hpos.map{|hp| hp.to_sym}}
puts("Obtaining MONDO profiles with HPOs ("+ mondo_profiles.length.to_s + ")") if options[:verbose] ### Verbose point
timeflag_mondoProfiles = Process.clock_gettime(Process::CLOCK_MONOTONIC) if options[:timer] # TIME FLAG

# Expand MONDO profiles if proceeds
if(!options[:expand].nil?)
  puts("Expanding MONDO profiles by disease relations ...") if options[:verbose] ### Verbose point
  # Import file
  expand_table = read_expand_file(options[:expand])
  # Prepare disease relations
  mondo.calc_dictionary(:xref, select_regex: /(OMIMPS:[0-9]*)/, store_tag: :OMIM, multiterm: true, substitute_alternatives: false)
  mondo.calc_dictionary(:xref, select_regex: /(Orphanet:[0-9]*)/, store_tag: :Orpha, multiterm: true, substitute_alternatives: false)
  # Find
  expandable_diseases = [mondo.dicts[:OMIM][:byValue].keys, mondo.dicts[:Orpha][:byValue].keys].flatten
  puts("\tDiseases linked to MONDO (" + expandable_diseases.length.to_s + ")") if options[:verbose] ### Verbose point
  expandable_diseases = expandable_diseases & expand_table.keys
  # Expand
  puts("\tDiseases linked to MONDO and HPOs (" + expandable_diseases.length.to_s + ")") if options[:verbose] ### Verbose point
  mondo_profiles_ids = mondo_profiles.keys
  expanded_profiles = [] if options[:verbose]
  expandable_diseases.each do |disease_ID|
    query = mondo.dicts[:OMIM][:byValue][disease_ID] # OMIM
    if query.nil? # Orphanet
      query = mondo.dicts[:Orpha][:byValue][disease_ID]
    end
    expandable_profiles = query
    expanded_profiles << expandable_profiles if options[:verbose]
    expandable_profiles.each do |prof|
      if mondo_profiles_ids.include?(prof) 
        mondo_profiles[prof] = [mondo_profiles[prof],expand_table[disease_ID]].flatten.uniq
      else
        mondo_profiles[prof] = expand_table[disease_ID]
        mondo_profiles_ids << prof        
      end
    end
  end
  timeflag_mondoExpanded = Process.clock_gettime(Process::CLOCK_MONOTONIC) if options[:timer] # TIME FLAG
  if options[:verbose]
    expanded_profiles = expanded_profiles.flatten.uniq
    new_profiles = expanded_profiles - mondo.dicts[:HP][:byTerm].keys
    puts "\tNew profiles (" + new_profiles.length.to_s + "). Expansion process affected profiles (" + expanded_profiles.length.to_s + "/" + mondo_profiles.length.to_s + ")"
  end
end

# Clean MONDO profiles
mondo_profiles_ids = mondo_profiles.keys
affected_ids = []
mondo_profiles_ids.each do |id|
  prof, rejected = hpo.check_ids(mondo_profiles[id], substitute: false)
  if prof.empty?
    mondo_profiles.delete(id)
  elsif !rejected.empty?
    mondo_profiles[id] = prof
    affected_ids << id
  end
end
timeflag_mondoCleaned = Process.clock_gettime(Process::CLOCK_MONOTONIC) if options[:timer] # TIME FLAG

puts "Cleaning MONDO profiles process has left (" + mondo_profiles.length.to_s + "/" + mondo_profiles_ids.length.to_s + ") and has modified (" + affected_ids.length.to_s + ") profiles" if options[:verbose]


# Compare
puts("Comparing patients against MONDO profiles. It can take a while ...") if options[:verbose] ### Verbose point
sims = hpo.compare_profiles(external_profiles: mondo_profiles, sim_type: :resnick, ic_type: :resnick, bidirectional: false, against_external: true) # Compare patients agains MONDO
timeflag_sims = Process.clock_gettime(Process::CLOCK_MONOTONIC) if options[:timer] # TIME FLAG

# Export
puts("Exporting results") if options[:verbose] ### Verbose point
sim_pairs_file = File.join(options[:output_file],"sim_pairs")
sims_pairs = write_profile_pairs(sims, sim_pairs_file)
if(options[:export])
  puts("Exporting ontologies to JSON") if options[:verbose] ### Verbose point
  mondo.write(File.join(options[:output_file],"mondo.json"))
  hpo.write(File.join(options[:output_file],"hpo.json")) 
end
timeflag_export = Process.clock_gettime(Process::CLOCK_MONOTONIC) if options[:timer] # TIME FLAG

# Plot
puts("Rendering plots ...") if options[:verbose] ### Verbose point
system("#{File.join(EXTERNAL_CODE, 'plot_heatmap.R')} -d #{sim_pairs_file} -o #{File.join(options[:output_file],"sim")} -m max -p -s -c 'MONDO terms' -r 'Patients'")    
timeflag_render = Process.clock_gettime(Process::CLOCK_MONOTONIC) if options[:timer] # TIME FLAG

if options[:timer]
  times = {}
  times[:loadOntologies] = timeflag_ontologyLoaded - timeflag_init
  times[:patientsLoadAndClean] = timeflag_patientsLoaded - timeflag_ontologyLoaded 
  times[:mondoProfilesCalc] = timeflag_mondoProfiles - timeflag_patientsLoaded
  if !options[:expand].nil?
    times[:mondoExpansion] = timeflag_mondoExpanded - timeflag_mondoProfiles
    times[:mondoCleaning] = timeflag_mondoCleaned - timeflag_mondoExpanded
  else
    times[:mondoCleaning] = timeflag_mondoCleaned - timeflag_mondoProfiles
  end
  times[:simCalc] = timeflag_sims - timeflag_mondoCleaned
  times[:export] = timeflag_export - timeflag_sims
  times[:plots] = timeflag_render - timeflag_export
  # Show
  puts " <<< TIME ELAPSED >>>"
  times.each{|k,v| puts("Time elapsed for (" + k.to_s + ") : " + v.round(2).to_s + " s")}
end


puts("Program finish") if options[:verbose] ### Verbose point
