#! /usr/bin/env ruby

# Script to compare a patient cohort against MONDO ontology defined diseases
# @author Fernando Moreno Jabato <jabato(at)uma(dot)es>

ROOT_PATH = File.dirname(__FILE__)
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'pets'))

require 'optparse'
require 'semtools'
require 'csv'
require 'constants.rb'
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


def read_infer_file(file, col_sep: "\t", comment_char: "#")
  infer_table = {}
  CSV.open(file, "r", { :col_sep => col_sep }).each do |row|
    if row[0][0] != comment_char
      dis_id = row[0]
      dis_id.gsub!('ORPHA','Orphanet')
      dis_id.gsub!('OMIM','OMIMPS')
      if !infer_table.include?(dis_id.to_sym)
        infer_table[dis_id] = [row[3].chomp.to_sym]
      else
        infer_table[dis_id] << row[3].chomp.to_sym
      end
    end
  end
  return infer_table
end


def read_add_file(file, col_sep: "\t", comment_char: "#")
  add_table = {}
  CSV.open(file, "r", { :col_sep => col_sep }).each do |row|
    if row[0][0] != comment_char
      mondo_id = row[0]
      mondo_id = mondo_id.to_sym
      if !add_table.include?(mondo_id)
        add_table[mondo_id] = [row[1].chomp.to_sym]
      else
        add_table[mondo_id] << row[1].chomp.to_sym
      end
    end
  end
  return add_table
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
  
  options[:add_file] = nil
  opts.on("-a", "--add_file PATH", "File with direct MONDO - HPO relations to be added") do |data|
    options[:add_file] = data
  end

  options[:infer] = nil
  opts.on("-N", "--infer FILE", "Infer file with HPO annotation for OMIM/Orphanet diseases") do |data|
    options[:infer] = data
  end

  options[:verbose] = false
  opts.on("-v", "--verbose", "Activate verbose mode") do
    options[:verbose] = true
  end

  options[:export] = false
  opts.on("-e", "--export", "Flag to export ontology items as JSON files") do 
    options[:export] = true
  end

  options[:expand] = false
  opts.on("-E", "--expand", "Flag to expand by parentals using HPO") do 
    options[:expand] = true
  end

  options[:timer] = false
  opts.on("-t", "--time", "Calculates time used for main steps") do 
    options[:timer] = true
  end

  options[:plots] = false
  opts.on("-p", "--plots", "Render plots") do 
    options[:plots] = true
  end

  allowed_sims = {"resnik" => :resnik, "lin" => :lin, "jiang_conrath" => :jiang_conrath}
  options[:similitude] = allowed_sims["resnik"]
  opts.on("-S", "--similitude TYPE", "Specify similitude. Allowed: resnik, resnik_observed, seco, zhou, sanchez. Default: resnik") do |data|
    sim = allowed_sims[data]
    if sim.nil?
      warn 'Specified similitude type is not allowed. Resnik will be used'
    else
      options[:similitude] = sim
    end      
  end

  allowed_ics = {"resnik" => :resnik, "resnik_observed" => :resnik_observed, 
                  "seco" => :seco, "zhou" => :zhou, "sanchez" => :sanchez}
  options[:ic] = allowed_sims["resnik"]
  opts.on("-I", "--ic TYPE", "Specify ic. Allowed: resnik, lin, jiang_conrath. Default: resnik") do |data|
    ic = allowed_ics[data]
    if ic.nil?
      warn 'Specified IC type is not allowed. Resnik will be used'
    else
      options[:ic] = ic
    end      
  end

  options[:hpo] = nil
  opts.on("-H", "--hpo FILE", "HPO JSON file") do |data|
    options[:hpo] = data
  end

  options[:mondo] = nil
  opts.on("-M", "--mondo FILE", "MONDO JSON file") do |data|
    options[:mondo] = data
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
  puts "\tInfer file :: " + options[:infer].inspect
  puts "\tAdd file :: " + options[:add_file].inspect
  puts "\tHPO JSON file :: " + options[:hpo].inspect
  puts "\tMONDO JSON file :: " + options[:mondo].inspect
  puts "\tIC type :: " + options[:ic].inspect
  puts "\tSimilitude type :: " + options[:similitude].inspect
  puts "\tVerbose mode :: " + options[:verbose].inspect
  puts "\tExport flag :: " + options[:export].inspect
  puts "\tExpand flag :: " + options[:expand].inspect
  puts "\tTimer flag :: " + options[:timer].inspect
end

timeflag_init = Process.clock_gettime(Process::CLOCK_MONOTONIC) if options[:timer] # TIME FLAG

# Load ontologies
puts("Loading MONDO ontology ...") if options[:verbose] ### Verbose point
if options[:mondo].nil?
  mondo = Ontology.new(file: MONDO_FILE, load_file: true)
else
  mondo = Ontology.new
  mondo.read(options[:mondo])
end
puts("MONDO ontology loaded\nLoading Human Phenotype Ontology ...") if options[:verbose] ### Verbose point
if options[:hpo].nil?
  hpo = Ontology.new(file: HPO_FILE, load_file: true)
else
  hpo = Ontology.new
  hpo.read(options[:hpo])
end
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
hpo.load_profiles(patients, reset_stored: true)
hpo.clean_profiles(store: true, remove_alternatives: false)
timeflag_patientsLoaded = Process.clock_gettime(Process::CLOCK_MONOTONIC) if options[:timer] # TIME FLAG

# Take MONDO profiles
mondo.calc_dictionary(:xref, select_regex: /(HP:[0-9]*)/, store_tag: :HP, multiterm: true, substitute_alternatives: false)
mondo_profiles = mondo.dicts[:HP][:byTerm].clone
mondo_profiles = mondo_profiles.each{|mondo,hpos| mondo_profiles[mondo] = hpos.map{|hp| hp.to_sym}}
puts("Obtaining MONDO profiles with HPOs ("+ mondo_profiles.length.to_s + ")") if options[:verbose] ### Verbose point
timeflag_mondoProfiles = Process.clock_gettime(Process::CLOCK_MONOTONIC) if options[:timer] # TIME FLAG

# Infer info to MONDO profiles if proceeds
if(!options[:infer].nil?)
  timeflag_init_infer = Process.clock_gettime(Process::CLOCK_MONOTONIC) if options[:timer] # TIME FLAG
  puts("Infering new info to MONDO profiles by disease relations ...") if options[:verbose] ### Verbose point
  # Import file
  infer_table = read_infer_file(options[:infer])
  # Prepare disease relations
  mondo.calc_dictionary(:xref, select_regex: /(OMIMPS:[0-9]*)/, store_tag: :OMIM, multiterm: true, substitute_alternatives: false)
  mondo.calc_dictionary(:xref, select_regex: /(Orphanet:[0-9]*)/, store_tag: :Orpha, multiterm: true, substitute_alternatives: false)
  # Find
  inferable_diseases = [mondo.dicts[:OMIM][:byValue].keys, mondo.dicts[:Orpha][:byValue].keys].flatten
  puts("\tDiseases linked to MONDO (" + inferable_diseases.length.to_s + ")") if options[:verbose] ### Verbose point
  inferable_diseases = inferable_diseases & infer_table.keys
  # Expand
  puts("\tDiseases linked to MONDO and HPOs (" + inferable_diseases.length.to_s + ")") if options[:verbose] ### Verbose point
  mondo_profiles_ids = mondo_profiles.keys
  infered_profiles = [] if options[:verbose]
  inferable_diseases.each do |disease_ID|
    query = mondo.dicts[:OMIM][:byValue][disease_ID] # OMIM
    if query.nil? # Orphanet
      query = mondo.dicts[:Orpha][:byValue][disease_ID]
    end
    inferable_profiles = query
    infered_profiles << inferable_profiles if options[:verbose]
    inferable_profiles.each do |prof|
      if mondo_profiles_ids.include?(prof) 
        mondo_profiles[prof] = [mondo_profiles[prof],infer_table[disease_ID]].flatten.uniq
      else
        mondo_profiles[prof] = infer_table[disease_ID]
        mondo_profiles_ids << prof        
      end
    end
  end
  timeflag_mondoExpanded = Process.clock_gettime(Process::CLOCK_MONOTONIC) if options[:timer] # TIME FLAG
  if options[:verbose]
    infered_profiles = infered_profiles.flatten.uniq
    new_profiles = infered_profiles - mondo.dicts[:HP][:byTerm].keys
    puts "\tNew profiles (" + new_profiles.length.to_s + "). Inferation process affected profiles (" + infered_profiles.length.to_s + "/" + mondo_profiles.length.to_s + ")"
  end
  timeflag_end_infer = Process.clock_gettime(Process::CLOCK_MONOTONIC) if options[:timer] # TIME FLAG
end

# Add direct MONDO profiles
if !options[:add_file].nil?
  timeflag_init_add = Process.clock_gettime(Process::CLOCK_MONOTONIC) if options[:timer] # TIME FLAG
  # Load
  add_table = read_add_file(options[:add_file])
  add_table.select!{|k,v| mondo.exists?(k)}
  puts("Adding MONDO relations by external file. A total of (" + add_table.length.inspect + ") profiles with external info are available") if options[:verbose] ### Verbose point
  new_items = add_table.select{|k,v| !mondo_profiles.keys.include?(k)}
  updatable_items = add_table.select{|k,v| mondo_profiles.keys.include?(k)}
  # Update
  mondo_profiles.merge!(new_items) if new_items.length > 0
  updatable_items.each{|k,v| mondo_profiles[k] = (mondo_profiles[k] + v).flatten.uniq} if updatable_items.length > 0
  puts("Adding process has added (" + new_items.length.inspect + ") new MONDOs and has updated ("+ updatable_items.length.inspect + ") profiles") if options[:verbose] ### Verbose point
  timeflag_end_add = Process.clock_gettime(Process::CLOCK_MONOTONIC) if options[:timer] # TIME FLAG
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

# Store
mondo.load_item_relations_to_terms(mondo_profiles, true)

# Expand by parentals
if options[:expand]
  timeflag_init_expand = Process.clock_gettime(Process::CLOCK_MONOTONIC) if options[:timer] # TIME FLAG
  old_length = mondo.items.length
  puts "Expanding items to parentals" if options[:verbose]
  mondo.expand_items_to_parentals(ontology: hpo)
  if mondo.items.length > old_length && options[:verbose]
    puts "New parentals without previous relations have been added (" + (mondo.items.length - old_length).inspect + ")"
  end
  timeflag_end_expand = Process.clock_gettime(Process::CLOCK_MONOTONIC) if options[:timer] # TIME FLAG
end

# Compare
puts("Comparing patients against MONDO profiles. It can take a while ...") if options[:verbose] ### Verbose point
sims = hpo.compare_profiles(external_profiles: mondo_profiles, sim_type: options[:similitude], ic_type: options[:ic], bidirectional: false, against_external: true) # Compare patients agains MONDO
timeflag_sims = Process.clock_gettime(Process::CLOCK_MONOTONIC) if options[:timer] # TIME FLAG

# Export
puts("Exporting results") if options[:verbose] ### Verbose point
sim_pairs_basename = "sim_pairs_ic"+ options[:ic].to_s + "_sim" + options[:similitude].to_s
sim_pairs_basename += "_added" if options[:add_file]
sim_pairs_basename += "_inferred" if options[:infer]
sim_pairs_basename += "_expanded" if options[:expand]
sim_pairs_file = File.join(options[:output_file],sim_pairs_basename)
sims_pairs = write_profile_pairs(sims, sim_pairs_file)
if(options[:export])
  puts("Exporting ontologies to JSON") if options[:verbose] ### Verbose point
  mondo.write(File.join(options[:output_file],"mondo.json"))
  hpo.write(File.join(options[:output_file],"hpo.json")) 
end
timeflag_export = Process.clock_gettime(Process::CLOCK_MONOTONIC) if options[:timer] # TIME FLAG

# Plot
if options[:plots]
  puts("Rendering plots ...") if options[:verbose] ### Verbose point
  system("#{File.join(EXTERNAL_CODE, 'plot_heatmap.R')} -d #{sim_pairs_file} -o #{File.join(options[:output_file],sim_pairs_basename)} -m max -p -s -c 'MONDO terms' -r 'Patients'")    
  timeflag_render = Process.clock_gettime(Process::CLOCK_MONOTONIC) if options[:timer] # TIME FLAG
end

if options[:timer]
  times = {}
  times[:loadOntologies] = timeflag_ontologyLoaded - timeflag_init
  times[:patientsLoadAndClean] = timeflag_patientsLoaded - timeflag_ontologyLoaded 
  times[:mondoProfilesCalc] = timeflag_mondoProfiles - timeflag_patientsLoaded
  times[:mondoAdd] = timeflag_end_add - timeflag_init_add if !options[:add_file].nil?
  times[:mondoInfer] = timeflag_end_infer - timeflag_init_infer if !options[:infer].nil?
  times[:mondoExpansion] = timeflag_end_expand - timeflag_init_expand if options[:expand]
  times[:mondoCleaning] = timeflag_mondoCleaned - timeflag_mondoProfiles
  times[:simCalc] = timeflag_sims - timeflag_mondoCleaned
  times[:export] = timeflag_export - timeflag_sims
  times[:plots] = timeflag_render - timeflag_export if options[:plots]
  # Show
  puts " <<< TIME ELAPSED >>>"
  times.each{|k,v| puts("Time elapsed for (" + k.to_s + ") : " + v.round(2).to_s + " s")}
end


puts("Program finish") if options[:verbose] ### Verbose point
