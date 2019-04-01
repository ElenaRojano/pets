#! /usr/bin/env ruby
# Rojano E. & Seoane P., March 2019
# Code to prepare data to get the associations between pathological phenotypes (HPO) and genomic regions (SOR)

ROOT_PATH = File.dirname(__FILE__)
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'gephepred'))

##############################
#LIBRARIES
##############################
require 'generalMethods.rb'
require 'optparse'
require "benchmark"

###############################
#METHODS
###############################

# def load_hpo_file(hpo_file)
# 	hpo_storage = []
# 	id = nil
# 	name = nil
# 	alt_id = []
# 	syn = []
# 	is_a = []
# 	File.open(hpo_file).each do |line|
# 		line.chomp!
# 		tag, info = line.split(': ')
# 		if tag == 'id' || tag == 'name' || tag == 'is_a' || tag == 'synonym' || tag == 'alt_id'
# 			if tag == 'id'
# 				hpo_storage << [id, alt_id.join('|'), name, syn.join('|')].concat(is_a) if !name.nil?  #if !temp[1].include?("obsolete") 
# 				id = info
# 				name = nil
# 				alt_id = []
# 				syn = []
# 				is_a = []
# 			end
# 			if tag == 'alt_id'
# 				alt_id << info
# 			elsif tag == 'is_a'
# 				is_a.concat(info.split(' ! '))
# 			elsif tag == 'synonym'
# 				syn << info.split('"')[1]
# 			else
# 				name = info
# 			end
# 		end
# 	end
# 	hpo_storage << [id, alt_id.join('|'), name, syn.join('|')].concat(is_a)
# 	return hpo_storage
# end

def create_hpo_dictionary(hpo_storage, hpo_black_list)
	storage = {}
	hpo_storage.each do |hpos|
		hpo_code = hpos.shift
		next if hpo_black_list.include?(hpo_code)
		alt_hpo_code = hpos.shift
		phenotype = hpos.shift
		synonyms = hpos.shift 
		relations = []
		hpos.each_slice(2) do |pair|
			#pair = HPO code, phenotype
			relations << pair
		end
		storage[phenotype] = [hpo_code, relations]
		if !synonyms.nil?
			synonyms.split('|').each do |syn|
				storage[syn] = [hpo_code, relations]
			end
		end
	end
	return storage
end

def load_hpo_black_list(excluded_hpo_file)
	excluded_hpos = []
	File.open(excluded_hpo_file).each do |line|
		line.chomp!
		excluded_hpos << line
	end
	return excluded_hpos
end

def loadPatientFile(patient_file, hpo_ontology, parents)
	patients = {}
	hpo_stats = {}
	not_found = []	
	chr_patients_genomic_region = {}
	File.open(patient_file).each do |line|
		line.chomp!
		next if line.include?("#")
		patient, chr, start, stop, phenotype_profile = line.split("\t", 5)
		next if phenotype_profile.nil? #For skipping patients without phenotypes
		phenotypes = phenotype_profile.split('|')
		phenotypes.each do |phenotype|
			get_all_hpos(patient, phenotype, patients, hpo_ontology, hpo_stats, not_found, parents)
		end
		query = chr_patients_genomic_region[chr]
		info = [patient, start.to_i, stop.to_i]
		if query.nil?
			chr_patients_genomic_region[chr] = [info]
		else
			query << info
		end

	end
	return patients, hpo_stats, not_found, chr_patients_genomic_region
end

def get_all_hpos(patient, phenotype, patients, hpo_ontology, hpo_stats, not_found, parents)
	query = hpo_ontology[phenotype]
	if !query.nil?
		hpo_code, relations = query
		query_stats = hpo_stats[hpo_code] # Do tracking of patients that have an hpo
		if query_stats.nil?
			hpo_stats[hpo_code] = [patient]
		elsif !query_stats.include?(patient)
			query_stats << patient
		end
		query_patient = patients[patient]
		if query_patient.nil?
		        patients[patient] = [hpo_code]
		else
		        query_patient << hpo_code
		end
		if !relations.nil? && parents # ADDING PARENTAL PHENOTYPES TO PATIENT
		    relations.each do |rel_code, rel_name|
		        get_all_hpos(patient, rel_name, patients, hpo_ontology, hpo_stats, not_found, parents)
	        end
	    end
	else
		not_found << phenotype if !not_found.include?(phenotype)
	end
end

def get_reference(genomic_ranges)
	#genomic_ranges = [patientID, mut_start, mut_stop]
	reference = []
	reference.concat(genomic_ranges.map{|gr| gr[1]})# get start
	reference.concat(genomic_ranges.map{|gr| gr[2]})# get stop
	reference.uniq!
	reference.sort!
	#Define overlap range
	final_reference = []
	reference.each_with_index do |coord,i|
		next_coord = reference[i + 1]
		final_reference << [coord, next_coord] if !next_coord.nil? 
	end
	return final_reference
end

def overlap_patients(genomic_ranges, reference)
	overlaps = []
	reference.each do |start, stop|
		patients = []
		genomic_ranges.each do |pt_id, pt_start, pt_stop|
			if (start <= pt_start && stop >= pt_stop) ||
				(start > pt_start && stop < pt_stop) ||
				(stop > pt_start && stop <= pt_stop) ||
				(start >= pt_start && start < pt_stop)
				patients << pt_id
			end
		end
		overlaps << patients.uniq
	end
	return overlaps
end

##############################
#OPTPARSE
##############################

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"

  options[:cluster_file] = 'cluster_coords.txt'
  opts.on("-c", "--cluster_file", "Cluster coords output file that will be used to translate SOR nodes") do |value|
    options[:cluster_file] = value
  end 

  options[:excluded_hpo] = nil
  opts.on("-e", "--excluded_hpo PATH", "List of HPO phenotypes to exclude (low informative)") do |excluded_hpo|
    options[:excluded_hpo] = excluded_hpo
  end

  options[:do_freq] = false
  opts.on("-f", "--do_freq", "Switch for calculate HPO frequency instead of IC") do
    options[:do_freq] = true
  end 

  options[:patient_file] = nil
  opts.on("-i", "--input_file PATH", "Input file with patients for parsing phenotypes to HPO codes") do |value|
    options[:patient_file] = value
  end

  options[:mutation_type] = 'A'
  opts.on("-m", "--mutation_type STRING", "Type of patient mutation, either it is a deletion (d) or duplication (D)") do |type|
    options[:mutation_type] = type
  end

  options[:nodes_file] = 'nodes.txt'
  opts.on("-n", "--nodes_file", "SOR-patient nodes output file for the tripartite network") do |value|
    options[:nodes_file] = value
  end 

  options[:output_file] = 'patient2hpo'
  opts.on("-o", "--output_file", "Output file with patients and their HPO codes") do |value|
    options[:output_file] = value
  end 

  options[:hpo_file] = nil
  opts.on("-p", "--hpo_file PATH", "Input HPO file for extracting HPO codes") do |value|
    options[:hpo_file] = value
  end

  options[:parents] = true
  opts.on("-r", "--no_parents", "Switch for not including HPO parents in results") do
    options[:parents] = false
  end 

  options[:thresold] = 1
  opts.on("-t", "--info_thresold FLOAT", "Thresold to discard non informative hpo") do |thresold|
    options[:thresold] = thresold.to_f
  end

end.parse!

###############################
#MAIN
###############################
hpo_black_list = load_hpo_black_list(options[:excluded_hpo])
hpo_storage = load_hpo_file(options[:hpo_file], false) #from generalMethods.rb, false to return only hpo_storage
hpoNameDictionary = create_hpo_dictionary(hpo_storage, hpo_black_list)
patients2hpo, hpo_stats, not_found, chr_patients_genomic_region = loadPatientFile(options[:patient_file], hpoNameDictionary, options[:parents])

patients_by_cluster = []
handler = File.open(options[:cluster_file], 'w')
chr_patients_genomic_region.each do |chrm, genomic_ranges|
	reference = get_reference(genomic_ranges) # Get putative overlap regions
	overlapping_patients = overlap_patients(genomic_ranges, reference) # See what patient has match with a overlap region
	clust_number = 0
	reference.each_with_index do |ref, i|
		current_patients = overlapping_patients[i]
		if current_patients.length > 1
			ref << chrm
			node_identifier = "#{chrm}.#{clust_number + 1}.#{options[:mutation_type]}.#{current_patients.length}"
			patients_by_cluster << [current_patients, node_identifier]
			ref << node_identifier
			handler.puts ref.join("\t")
			clust_number += 1
		end
	end
end
handler.close

handler = File.open(options[:nodes_file], 'w')
patients_by_cluster.each do |patients, node_id|
	patients.each do |patient|
		handler.puts "#{node_id}\t#{patient}" 
	end
end
not_found = not_found - hpo_black_list
File.open('missing_hpo_names', 'w'){|f| f.puts not_found}
number_patients = patients2hpo.length
handler.close

handler = File.open(options[:output_file], 'w')
patients2hpo.each do |patient, code|
	stat = nil
	result = nil
  code.uniq.each do |c|
    stat = hpo_stats[c].length.to_f / number_patients #hpo frequency in patients
    result = -Math.log10(stat)
    if result >= options[:thresold]
    	if options[:do_freq]
    		handler.puts "#{c}\t#{stat}"
    	else
    		handler.puts "#{patient}\t#{c}\t#{result}"  
    	end
    end
  end
end
handler.close
