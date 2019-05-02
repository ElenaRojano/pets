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

def loadPatientFile(patient_file, hpo_storage, hpo_dictionary, add_parents)
	patient2phenotype = {}
	hpo_count = {}
	not_found = []	
	patients_genomic_region_by_chr = {}
	File.open(patient_file).each do |line|
		line.chomp!
		next if line.include?("#")
		patient, chr, start, stop, phenotype_profile = line.split("\t", 5)
		next if phenotype_profile.nil? #For skipping patients without phenotypes
		phenotypes = phenotype_profile.split('|')
		phenotypes.each do |hpo_name|
			hpo_code = hpo_dictionary[hpo_name]
			if hpo_code.nil?
				not_found << hpo_name if !not_found.include?(hpo_name)
			else
				get_all_hpos(patient, hpo_code, patient2phenotype, hpo_storage, hpo_count, add_parents)
			end
		end
		info = [patient, start.to_i, stop.to_i]
		add_record(patients_genomic_region_by_chr, chr, info)

	end
	return patient2phenotype, hpo_count, not_found, patients_genomic_region_by_chr
end

def add_record(hash, key, record)
	query = hash[key]
	if query.nil?
		hash[key] = [record]
	elsif !query.include?(record)
		query << record
	end
end

# 08/04/19
def get_all_hpos(patient, hpo_code, patient2phenotype, hpo_storage, hpo_count, add_parents)
	add_record(hpo_count, hpo_code, patient)
	add_record(patient2phenotype, patient, hpo_code)
	if add_parents
		hpo_parent_codes = hpo_storage[hpo_code][2]
    	hpo_parent_codes.each do |parent_code|
			get_all_hpos(patient, parent_code, patient2phenotype, hpo_storage, hpo_count, add_parents)
    	end	
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

def generate_cluster_regions(patients_genomic_region_by_chr, mutation_type)
	patients_by_cluster = {}
	sors = []
	patients_genomic_region_by_chr.each do |chrm, genomic_ranges|
		reference = get_reference(genomic_ranges) # Get putative overlap regions
		overlapping_patients = overlap_patients(genomic_ranges, reference) # See what patient has match with a overlap region
		clust_number = 1
		reference.each_with_index do |ref, i|
			current_patients = overlapping_patients[i]
			if current_patients.length > 1
				ref << chrm
				node_identifier = "#{chrm}.#{clust_number}.#{mutation_type}.#{current_patients.length}"
				ref << node_identifier
				save_sor(current_patients, node_identifier, patients_by_cluster)
				sors << ref
				clust_number += 1
			end
		end
	end
	return patients_by_cluster, sors
end

def save_sor(current_patients, node_identifier, patients_by_cluster)
	current_patients.each do |patient|
		add_record(patients_by_cluster, patient, node_identifier)
	end
end

def build_tripartite_network(patients2hpo, hpo_stats, ic_threshold, patients_by_cluster)
	tripartite_network = []
	patients_by_cluster.each do |patient, node_ids|
		node_ids.each do |node_id|
			tripartite_network << [node_id, patient] 
		end
	end
	patients_list = patients_by_cluster.keys
	patients2hpo.each do |patient, code|
		if patients_list.include?(patient)
		  code.each do |c|
			tripartite_network << [c, patient] if hpo_stats[c].last >= ic_threshold 
		  end
		end
	end
	return tripartite_network
end

def compute_hpo_stats(hpo_count, threshold_ic, patient_number)
	hpo_stats = {}
	patient_hpo_ic = []
	hpo_count.each do |hpo_code, patient_ids|
	    hpo_freq = patient_ids.length.fdiv(patient_number) #hpo frequency in patients
	    hpo_ic = -Math.log10(hpo_freq)
	    if hpo_ic >= threshold_ic
	    		hpo_stats[hpo_code] = [hpo_freq, hpo_ic]
	    end
	    patient_ids.each do |patient_id|
	    	patient_hpo_ic << [patient_id, hpo_code, hpo_ic]
	    end
	end
	return hpo_stats, patient_hpo_ic.sort{|a,b| a.first.to_i <=> b.first.to_i}
end

def write_hash(hash, file_path, header = [])
	File.open(file_path, 'w') do |handler|
 		handler.puts header.join("\t") if !header.empty?
 		hash.each do |key, array|
 			handler.puts "#{key}\t#{array.join("\t")}"
 		end
 	end
end

def write_array(array, file_path)
	File.open(file_path, 'w') do |handler|
		array.each do |record|
			if record.class == String
				line = record
			else
				line = record.join("\t")
			end
			handler.puts line  
		end
	end
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

  options[:patient_file] = nil
  opts.on("-i", "--input_file PATH", "Input file with patients for parsing phenotypes to HPO codes") do |value|
    options[:patient_file] = value
  end

  options[:mutation_type] = 'A'
  opts.on("-m", "--mutation_type STRING", "Type of patient mutation, either it is a deletion (d) or duplication (D)") do |type|
    options[:mutation_type] = type
  end

  options[:output_file] = 'tripartite_network.txt'
  opts.on("-o", "--output_file", "Output file for the tripartite network") do |value|
    options[:output_file] = value
  end 

  options[:hpo_stat_file] = 'hpo_stats.txt'
  opts.on("-s", "--hpo_stat_file", "Output file with HPO codes, their frequency and CI") do |value|
    options[:hpo_stat_file] = value
  end 

  options[:hpo_file] = nil
  opts.on("-p", "--hpo_file PATH", "Input HPO file for extracting HPO codes") do |value|
    options[:hpo_file] = value
  end

  options[:add_parents] = true
  opts.on("-r", "--no_parents", "Switch for not including HPO parents in results") do
    options[:add_parents] = false
  end 

  options[:thresold] = 1
  opts.on("-t", "--info_thresold FLOAT", "IC thresold to discard non informative hpo") do |thresold|
    options[:thresold] = thresold.to_f
  end

end.parse!

###############################
#MAIN
###############################
hpo_black_list = load_hpo_black_list(options[:excluded_hpo])
hpo_storage = load_hpo_file(options[:hpo_file], hpo_black_list)
hpo_dictionary = create_hpo_dictionary(hpo_storage)
patients2hpo, hpo_count, not_found, chr_patients_genomic_region = loadPatientFile(options[:patient_file], hpo_storage, hpo_dictionary, options[:add_parents])
hpo_stats, patient_hpo_ic = compute_hpo_stats(hpo_count, options[:thresold], patients2hpo.length)
patients_by_cluster, sors = generate_cluster_regions(chr_patients_genomic_region, options[:mutation_type])
tripartite_network = build_tripartite_network(patients2hpo, hpo_stats, options[:thresold], patients_by_cluster)

write_array(not_found - hpo_black_list, 'missing_hpo_names')
write_array(sors, options[:cluster_file])
write_hash(hpo_stats, options[:hpo_stat_file], %w[HPOcode Frequency IC])
write_array(tripartite_network, options[:output_file])
write_array(patient_hpo_ic, 'filtered_hpo.txt')

