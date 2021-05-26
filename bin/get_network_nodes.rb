#! /usr/bin/env ruby
# Rojano E. & Seoane P., March 2019
# Code to prepare data to get the associations between pathological phenotypes (HPO) and genomic regions (SOR)

ROOT_PATH = File.dirname(__FILE__)
EXTERNAL_DATA = File.expand_path(File.join(ROOT_PATH, '..', 'external_data'))
HPO_FILE = File.join(EXTERNAL_DATA, 'hp.obo')
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'pets'))

##############################
#LIBRARIES
##############################
require 'generalMethods.rb'
require 'optparse'
require 'semtools'

###############################
#METHODS
###############################

def loadPatientFile(patient_file, hpo, add_parents)
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
		# phenotypes, rejected = hpo.translate_names2codes(phenotypes)
		phenotypes, rejected = hpo.translate_names(phenotypes)
		not_found = not_found | rejected
		phenotypes.each do |hpo_code|
			get_all_hpos(patient, hpo_code, patient2phenotype, hpo, hpo_count, add_parents) if !hpo.is_removable(hpo_code)
		end
		info = [patient, start.to_i, stop.to_i]
		add_record(patients_genomic_region_by_chr, chr, info)
	end
	return patient2phenotype, hpo_count, not_found, patients_genomic_region_by_chr
end


def get_all_hpos(patient, hpo_code, patient2phenotype, hpo, hpo_count, add_parents)
	add_record(hpo_count, hpo_code, patient)
	add_record(patient2phenotype, patient, hpo_code)
	if add_parents == 'root'
		# hpo_parent_codes = hpo.get_parents(hpo_code)
		hpo_parent_codes = hpo.get_ancestors(hpo_code)
    	hpo_parent_codes.each do |parent_code|
			add_record(hpo_count, parent_code, patient)
			add_record(patient2phenotype, patient, parent_code)
    	end	
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

def compute_hpo_stats(hpo_count, patient_number)
	hpo_stats = {}
	patient_hpo_ic = []
	hpo_count.each do |hpo_code, patient_ids|
	    hpo_freq = patient_ids.length.fdiv(patient_number) #hpo frequency in patients
	    hpo_ic = -Math.log10(hpo_freq)
	    hpo_stats[hpo_code] = [hpo_freq, hpo_ic]
	    patient_ids.each do |patient_id|
	    	patient_hpo_ic << [patient_id, hpo_code, hpo_ic]
	    end
	end
	return hpo_stats, patient_hpo_ic.sort{|a,b| a.first <=> b.first}
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
  opts.on("-c", "--cluster_file PATH", "Cluster coords output file that will be used to translate SOR nodes") do |value|
    options[:cluster_file] = File.basename(value)
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
  opts.on("-o", "--output_file PATH", "Output file for the tripartite network") do |value|
    options[:output_file] = value
  end 

  options[:hpo_file] = nil
  opts.on("-p", "--hpo_file PATH", "Input HPO file for extracting HPO codes") do |value|
    options[:hpo_file] = value
  end

  options[:add_parents] = nil
  opts.on("-r", "--parents STRING", "'root' to add all parents until the ontology root. 'coh' to add parents until the most general term in the cohort(TODO).") do |value|
    options[:add_parents] = value
  end 

  options[:hpo_stat_file] = 'hpo_stats.txt'
  opts.on("-s", "--hpo_stat_file PATH", "Output file with HPO codes, their frequency and CI") do |value|
    options[:hpo_stat_file] = File.basename(value)
  end 

  options[:thresold] = 0
  opts.on("-t", "--info_thresold FLOAT", "IC thresold to discard non informative hpo. Default: 0.") do |thresold|
    options[:thresold] = thresold.to_f
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
hpo_file = ENV['hpo_file'] if hpo_file.nil?
hpo_file = HPO_FILE if hpo_file.nil?

# hpo = Ontology.new
# hpo.load_black_list(options[:excluded_hpo]) if !options[:excluded_hpo].nil?
# hpo.load_data(hpo_file)
if !options[:excluded_hpo].nil?
	hpo = Ontology.new(file: hpo_file, load_file: true, removable_terms: read_excluded_hpo_file(options[:excluded_hpo]))
else
	hpo = Ontology.new(file: hpo_file, load_file: true)
end
patients2hpo, hpo_count, not_found, chr_patients_genomic_region = loadPatientFile(options[:patient_file], hpo, options[:add_parents])

hpo_stats, patient_hpo_ic = compute_hpo_stats(hpo_count, patients2hpo.length)
patients_by_cluster, sors = generate_cluster_regions(chr_patients_genomic_region, options[:mutation_type])

tripartite_network = build_tripartite_network(patients2hpo, hpo_stats, options[:thresold], patients_by_cluster)

# write_array(not_found - hpo.excluded_codes, File.join(output_folder, 'missing_hpo_names'))
write_array(not_found - hpo.removable_terms, File.join(output_folder, 'missing_hpo_names'))
write_array(sors, File.join(output_folder, options[:cluster_file]))
write_hash(hpo_stats.select{|hp_code, stats| stats.last > options[:thresold]}, File.join(output_folder, options[:hpo_stat_file]), %w[HPOcode Frequency IC])
write_array(tripartite_network, options[:output_file])
write_array(patient_hpo_ic, File.join(output_folder, 'filtered_hpo.txt'))