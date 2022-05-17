#! /usr/bin/env ruby

ROOT_PATH = File.dirname(__FILE__)
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'pets'))

require 'fileutils'
require 'optparse'
require 'report_html'
require 'semtools'
require 'pets'

#############################################################################################
## METHODS
############################################################################################
def get_evidence_coordinates(entity, genomic_coordinates, candidates_ids)
	all_coordinates = genomic_coordinates[entity]
	coords = all_coordinates.select{|id, coordinates| candidates_ids.include?(id.to_sym)}
	return coords
end

def make_report(profile_id, all_candidates, all_genomic_coordinates, similarity_matrixs, evidences, prof_vars, template, output)
	container = {
		profile_id: profile_id,
		candidates: all_candidates.each{|c| c[0] = c.first.to_s},
		genomic_coordinates: all_genomic_coordinates.transform_values{|c| c.first(2) },
		similarity_matrixs: similarity_matrixs,
		evidences: evidences,
		var_ids: prof_vars.nil? ? nil : prof_vars.keys.map{|i| [i, 0]},
		var_coordinates: prof_vars
	}
	report = Report_html.new(container, 'Evidence profile report')
	report.build(template)
	report.write(File.join(output, profile_id.to_s + '.html'))
end

def get_genome_hotspots(similarity_matrixs, all_genomic_coordinates)
	regions = Genomic_Feature.new(all_genomic_coordinates.values.map{|g| g[0..2]})
	candidates_by_window, genome_windows = regions.generate_cluster_regions(:reg_overlap, 'A', 1)
	# TODO: COMPLETE UNTIL FULL PREDICTOR
end


#############################################################################################
## OPTPARSE
############################################################################################

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"

  options[:profiles_file] = nil
  opts.on("-p", "--profiles_file PATH", "Path to profiles file. One profile per line and HP terms must be comma separated") do |item|
    options[:profiles_file] = item
  end

  options[:evidence_file] = nil
  opts.on("-e", "--evidence_file PATH", "Path to evidence file. The file must have a column with the evidence id and a second column with the profile with HP terms comma separated") do |item|
    options[:evidence_file] = item
  end

  options[:coordinates_file] = nil
  opts.on("-g", "--genomic_coordinates_file PATH", "Path to file with genomic coordinates for each genomic element in evidence file. One genomic element per line with format id, chr and start position.") do |item|
    options[:coordinates_file] = item
  end

  options[:evidences] = nil
  opts.on("-E", "--evidence_folder PATH", "Path to evidence folder.") do |item|
    options[:evidences] = item
  end

  options[:output_folder] = 'evidence_reports'
  opts.on("-o", "--output_folder PATH", 'Folder to save reports from profiles') do |item|
    options[:output_folder] = item
  end

  options[:variant_data] = nil
  opts.on("-V", "--variant_data PATH", 'Folder to tables of patient variants') do |item|
    options[:variant_data] = item
  end

 opts.on_tail("-h", "--help", "Show this message") do
    puts opts
    exit
  end

end.parse!

#############################################################################################
## MAIN
############################################################################################

hpo_file = !ENV['hpo_file'].nil? ? ENV['hpo_file'] : HPO_FILE
hpo = Ontology.new
hpo.read(hpo_file)

profiles = load_profiles(options[:profiles_file], hpo)
profile_variants = options[:variant_data].nil? ? {} : load_variants(options[:variant_data])
evidences, genomic_coordinates = load_evidences(options[:evidences], hpo)

hpo.load_profiles(profiles)
evidences_similarity = {}
evidences.each do |pair, data|
	entity, profile_type = pair.split('_')
	if profile_type == 'HP'
		evidence_profiles = data[:prof]
		evidence_profiles.transform_keys!{|prof_id, terms| prof_id.to_sym}
		evidences_similarity[pair] = hpo.compare_profiles(external_profiles: evidence_profiles, sim_type: :lin, bidirectional: false)
	end
end 

template = File.open(File.join(REPORT_FOLDER, 'evidence_profile.erb')).read
FileUtils.mkdir_p(options[:output_folder])
profiles.each do |profile_id, reference_prof|
	all_candidates = []
	all_genomic_coordinates = {}
	similarity_matrixs = {}
	evidences_similarity.each do |pair, ev_profiles_similarity|
		entity = pair.split('_').first
		similarities = ev_profiles_similarity[profile_id.to_sym]
		candidate_sim_matrix, candidates, candidates_ids = get_similarity_matrix(reference_prof, similarities, evidences[pair][:prof], hpo, 40, 40)
		candidate_sim_matrix.unshift(['HP'] + candidates_ids)	
		all_candidates.concat(candidates)
		similarity_matrixs[pair] = candidate_sim_matrix
		coords = get_evidence_coordinates(entity, genomic_coordinates, candidates_ids)
		all_genomic_coordinates.merge!(coords)
	end
	get_genome_hotspots(similarity_matrixs, all_genomic_coordinates)
	prof_vars = profile_variants[profile_id]
	make_report(profile_id, all_candidates, all_genomic_coordinates, similarity_matrixs, evidences, prof_vars, template, options[:output_folder])
end