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
def load_pathogenic_scores(path)
	scores = {}
	File.open(path).each do |line|
		feature, score = line.split("\t")
		scores[feature] = score.to_f
	end
	return scores
end

def get_evidence_coordinates(entity, genomic_coordinates, candidates_ids)
	coords = nil
	all_coordinates = genomic_coordinates[entity]
	coords = all_coordinates.select{|id, coordinates| candidates_ids.include?(id.to_sym)} if !all_coordinates.nil?
	return coords
end

def make_report(profile_id, all_candidates, all_genomic_coordinates, similarity_matrixs, 
							evidences, prof_vars, hotspots_with_pat_vars, template, output)
	var_ids, var_coors = format_variants4report(prof_vars)
	container = {
		profile_id: profile_id,
		candidates: all_candidates.each{|c| c[0] = c.first.to_s},
		genomic_coordinates: all_genomic_coordinates.transform_values{|c| c.first(2) },
		similarity_matrixs: similarity_matrixs,
		evidences: evidences,
		var_ids: var_ids,
		var_coordinates: var_coors,
		hotspot_table: hotspots_with_pat_vars
	}
	report = Report_html.new(container, 'Evidence profile report')
	report.build(template)
	report.write(File.join(output, profile_id.to_s + '.html'))
end

def format_variants4report(var_data)
	if var_data.nil?
		var_ids, var_coors = nil
	else 
		var_ids = []
		var_coors = {}
		count = 0
		var_data.each do |chr, reg| 
			var_id = "var_#{count}"
			var_ids << [var_id, 0]
			var_coors[var_id] = [chr.to_s, reg[:start]]
			count += 1
		end
	end
	return var_ids, var_coors
end


def generate_prediction(similarity_matrixs, all_genomic_coordinates, prof_vars)
	hotspots_with_pat_vars = []
	if !prof_vars.nil?
		phen_regions = Genomic_Feature.hash2genomic_feature(all_genomic_coordinates){|k, v| v[0..2].concat([k])}
		phen_candidates_by_hotspot, phen_genome_hotspots = phen_regions.generate_cluster_regions(:reg_overlap, 'A', 0, true)
		genome_matches = phen_genome_hotspots.match(prof_vars)
		hotspot_with_phen_candidates = invert_hash(phen_candidates_by_hotspot) 
		genome_matches.each do |hotspot_id, pat_vars|
			reg = phen_genome_hotspots.region_by_to(hotspot_id)
			coords = [reg[:chr], reg[:start], reg[:stop]]
			hotspots_with_pat_vars << [hotspot_id, coords, hotspot_with_phen_candidates[hotspot_id], pat_vars]
		end
		# TODO: see to use original similarities without use top candidates in similarity_matrixs
		# TODO: COMPLETE UNTIL FULL PREDICTOR
	end
	return hotspots_with_pat_vars
end

def invert_hash(h)
	new_h = {}
	h.each do |k, vals|
		vals.each do |v|
			query = new_h[v]
			if query.nil?
				new_h[v] = [k]
			else
				query << k
			end
		end
	end
	return new_h
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

  options[:pathogenic_scores] = nil # TODO: Generalize to a folder with a table per patient
  opts.on("-P", "--pathogenic_scores PATH", 'File with genome features an their pathogenic scores') do |item|
    options[:pathogenic_scores] = item
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
hpo = Ontology.new(file: hpo_file, load_file: true)

profiles = load_profiles(options[:profiles_file], hpo)
profile_variants = options[:variant_data].nil? ? {} : load_variants(options[:variant_data])
evidences, genomic_coordinates = load_evidences(options[:evidences], hpo)
pathogenic_scores = options[:pathogenic_scores].nil? ? {} : load_pathogenic_scores(options[:pathogenic_scores])

hpo.load_profiles(profiles)
evidences_similarity = {}
evidences.each do |pair, data|
	entity, profile_type = pair.split('_')
	if profile_type == 'HP'
		evidence_profiles = data[:prof]
		evidence_profiles.transform_keys!{|prof_id, terms| prof_id.to_sym}
		similarities = hpo.compare_profiles(external_profiles: evidence_profiles, sim_type: :lin, bidirectional: false)
		evidences_similarity[pair] = similarities if !similarities.empty?
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
		coords = get_evidence_coordinates(entity, genomic_coordinates, candidates_ids)
		candidate_sim_matrix.unshift(['HP'] + candidates_ids)
		if !pathogenic_scores.empty? # priorize by pathogenic scores
			candidate_sim_matrix_patho, candidates_patho, candidates_ids_patho = get_similarity_matrix(
				reference_prof, similarities, 
				evidences[pair][:prof], hpo, 40, 40, 
				other_scores = pathogenic_scores, id2label = evidences[pair][:id2lab])
			if !candidate_sim_matrix_patho.empty?
				candidate_sim_matrix_patho.unshift(['HP'] + candidates_ids_patho)	
				similarity_matrixs[pair + '_path_vars'] = candidate_sim_matrix_patho
				evidences[pair + '_path_vars'] = evidences[pair]
			end
		end
		next if coords.nil?
		all_candidates.concat(candidates)
		similarity_matrixs[pair] = candidate_sim_matrix
		all_genomic_coordinates.merge!(coords)
	end
	prof_vars = profile_variants[profile_id]
	hotspots_with_pat_vars = generate_prediction(similarity_matrixs, all_genomic_coordinates, prof_vars)
	make_report(
		profile_id, 
		all_candidates, 
		all_genomic_coordinates, 
		similarity_matrixs, 
		evidences, prof_vars,
		hotspots_with_pat_vars,
		template, options[:output_folder]
	)
end