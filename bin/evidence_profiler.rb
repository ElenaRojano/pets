#! /usr/bin/env ruby

require 'fileutils'
require 'optparse'
require 'report_html'
require 'semtools'

ROOT_PATH = File.dirname(__FILE__)
REPORT_FOLDER = File.expand_path(File.join(ROOT_PATH, '..', 'templates'))
EXTERNAL_DATA = File.expand_path(File.join(ROOT_PATH, '..', 'external_data'))
EXTERNAL_CODE = File.expand_path(File.join(ROOT_PATH, '..', 'external_code'))
HPO_FILE = File.join(EXTERNAL_DATA, 'hp.json')

class Report_html 
	def circular_genome(user_options = {}, &block)
		default_options = {}.merge!(user_options)
		coordinates = user_options[:genomic_coordinates]
		html_string = canvasXpress_main(default_options, block) do |options, config, samples, vars, values, object_id, x, z|
			config['graphType'] = 'Circular'
			config["arcSegmentsSeparation"] = 3
        	config["colorScheme"] = "Tableau"
        	config["colors"] = ["#332288","#6699CC","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#661100","#CC6677","#AA4466","#882255","#AA4499"]
			config["showIdeogram"] = true
			chr = []
			pos = []
			tags2remove = []
			vars.each_with_index do |var, i|
				coord = coordinates[var]
				if !coord.nil?
					tag = coord.first.gsub(/[^\dXY]/,'')
					if tag == 'X' || tag == 'Y' || (tag.to_i > 0 && tag.to_i <= 22)
						chr << coord.first.gsub(/[^\dXY]/,'')
						pos << coord.last - 1
					else
						tags2remove << i
					end
				else
					tags2remove << i
				end
			end
			tags2remove.reverse_each{|i| ent = vars.delete_at(i); warn("Feature #{ent} has not valid coordinates")} # Remove entities with invalid coordinates
			z['chr'] = chr
			z['pos'] = pos
		end
		return html_string
	end
end

#############################################################################################
## METHODS
############################################################################################
def load_profiles(file_path, hpo)
	profiles = {}
	#count = 0
	File.open(file_path).each do |line|
		id, profile = line.chomp.split("\t")
		hpos = profile.split(',').map{|a| a.to_sym}
		hpos, rejected_hpos = hpo.check_ids(hpos)
		if !hpos.empty?
			hpos = hpo.clean_profile(hpos)
			profiles[id] = hpos if !hpos.empty?
		end
	end
	return profiles
end

def load_coordinates(file_path)
	coordinates = {}
	File.open(file_path).each do |line|
		id, chr, start = line.chomp.split("\t")
		coordinates[id] = [chr, start.to_i]
	end
	return coordinates
end

def get_detailed_similarity(profile, candidates, evidences, hpo)
	profile_length = profile.length
	matrix = []
	profile_length.times{ matrix << Array.new(candidates.length, 0)}
	cand_number = 0
	candidates.each do |candidate_id, similarity|
		local_sim = []
		candidate_evidence = evidences[candidate_id]
		profile.each do |profile_term|
			candidate_evidence.each do |candidate_term|
				term_sim = hpo.compare([candidate_term], [profile_term], sim_type: :lin, bidirectional: false)
				local_sim << [profile_term, candidate_term, term_sim]
			end
		end
		local_sim.sort!{|s1, s2| s2.last <=> s1.last}
		final_pairs = []
		processed_profile_terms = []
		processed_candidate_terms = []
		local_sim.each do |pr_term, cd_term, sim|
			if !processed_profile_terms.include?(pr_term) && !processed_candidate_terms.include?(cd_term)
				final_pairs << [pr_term, cd_term, sim]
				processed_profile_terms << pr_term
				processed_candidate_terms << cd_term
			end
			break if profile_length == processed_profile_terms.length
		end
		final_pairs.each do |pr_term, cd_term, similarity|
			matrix[profile.index(pr_term)][cand_number] = similarity
		end
		cand_number += 1
	end
	return matrix
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

  options[:output_folder] = 'evidence_reports'
  opts.on("-o", "--output_folder PATH", 'Folder to save reports from profiles') do |item|
    options[:output_folder] = item
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

genomic_coordinates = load_coordinates(options[:coordinates_file])
profiles = load_profiles(options[:profiles_file], hpo)
evidences = load_profiles(options[:evidence_file], hpo)

hpo.load_profiles(profiles)
FileUtils.mkdir_p(options[:output_folder])
template = File.open(File.join(REPORT_FOLDER, 'evidence_profile.erb')).read

profiles_similarity = hpo.compare_profiles(external_profiles: evidences, sim_type: :lin, bidirectional: false)
profiles_similarity.each do |profile_id, similarities|
	candidates = similarities.to_a.sort{|s1, s2| s2.last <=> s1.last}.first(40)
	candidates_ids = candidates.map{|c| c.first}
	profile = profiles[profile_id.to_s] #compare_profiles return string id when al data is loaded with ids as symbols
	candidate_similarity_matrix = get_detailed_similarity(profile, candidates, evidences, hpo)
	candidate_similarity_matrix.each_with_index do |row, i|
		row.unshift(hpo.translate_id(profile[i]))
	end
	candidate_similarity_matrix.unshift(['HP'] + candidates_ids)
	container = {
		profile_id: profile_id,
		similarity_matrix: candidate_similarity_matrix,
		genomic_coordinates: genomic_coordinates.select{|id, coordinates| candidates_ids.include?(id)},
		candidates: candidates
	}
	report = Report_html.new(container, 'Evidence profile report')
	report.build(template)
	report.write(File.join(options[:output_folder], profile_id.to_s + '.html'))
end
