require 'uri'
require 'net/ftp'
require 'net/http'
require 'zlib'
require 'json'
require 'benchmark'

def system_call(code_folder, script, args_string)
	cmd = File.join(code_folder, script) + ' ' + args_string
	puts "==> #{cmd}"
	Benchmark.bm do |x|
		x.report{	system(cmd) }
	end
end

def add_record(hash, key, record, uniq=false)
	query = hash[key]
	if query.nil?
		hash[key] = [record]
	elsif !uniq # We not take care by repeated entries
		query << record
	elsif !query.include?(record) # We want uniq entries
		query << record
	end
end


def compute_IC_values(patient_data, total_patients)
	patients_per_hpo = Hash.new(0)
	last_patient_ID = ''
	patient_data.each do |patient_ID, metadata|
		patient, count = patient_ID.split('_i')
		if patient != last_patient_ID
			hpos, chr, start, stop = metadata
			hpos.each do |h|
				patients_per_hpo[h] += 1
			end
		end
		last_patient_ID = patient
	end
	# STDERR.puts patients_per_hpo.inspect
	# Process.exit
	patients_per_hpo.each do |hpo, patient_number|
		patients_per_hpo[hpo] = -Math.log10(patient_number.fdiv(total_patients))
	end
	return patients_per_hpo
end


def parse_kegg_from_biosystems(biosystems_gene_path, biosystems_info_path)
	kegg_data = {}
	gene2biosystems = load_biosystem2gene_dictionary(biosystems_gene_path)
	keggAttributes = loadBiosistemsInfo(biosystems_info_path, 'KEGG')
	keggAttributes.select!{|kegg_id, data| data.first =~ /^hsa/}

	gene2biosystems.each do |geneID, pathways|
		kegg_pathways = []
		pathways.each do |biosystem|
			kAttrib = keggAttributes[biosystem]
			kegg_pathways << kAttrib if !kAttrib.nil?
		end
		kegg_data[geneID] = kegg_pathways
	end
	return kegg_data
end

def loadBiosistemsInfo(biosystems_info_path, filterDB)
	bsid2attributes = {}
	infile = open(biosystems_info_path)
	gz = Zlib::GzipReader.new(infile)
	gz.each_line do |line|
		line.chomp!
		#STDERR.puts line.inspect
		fields = line.encode('UTF-8', 'binary', invalid: :replace, undef: :replace, replace: '').split("\t")
		bsid = fields.shift
		bsid2attributes[bsid] = [fields[1], fields[2]] if filterDB == fields[0]
	end
	return bsid2attributes
end

def load_biosystem2gene_dictionary(biosystems_gene_path)
	gene2kegg = {}
	infile = open(biosystems_gene_path)
	gz = Zlib::GzipReader.new(infile)
	gz.each_line do |line|
		line.chomp!
		biosystem, gene, score = line.split("\t")
		query = gene2kegg[gene]
		if query.nil?
			gene2kegg[gene] = [biosystem]
		else
			query << biosystem
		end
	end
	return gene2kegg
end

def merge_genes_with_kegg_data(gene_list, kegg_data)
	merged_data = {}
	gene_list.each do |geneID, values|
		geneName, geneSyn, description = values
		kegg_entry = kegg_data[geneID]	
		kegg_entry = [] if kegg_entry.nil?
		merged_data[geneID] = [geneName, description, kegg_entry, geneSyn]
	end
	return merged_data
end

def compute_pathway_enrichment(genes_clusters, genes_with_kegg)
	pathways_genes_in_predictions = {}
	genes_in_predictions = []
	genes_clusters.each do |cluster|
		cluster.each do |geneID, data|
			geneName, description, pathways, geneSyns = data
			pathways.each do |pathway|
				query = pathways_genes_in_predictions[pathway]
				if query.nil?
					pathways_genes_in_predictions[pathway] = [geneID]
				else
					query << geneID if !query.include?(geneID)
				end
			end
			genes_in_predictions << geneID if !genes_in_predictions.include?(geneID)
		end
	end
	genes_out_of_predictions = genes_with_kegg.keys - genes_in_predictions
	gene_number = genes_with_kegg.length
	stats = []
	pathways_genes_in_predictions.each do |pathway, pathway_predicted_genes|
		pathway_id, pathway_name = pathway
		no_pathway_predicted_genes = genes_in_predictions - pathway_predicted_genes 
		pathway_no_predicted_genes_count = 0
		no_pathway_no_predicted_genes_count = 0
		genes_out_of_predictions.each do |geneID|
			query = genes_with_kegg[geneID]
			if query[2].map{|pathway_info| pathway_info.first}.include?(pathway_id)
				pathway_no_predicted_genes_count += 1
			else
				no_pathway_no_predicted_genes_count += 1
			end
		end
		#Fisher => http://www.biostathandbook.com/fishers.html
		no_pathway_predicted_genes_count = no_pathway_predicted_genes.length
		pathway_predicted_genes_count = pathway_predicted_genes.length
		accumulated_prob = 0
		pathway_no_predicted_genes_count.times do |n|
			no_pathway_predicted_genes_count_shifted = no_pathway_predicted_genes_count - n
			pathway_predicted_genes_count_shifted = pathway_predicted_genes_count - n
			if no_pathway_predicted_genes_count_shifted >= 0 && pathway_predicted_genes_count_shifted >= 0
				accumulated_prob += compute_hyper_prob(
					n, 
					no_pathway_predicted_genes_count_shifted, 
					pathway_predicted_genes_count_shifted, 
					no_pathway_no_predicted_genes_count + n, 
					gene_number
				)
			else
				break
			end
		end
		contigency = [pathway_no_predicted_genes_count, no_pathway_predicted_genes_count, pathway_predicted_genes_count, no_pathway_no_predicted_genes_count]
		stats << [pathway, pathway_predicted_genes, contigency, accumulated_prob]
	end
	return stats
end

def compute_hyper_prob(a, b, c, d, n)
	binomA = binom(a + b, a)
	binomC = binom(c + d, c)
	divisor = binom(n, a + c)
	return (binomA * binomC).fdiv(divisor)
end

def binom(n,k)
	if k > 0 && k < n
		res = (1+n-k..n).inject(:*)/(1..k).inject(:*)
	else
		res = 1
	end
end




def get_and_parse_external_data(all_paths)
	sources = [
	  ['ftp.ncbi.nlm.nih.gov', 'genomes/H_sapiens/ARCHIVE/ANNOTATION_RELEASE.105/GFF/ref_GRCh37.p13_top_level.gff3.gz', all_paths[:gene_data]],
	  ['ftp.ncbi.nlm.nih.gov', 'pub/biosystems/CURRENT/biosystems_gene.gz', all_paths[:biosystems_gene]],
	  ['ftp.ncbi.nlm.nih.gov', 'pub/biosystems/CURRENT/bsid2info.gz', all_paths[:biosystems_info]]
	]
	sources.each do |server, path, output|
	  download(server, path, output) if !File.exists?(output)
	end

	genes_with_kegg = {}
	gene_location = {}
	if !File.exists?(all_paths[:gene_data_with_pathways]) || !File.exists?(all_paths[:gene_location])
	  gene_list, gene_location = load_gene_data(all_paths[:gene_data])
	  ### kegg_data = parse_kegg_data(genes_found_attributes.keys)
	  kegg_data = parse_kegg_from_biosystems(all_paths[:biosystems_gene], all_paths[:biosystems_info])
	  genes_with_kegg = merge_genes_with_kegg_data(gene_list, kegg_data)
	  write_compressed_plain_file(genes_with_kegg, all_paths[:gene_data_with_pathways])
	  write_compressed_plain_file(gene_location, all_paths[:gene_location])
	else
	  gene_location = read_compressed_json(all_paths[:gene_location])
	  genes_with_kegg = read_compressed_json(all_paths[:gene_data_with_pathways])
	end
	return gene_location, genes_with_kegg
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

def get_similarity_matrix(reference_prof, similarities, evidence_profiles, hpo, term_limit, candidate_limit, other_scores = {}, id2label = {})
		candidates = similarities.to_a
		if other_scores.empty?
			candidates.sort!{|s1, s2| s2.last <=> s1.last}
			candidates = candidates.first(candidate_limit)
		else # Prioritize first by the external list of scores, select the candidates and then rioritize by similarities
			selected_candidates = []
			candidates.each do |cand|
				cand_id = cand[0]
				cand_lab = id2label[cand_id.to_s]
				next if !cand_lab.nil?
				other_score = other_scores[cand_lab]
				next if !other_score.nil?
				cand << other_score
				selected_candidates << cand
			end
			selected_candidates.sort!{|e1, e2| e2[2] <=> e1[2]}
			candidates = selected_candidates.first(candidate_limit)
			candidates.sort!{|e1, e2| e2[1] <=> e1[1]}
			candidates.each{|c| c.pop}
		end
		candidates_ids = candidates.map{|c| c.first}
		candidate_similarity_matrix = get_detailed_similarity(reference_prof, candidates, evidence_profiles, hpo)
		candidate_similarity_matrix.each_with_index do |row, i|
			row.unshift(hpo.translate_id(reference_prof[i]))
		end
		candidate_similarity_matrix.sort!{|r1,r2| r2[1..r2.length].inject(0){|sum,n| sum +n} <=> r1[1..r1.length].inject(0){|sum,n| sum +n}}
		candidate_similarity_matrix = candidate_similarity_matrix.first(term_limit)
		return candidate_similarity_matrix, candidates, candidates_ids
end