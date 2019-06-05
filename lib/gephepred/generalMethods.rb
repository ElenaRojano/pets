require 'uri'
#Common methods for predictors
#Training file example = 9  131371492   131375954   HP:0010974  2.41161970596	9.3.A.5
#1. Indexing by chr (region)

def load_training_file4regions(training_file)
	training_set = {}
	posInfo = loadFile(training_file)
	posInfo.each do |info|
		chr = info.shift
		query = training_set[chr]
		if query.nil?
			training_set[chr] = [info]
		else
			query << info
		end
	end
	return training_set
end

#2. Indexing by hpo (code)
#prepare training file for analysis using phenotype2region prediction
def load_training_file4HPO(training_file, thresold=0)
	training_set = {}
	information = loadFile(training_file, thresold)
	information.each do |info|
		hpoCode = info.delete_at(4)
		query = training_set[hpoCode]
		if query.nil?
			training_set[hpoCode] = [info]  
		else
			query << info
		end
	end
	# STDERR.puts training_set.keys.inspect
	return training_set
end


#3. Load training info file:
#Chr;Start;Stop;HPO;Association;node
def loadFile(file, thresold=0)
	information = []
	File.open(file).each do |line|
		line.chomp!
		allInfo = line.split("\t")
		associationValue = allInfo[4].to_f
		if associationValue >= thresold
			chr = allInfo[0]
			startPos = allInfo[1].to_i
			stopPos = allInfo[2].to_i
			hpoCode = allInfo[3]
			nodeID = allInfo[5]
			information << [chr, startPos, stopPos, nodeID, hpoCode, associationValue]
		end
	end
	return information
end

def add_record2storage(hpo_storage, id, name, is_a, syn, alt_ids, hpo_black_list)
	if !hpo_black_list.include?(id)
		attributes = [id, name, is_a - hpo_black_list, syn]
		hpo_storage[id] = attributes
		alt_ids.each do |altid|
			hpo_storage[altid] = attributes
		end
	end 
end

def load_hpo_file(hpo_file, hpo_black_list=[])
	hpo_storage = {}
	id = nil
	name = nil
	alt_ids = []
	syn = []
	is_a = []
	File.open(hpo_file).each do |line|
		line.chomp!
		tag, info = line.split(': ')
		if tag == 'id' || tag == 'name' || tag == 'is_a' || tag == 'synonym' || tag == 'alt_id'
			if tag == 'id'
				add_record2storage(hpo_storage, id, name, is_a, syn, alt_ids, hpo_black_list) if !name.nil?
				id = info
				name = nil
				alt_id = []
				syn = []
				is_a = []
			end
			if tag == 'alt_id'
				alt_ids << info
			elsif tag == 'is_a'
				is_a << info.split(' ! ')[0]
			elsif tag == 'synonym'
				syn << info.split('"')[1] #to keep only the name of the synonym
			else
				name = info
			end
		end
	end
	add_record2storage(hpo_storage, id, name, is_a, syn, alt_ids, hpo_black_list)
	# STDERR.puts hpo_storage.inspect
	# Process.exit
	return hpo_storage
end

def load_hpo_black_list(excluded_hpo_file)
	excluded_hpos = []
	File.open(excluded_hpo_file).each do |line|
		line.chomp!
		excluded_hpos << line
	end
	return excluded_hpos
end

def create_hpo_dictionary(hpo_storage)
	hpo_dictionary = {}
	hpo_storage.each do |hpo, metadata|
		hpo_code, hpo_name, hpo_parents, hpo_synonyms = metadata
		hpo_dictionary[hpo_name] = hpo_code
		hpo_synonyms.each do |syn|
			hpo_dictionary[syn] = hpo_code
		end		
	end
	return hpo_dictionary
end

def get_child_parent_relations(hpo_storage)
	# for getting hpo childs
	storage_child = {}
	hpo_storage.each do |hpo_code, hpo_data|
		id, name, is_a, syn = hpo_data
		hpo_child = [id, name]
		is_a.each do |par_hpo_code|
			query = storage_child[par_hpo_code]
			if query.nil?
				storage_child[par_hpo_code] = [hpo_child]
			else
				query << hpo_child
			end
		end
	end
	return storage_child
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

# def get_child_parent_relations(hpo_storage)
# 	# for getting hpo childs
# 	storage_child = {}
# 	hpo_storage.each do |hpo_code, hpo_data|
# 		STDERR.puts hpo_data[3].inspect
# 		Process.exit
# 		main_code, hpo_name, synonyms, parents = hpo_data
# 		parents.each do |par_hpo_code, par_hpo_name|
# 			query = storage_child[par_hpo_code]
# 			hpo_child = [main_code, hpo_name]
# 			if query.nil?
# 				storage_child[par_hpo_code] = [par_hpo_name, [hpo_child]]
# 			else
# 				query.last << hpo_child
# 			end
# 		end
# 	end
	
# 	return storage_child
# end


def load_hpo_ci_values(information_coefficient_file)
	hpos_ci_values = {}
	File.open(information_coefficient_file).each do |line|
		line.chomp!
		hpo_code, ci = line.split("\t")
		hpos_ci_values[hpo_code] = ci.to_f
	end
	return hpos_ci_values
end

def load_clustered_patients(file)
	clusters = {}
	File.open(file).each do |line|
		line.chomp!
		pat_id, cluster_id = line.split("\t")
		query = clusters[cluster_id]
		if query.nil?
			clusters[cluster_id] = [pat_id]
		else
			query << pat_id
		end
	end
	return clusters
end

def load_gene_data(gene_data_path)
	gene_list = {} #geneID => attr
	gene_location = {} # chr => gene
	infile = open(gene_data_path)
	gz = Zlib::GzipReader.new(infile)
	current_chr = nil
	genes = []
	gz.each_line do |line|
		line.chomp!
		next if line =~ /^#/ 
		fields = line.split("\t")
		if fields[8].include?('genome=chromosome')
			chr = fields[8].split(';')[1].split('=').last
			gene_location[current_chr] = genes
			genes = []
			current_chr = chr
		elsif fields[2] == 'gene'
			attributes = {}
			fields[8].split(';').each do |pair| 
				key, value = pair.split('=')
				attributes[key] = value
			end
			geneNames = []
			geneNames << attributes['gene'] if !attributes['gene'].nil?
			geneNames.concat(attributes['gene_synonym'].split(',')) if !attributes['gene_synonym'].nil?
			description = attributes['description']
			description = URI.unescape(description) if !description.nil?
			attributes['Dbxref'] =~ /GeneID:(\d+)/
			gene_list[$1] = [geneNames, description]
			genes << [$1, fields[3].to_i, fields[4].to_i]
		end
	end
	gene_location[current_chr] = genes
	return gene_list, gene_location
end

def parse_kegg_data(query_genes)
	kegg_data = {} #gene => attb
    while !query_genes.empty?
    	gene_set = query_genes.shift(10)
	    url = "http://rest.kegg.jp/get/#{gene_set.map{|qg| "hsa:#{qg}"}.join('+')}"
	    uri = URI(url)
	    response = Net::HTTP.get(uri) 
		geneID = nil
		gene_names = []
		definition = nil
		pathways = []
		parsing_pathway_field = false
		response.squeeze(' ').each_line do |line|
			line.chomp!
			if line =~ /^ENTRY/
				geneID = line.split(' ')[1]
			elsif line =~ /^NAME/
				gene_names = line.split(' ', 2).last.split(', ')
			elsif line =~ /^DEFINITION/
				definition = line.split(' ', 2)[1]
			elsif line =~ /^PATHWAY/
				pathways << line.split(' ', 3)[1..2]
				parsing_pathway_field = true
			elsif line =~ /^BRITE/ || line =~ /^POSITION/ || line =~ /^DISEASE/ || line =~ /^MODULE/ || line =~ /^DRUG_TARGET/ || line =~ /^NETWORK/
				parsing_pathway_field = false
			elsif parsing_pathway_field
				pathways << line.strip.split(' ', 2)
			elsif line == '///'
				parsing_pathway_field = false
				kegg_data[geneID] = [gene_names, definition, pathways]
				pathways = []
				gene_names = []
			end
		end
	end
	return kegg_data
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
	gene_list.each do |geneID, attributes|
		query = kegg_data[geneID]
		if query.nil?
			attributes << []
		else
			attributes << query
		end
		merged_data[geneID] = attributes
	end
	return merged_data
end

def write_compressed_plain_file(data, path)
	File.open(path, 'w') do |f|
		gz = Zlib::GzipWriter.new(f)
		gz.write data.to_json
		gz.close
	end
end

def read_compressed_json(path)
	infile = open(path)
	gz = Zlib::GzipReader.new(infile)
	object = JSON.parse(gz.read)
	return object
end

def compute_pathway_enrichment(genes_clusters, genes_with_kegg)
	pathways_genes_in_predictions = {}
	genes_in_predictions = []
	genes_clusters.each do |cluster|
		cluster.each do |geneID, data|
			geneNames, description, pathways = data
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

def generate_cluster_regions(patients_genomic_region_by_chr, mutation_type, pat_per_reg = 1)
	patients_out_of_cluster = 0
	patients_by_cluster = {}
	sors = []
	patients_genomic_region_by_chr.each do |chrm, genomic_ranges|
		reference = get_reference(genomic_ranges) # Get putative overlap regions
		overlapping_patients = overlap_patients(genomic_ranges, reference) # See what patient has match with a overlap region
		clust_number = 1
		reference.each_with_index do |ref, i|
			current_patients = overlapping_patients[i]
			if current_patients.length > pat_per_reg
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

def add_record(hash, key, record)
	query = hash[key]
	if query.nil?
		hash[key] = [record]
	elsif !query.include?(record)
		query << record
	end
end

def load_patient_cohort(options)
	patient_data = {}
	count = 0
	fields2extract = get_fields2extract(options)
	field_numbers = fields2extract.values
  original_ids = []
  File.open(options[:input_file]).each do |line|
    line.chomp!
    if options[:header] && count == 0
      line.gsub!(/#\s*/,'') # correct comment like  headers
      field_names = line.split("\t")
      get_field_numbers2extract(field_names, fields2extract)
      field_numbers = fields2extract.values
    else
      fields = line.split("\t")
      pat_record = field_numbers.map{|n| fields[n]}
      if fields2extract[:pat_id_col].nil?
        pat_id = "pat_#{count}" #generate ids
      else
        original_id = pat_record.shift
        original_ids << original_id
        pat_id = original_id + "_i#{count}" # make sure that ids are uniq
      end
      patient_data[pat_id] = pat_record
    end
    count +=1
  end
  fields2extract[:pat_id_col].nil? ? patient_number = count : patient_number = original_ids.uniq.length
  options[:pat_id_col] = 'generated' if fields2extract[:pat_id_col].nil?
  return patient_data, patient_number
end 

def get_fields2extract(options)
	fields2extract = {}
	[:pat_id_col, :hpo_col, :chromosome_col, :start_col, :end_col].each do |field|
		col = options[field]
		if !col.nil?
			col = col.to_i if !options[:header]
			fields2extract[field] = col
		end
	end
	return fields2extract
end

def get_field_numbers2extract(field_names, fields2extract)
  fields2extract.each do |field, name|
    fields2extract[field] = field_names.index(name)
  end
end