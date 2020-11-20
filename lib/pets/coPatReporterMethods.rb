require 'numo/narray'

HPOS = 0
CHR = 1
START = 2
STOP = 3

def format_patient_data(patient_data, options, hpo)
  all_hpo = []
  rejected_hpos = []
  suggested_childs = {}
  total_terms = 0
  terms_with_more_specific_childs = 0
  patient_data.each do |pat_id, patient_record|
    hpos, chr, start, stop = patient_record

    if options[:hpo_names]
      # hpos, pat_rejected_hpos = hpo.translate_names2codes(hpos)
      hpos, pat_rejected_hpos = hpo.translate_names(hpos)
      if !pat_rejected_hpos.empty?
        STDERR.puts "WARNING: patient #{pat_id} has the unknown hpo NAMES '#{pat_rejected_hpos.join(',')}'. Rejected."
        rejected_hpos.concat(pat_rejected_hpos)
      end
    end

    # hpos, pat_rejected_hpos = hpo.check_codes(hpos)
    hpos, pat_rejected_hpos = hpo.check_ids(hpos.map{|a| a.to_sym})
    if !pat_rejected_hpos.empty?
      STDERR.puts "WARNING: patient #{pat_id} has the unknown hpo CODES '#{pat_rejected_hpos.join(',')}'. Rejected."
      rejected_hpos.concat(pat_rejected_hpos)
    end
    total_terms += hpos.length
    # more_specific_childs = hpo.get_more_specific_childs_table(hpos)
    more_specific_childs = hpo.get_childs_table(hpos, true)
    terms_with_more_specific_childs += more_specific_childs.select{|hpo_record| !hpo_record.last.empty?}.length #Exclude phenotypes with no childs
    suggested_childs[pat_id] = more_specific_childs  
    all_hpo.concat(hpos)
    patient_record[HPOS] = hpos
  end
  return all_hpo.uniq, suggested_childs, rejected_hpos.uniq, terms_with_more_specific_childs.fdiv(total_terms)
end

def generate_patient_hpo_matrix(patient_data, cohort_hpos)
  matrix = []
  n = cohort_hpos.length
  patient_data.each do |pat_id, patient_record|
    pat_hpos = patient_record[HPOS]
    vector = Array.new(n, 0)
    pat_hpos.each do |hpo|
      vector[cohort_hpos.index(hpo)] = 1
    end
    matrix << vector
  end
  return matrix
end

def write_matrix_for_R(matrix, x_names, y_names, file)
  File.open(file, 'w') do |f|
    f.puts x_names.join("\t")
    matrix.each_with_index do |row, i|
      f.puts [y_names[i]].concat(row).join("\t")
    end
  end
end

def process_clustered_patients(options, clustered_patients, patient_data, hpo, onto_ic, freq_ic, patient_id_type) # get ic and chromosomes
  if options[:ic_stats] == 'freq_internal'
    ic_file = ENV['ic_file']
    ic_file = IC_FILE if ic_file.nil?
    phenotype_ic = load_hpo_ci_values(ic_file)
  elsif options[:ic_stats] == 'freq'
    phenotype_ic = freq_ic
  elsif options[:ic_stats] == 'onto'
    phenotype_ic = onto_ic
  end
  all_ics = []
  top_cluster_phenotypes = []
  cluster_data_by_chromosomes = []
  multi_chromosome_patients = 0
  processed_clusters = 0
  clustered_patients.sort_by{|cl_id, pat_ids| pat_ids.length }.reverse.each do |cluster_id, patient_ids|
    next if patient_ids.length == 1
    chrs = Hash.new(0)
    all_phens = []
    profile_ics = []
    processed_patients = []
    patient_ids.each do |pat_id|
      patient = patient_data[pat_id]  
      pat_id = pat_id.gsub(/_i\d+$/,'') if patient_id_type != 'generated'
      if !processed_patients.include?(pat_id) # Check that current cluster member is not an additional mutation of a previous patient
        processed_patients << pat_id 
        phenotypes = patient[HPOS]
        profile_ics << get_profile_ic(phenotypes, phenotype_ic)
        if processed_clusters < options[:clusters2show_detailed_phen_data]
          # phen_names, rejected_codes = hpo.translate_codes2names(phenotypes) #optional
          phen_names, rejected_codes = hpo.translate_ids(phenotypes) #optional
          all_phens << phen_names
        end
      end
      chrs[patient[CHR]] += 1 if !options[:chromosome_col].nil?
    end
    num_of_patients = processed_patients.length
    next if num_of_patients == 1 # Check that current cluster only has one patient with several mutations
    top_cluster_phenotypes << all_phens if processed_clusters < options[:clusters2show_detailed_phen_data]
    all_ics << profile_ics
    if !options[:chromosome_col].nil?
      multi_chromosome_patients += num_of_patients if chrs.length > 1
      chrs.each do |chr, count|
        cluster_data_by_chromosomes << [cluster_id, num_of_patients, chr, count]
      end
    end
    processed_clusters += 1
  end
  # STDERR.puts cluster_data_by_chromosomes.inspect
  return all_ics, cluster_data_by_chromosomes, top_cluster_phenotypes, multi_chromosome_patients
end

def get_profile_ic(hpo_names, phenotype_ic)
  ic = 0
  profile_length = 0
  hpo_names.each do |hpo_id|
    hpo_ic = phenotype_ic[hpo_id]
    # STDERR.puts phenotype_ic.inspect
    ic += hpo_ic if !hpo_ic.nil?
    profile_length += 1
  end
  profile_length = 1 if profile_length == 0
  return ic.fdiv(profile_length)
end

def write_cluster_ic_data(all_ics, cluster_ic_data_file, limit)
  File.open(cluster_ic_data_file, 'w') do |f|
    f.puts %w[cluster_id ic].join("\t")
    all_ics.each_with_index do |cluster_ics, i|
      break if i == limit
      cluster_length = cluster_ics.length
      cluster_ics.each do |clust_ic|
        f.puts "#{cluster_length}_#{i}\t#{clust_ic}"
      end
    end
  end
end

def write_cluster_chromosome_data(cluster_data, cluster_chromosome_data_file, limit)
  File.open(cluster_chromosome_data_file, 'w') do |f|
    f.puts %w[cluster_id chr count].join("\t")
    index = 0
    last_id = cluster_data.first.first unless cluster_data.empty?
    cluster_data.each do |cluster_id, patient_number, chr, count|
      index += 1 if cluster_id != last_id 
      break if index == limit
      f.puts ["#{patient_number}_#{index}", chr, count].join("\t")
      last_id = cluster_id
    end
  end
end

def write_coverage_data(coverage_to_plot, coverage_to_plot_file)
  File.open(coverage_to_plot_file, 'w') do |f|
    coverage_to_plot.each do |chr, position, freq|
     f.puts "#{chr}\t#{position}\t#{freq}"
   end
  end
end

def get_uniq_hpo_profiles(patient_data) # To avoid duplications due to more one mutation in the same patient
  #STDERR.puts patient_data.keys.inspect
  #Process.exit
  #transformar a hash
  hpo_profiles = {}
  #hpo_profiles = []
  parsed_pats = []
  patient_data.each do |variant_id, patient_rec|
    pat_id, count = variant_id.split('_i')
    if parsed_pats.include?(pat_id)
      next
    else
      parsed_pats << pat_id
      #hpo_profiles << patient_rec[HPOS]
      hpo_profiles[pat_id] = patient_rec[HPOS]
    end
  end
  return hpo_profiles
end

def get_patient_ids(patient_data) # To aviod duplications due to more one mutation in the same patient
  ids = []
  patient_ids = patient_data.keys
  patient_ids.each do |pat_id|
    id, count = pat_id.split('_i')
    ids << id
  end
  return  ids.uniq
end

def get_summary_stats(patient_data, cohort_hpos, hpo)
  stats = []
  stats << ['Unique HPOs', cohort_hpos.length]
  stats << ['Number of patients in the cohort', get_patient_ids(patient_data).length]
  # stats << ['HPOs per patient (average)', hpo.get_profile_mean_length]
  stats << ['HPOs per patient (average)', hpo.get_profiles_mean_size]
  stats << ['HPOs for patient in percentile 90', hpo.get_profile_length_at_percentile(perc=90)]
  return stats
end

def cluster_patients(patient_data, cohort_hpos, matrix_file, clustered_patients_file)
  pat_hpo_matrix = generate_patient_hpo_matrix(patient_data, cohort_hpos)
  write_matrix_for_R(pat_hpo_matrix, cohort_hpos, patient_data.keys, matrix_file)
  system("#{File.join(EXTERNAL_CODE, 'get_clusters.R')} #{matrix_file} #{clustered_patients_file}") if !File.exists?(clustered_patients_file)
  clustered_patients = load_clustered_patients(clustered_patients_file)
  return(clustered_patients)
end

def get_profile_ontology_distribution_tables(hpo)
  cohort_ontology_levels = hpo.get_ontology_levels_from_profiles(uniq=false)
  uniq_cohort_ontology_levels = hpo.get_ontology_levels_from_profiles
  # hpo_ontology_levels = hpo.term_level
  hpo_ontology_levels = hpo.get_ontology_levels
  total_ontology_terms = hpo_ontology_levels.values.flatten.length
  total_cohort_terms = cohort_ontology_levels.values.flatten.length
  total_uniq_cohort_terms = uniq_cohort_ontology_levels.values.flatten.length

  ontology_levels = []
  distribution_percentage = []
  hpo_ontology_levels.each do |level, terms|
    cohort_terms = cohort_ontology_levels[level]
    uniq_cohort_terms = uniq_cohort_ontology_levels[level]
    if cohort_terms.nil? || uniq_cohort_terms.nil?
      num = 0
      u_num = 0
    else
      num = cohort_terms.length
      u_num = uniq_cohort_terms.length
    end
    ontology_levels << [level, terms.length, num]
    distribution_percentage << [
      level,
      (terms.length.fdiv(total_ontology_terms)*100).round(3),
      (num.fdiv(total_cohort_terms)*100).round(3),
      (u_num.fdiv(total_uniq_cohort_terms)*100).round(3)
    ]
  end
  ontology_levels.sort! { |x,y| x.first <=> y.first }
  distribution_percentage.sort! { |x,y| x.first <=> y.first }
  ontology_levels.unshift ["level", "ontology", "cohort"]
  distribution_percentage.unshift ["level", "ontology", "weighted cohort", "uniq terms cohort"]
  return ontology_levels, distribution_percentage
end


def write_detailed_hpo_profile_evaluation(suggested_childs, detailed_profile_evaluation_file, summary_stats)
  CSV.open(detailed_profile_evaluation_file, "wb") do |csv|
    suggested_childs.each do |pat_id, suggestions|
      warning = nil
      warning = 'WARNING: Very few phenotypes' if suggestions.length < 4
      csv << ["PATIENT #{pat_id}", "#{warning}"]
      csv << ["CURRENT PHENOTYPES", "PUTATIVE MORE SPECIFIC PHENOTYPES"]
      suggestions.each do |parent, childs|
        parent_code, parent_name = parent
        if childs.empty?
          csv << ["#{parent_name} (#{parent_code})", '-']
        else
          parent_writed = false
          childs.each do |child_code, child_name|
            if !parent_writed
              parent_field = "#{parent_name} (#{parent_code})"
              parent_writed = true
            else
              parent_field = ""
            end
            csv << [parent_field, "#{child_name} (#{child_code})"]
          end
        end
      end
      csv << ["", ""]
    end
  end
end

def write_arrays4scatterplot(x_axis_value, y_axis_value, filename, x_axis_name, y_axis_name)
  File.open(filename, 'w') do |f|
    f.puts "#{x_axis_name}\t#{y_axis_name}"
    x_axis_value.each_with_index do |value,i|
        f.puts [value, y_axis_value[i]].join("\t")
    end
  end
end

def process_patient_data(patient_data)
	parsed_patient_data = {}
	patient_data.each do |patientID, metadata|
		phenotypes, chr, start, stop = metadata
		info = [patientID, start.to_i, stop.to_i]
		query = parsed_patient_data[chr]
		if query.nil?
			parsed_patient_data[chr] = [info]
		else
			query << info
		end
	end
	return parsed_patient_data
end

def get_final_coverage(raw_coverage, bin_size)
	coverage_to_plot = []
	raw_coverage.each do |chr, coverages|
		coverages.each do |start, stop, coverage|
			bin_start = start - start % bin_size
			bin_stop = stop - stop%bin_size
			while bin_start < bin_stop
				coverage_to_plot << [chr, bin_start, coverage]
				bin_start += bin_size
			end
		end
	end
	return coverage_to_plot
end

def get_sor_length_distribution(raw_coverage)
	all_cnvs_length = []
	cnvs_count = []
	raw_coverage.each do |chr, coords_info|
		coords_info.each do |start, stop, pat_records|
			region_length = stop - start + 1
			all_cnvs_length << [region_length, pat_records]
		end
	end
	all_cnvs_length.sort!{|c1, c2| c1[1] <=> c2[1]}
	return all_cnvs_length
end

def get_cnvs_length(patient_data)
	length_stats = Hash.new(0)
	patient_data.each do |pat_id, patient_record|
    	string_hpos, chr, start, stop = patient_record
    	length_stats[stop - start] += 1
    end
    return length_stats.to_a.sort!{|stat| stat[1] <=> stat[1] }
end


def calculate_coverage(regions_data, delete_thresold = 0)
	raw_coverage = {}
	n_regions = 0
	patients = 0
	nt = 0
	regions_data.each do |start, stop, chr, node|
		number_of_patients = node.split('.').last.to_i
		if number_of_patients <= delete_thresold
			number_of_patients = 0
		else
			n_regions += 1
			patients += number_of_patients
			nt += stop - start			
		end
		coords = [start, stop, number_of_patients]
		query = raw_coverage[chr]
		if query.nil?
			raw_coverage[chr] = [coords]
		else
			query << coords
		end
	end
	return raw_coverage, n_regions, nt, patients.fdiv(n_regions)
end

def get_profile_redundancy(hpo)
  #TODO: sort both arrays consequently
  #TODO: bear in mind join both arrays with zip and sort by one, to get an [a[a]]
  # profile_sizes = hpo.get_profile_sizes
  profile_sizes = hpo.get_profiles_sizes
  # parental_hpos_per_profile = hpo.compute_redundant_parental_terms_per_profile
  parental_hpos_per_profile = hpo.parentals_per_profile# clean_profiles
  parental_hpos_per_profile = parental_hpos_per_profile.map{|item| item[0]}
  profile_sizes, parental_hpos_per_profile = profile_sizes.zip(parental_hpos_per_profile).sort_by{|i| i.first}.reverse.transpose
  return profile_sizes, parental_hpos_per_profile
end

def format_profiles_similarity_data(profiles_similarity)
  matrix = []
  element_names = profiles_similarity.keys
  matrix << element_names
  profiles_similarity.each do |elementA, relations|
    row = [elementA]
    element_names.each do |elementB|
      if elementA == elementB
        row << 'NA'
      else
        query = relations[elementB]
        if !query.nil?
          row << query
        else
          row << profiles_similarity[elementB][elementA]
        end
      end
    end
    matrix << row
  end
  matrix[0].unshift('pat')
  return matrix
end

def format_profiles_similarity_data_pairs(profiles_similarity)
  pairs = []
  element_names = profiles_similarity.keys
  profiles_similarity.each do |elementA, relations|
    element_names.each do |elementB|
      if elementA != elementB
        pair = [elementA, elementB]
        query = relations[elementB]
        if !query.nil?
          pair << query
        else
          pair << profiles_similarity[elementB][elementA]
        end
        pairs << pair
      end
    end
  end
  return pairs
end

def format_profiles_similarity_data_numo(profiles_similarity)
  element_names = profiles_similarity.keys
  matrix = Numo::DFloat.zeros(element_names.length, element_names.length)
  i = 0
  profiles_similarity.each do |elementA, relations|
    element_names.each_with_index do |elementB, j|
      if elementA != elementB
        query = relations[elementB]
        if !query.nil?
          matrix[i, j] = query
        else
          matrix[i, j] = profiles_similarity[elementB][elementA]
        end
      end
    end
    i += 1
  end
  return matrix, element_names
end

def write_similarity_matrix(similarity_matrix, similarity_matrix_file)  
  File.open(similarity_matrix_file, 'w') do |f|
    similarity_matrix.each do |row|
      f.puts row.join("\t")
    end
  end
end

def write_profile_pairs(similarity_pairs, filename)
  File.open(filename, 'w') do |f|
    similarity_pairs.each do |pairsA, pairsB_and_values|
      pairsB_and_values.each do |pairsB, values|
        f.puts "#{pairsA}\t#{pairsB}\t#{values}"
      end
    end
  end
end

def parse_clusters_file(clusters_file, patient_profiles)
  clusters_info = {}
  clusters_table = []
  File.open(clusters_file).each do |line|
    line.chomp!
    patientID, clusterID = line.split("\t")
    patientHPOProfile = patient_profiles[patientID]
    query = clusters_info[clusterID]
    if query.nil? 
      clusters_info[clusterID] = {patientID => patientHPOProfile}
    else
      query[patientID] = patientHPOProfile
    end
  end
  clusters_info.each do |clusterID, patients_info|
    patients_per_cluster = patients_info.keys.length
    clusters_table << [clusterID, patients_per_cluster, patients_info.keys, patients_info.values]
  end
  return clusters_table, clusters_info
end

def get_patient_hpo_frequency(patient_uniq_profiles, hpo_frequency_file)
  hpo_frequency = Hash.new(0)
  patient_uniq_profiles.values.each do |hpos|
    hpos.each do |hpo|
      hpo_frequency[hpo] += 1
    end
  end
  File.open(hpo_frequency_file, 'w') do |f|
    hpo_frequency.each do |hpo_code, freq|
      f.puts "#{hpo_code.to_s}\t#{freq}"
    end
  end
end

def get_cluster_metadata(clusters_info, output_file)
  distribution_average_phenotypes_by_patients_cluster = []
  tmp = []
  total_clusters = clusters_info.keys.length
  average_patients_per_cluster = clusters_info.values.map{|a| a.keys.length}.inject(0){|i, sum| i + sum}.fdiv(total_clusters)
  clusters_info.values.each do |i|
      patients_by_cluster = i.keys.length
      num_phenotypes_in_cluster_by_patients_ary = i.values.map{|a| a.length}
      num_phenotypes_in_cluster_by_patients_ave = num_phenotypes_in_cluster_by_patients_ary.sum.fdiv(num_phenotypes_in_cluster_by_patients_ary.length)
      distribution_average_phenotypes_by_patients_cluster << [patients_by_cluster, num_phenotypes_in_cluster_by_patients_ave]
      tmp << num_phenotypes_in_cluster_by_patients_ary
  end
  average_phenotypes_by_cluster = tmp.flatten.sum.fdiv(total_clusters)
  File.open(output_file, 'w') do |f|
    f.puts "#{'PatientsNumber'}\t#{'HPOAverage'}"
    distribution_average_phenotypes_by_patients_cluster.each do |patient_num, ave|
      f.puts "#{patient_num}\t#{ave}"
    end
  end
  #return distribution_average_phenotypes_by_patients_cluster
end