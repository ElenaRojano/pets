require "statistics2"
require "terminal-table"
#require "report_html"
#require 'bigdecimal'

def search4HPO(info2predict, trainingData)
	#search if there are profile HPOs within the association file 
	hpo_regions = {}
	info2predict.each do |hpo|
		regions = trainingData[hpo]
		if !regions.nil?
			hpo_regions[hpo] = regions
		end
	end
	return hpo_regions
end

def group_by_region(hpo_regions)
	#hpo_regions-> hpo => [[chr, start, stop, regID, score], [...]]
	region2hpo = {}
	regionAttributes = {}
	association_scores = {}
	hpo_regions.each do |hpo, regions|
		regions.each do |chr, start, stop, regionID, association_score|
			query = region2hpo[regionID]
			if query.nil?
				region2hpo[regionID] = [hpo]
			else
				query << hpo
			end
			query = regionAttributes[regionID]
			if query.nil?
				total_patients_in_region = regionID.split('.')[3].to_i
				region_length = stop - start
				regionAttributes[regionID] = [chr, start, stop, total_patients_in_region, region_length]
			end
			query = association_scores[regionID]
			if query.nil?
				association_scores[regionID] = {hpo => association_score}
			else
				query[hpo] = association_score 
			end
		end		
	end
	return region2hpo, regionAttributes, association_scores
end

def add_parentals_of_not_found_hpos_in_regions(
	patient_hpo_profile, 
	trainingData, 
	region2hpo, 
	regionAttributes, 
	association_scores,
	hpo_metadata # hpo_code => [phenotype, relations]
	)
	new_hpos = []
	region2hpo.each do |regionID, hpos|
		hpos_not_found = patient_hpo_profile - hpos
		parental_hpos = []
		hpo_scores = {}
		hpos_not_found.each do |hpo|
			region, parental_hpo = get_region_with_parental_hpo(hpo, regionID, trainingData , hpo_metadata)
			if !region.nil? && 
				!parental_hpos.include?(parental_hpo) && 
				!patient_hpo_profile.include?(parental_hpo)
				parental_hpos << parental_hpo
				hpo_scores[parental_hpo] = region.last
			end
		end
		hpos.concat(parental_hpos)
		new_hpos.concat(parental_hpos)
		association_scores[regionID].merge!(hpo_scores)
	end
	patient_hpo_profile.concat(new_hpos.uniq)
end

def get_region_with_parental_hpo(hpo, regionID, trainingData , hpo_metadata)
	region = nil
	final_hpo = nil
	hpos = [hpo]
	while !hpos.empty?
		temp = []
		hpos.each do |hp| 
			hpo_data = hpo_metadata[hp]
			if !hpo_data.nil?
				main_hpo_code, phenotype, relations = hpo_data
				temp.concat(relations.map{|rel| rel.first})
			end
		end
		temp.each do |temp_hpo|
			regions = trainingData[temp_hpo]
			if !regions.nil?
				final_reg = regions.select{|reg| reg[3] == regionID}
				if !final_reg.empty?
					region = final_reg.first
					final_hpo = temp_hpo
					temp = []
					break
				end
			end
		end
		hpos = temp
	end
	return region, final_hpo
end

def generate_hpo_region_matrix(region2hpo, association_scores, info2predict, null_value=0)
	# #method for creating the hpo to region matrix for plotting
	# #info2predict = hpo list from user
	# #hpo_associated_regions = [[chr, start, stop, [hpos_list], [weighted_association_scores]]]
	hpo_region_matrix = []
	region2hpo.each do |regionID, hpos_list|
		row = []
		info2predict.each do |user_hpo|
			value = association_scores[regionID][user_hpo]
			if value.nil?
				row << null_value
			else
				row << value
			end
		end
		hpo_region_matrix << row
	end
	return hpo_region_matrix
end

def scoring_regions(regionAttributes, hpo_region_matrix, scoring_system, pvalue_cutoff, freedom_degree, null_value=0)
	#hpo_associated_regions = [[chr, start, stop, [hpos_list], [weighted_association_scores]]]
	#hpo_region_matrix = [[0, 0.4, 0, 0.4], [0, 0, 0.5, 0.4]]
	# STDERR.puts "EH"
	regionAttributes_array = regionAttributes.values
	max_cluster_length = hpo_region_matrix.map{|x| x.count {|i| i != 0}}.max if freedom_degree == 'maxnum'
	hpo_region_matrix.each_with_index do |associations, i|
		sample_length = nil
		if freedom_degree == 'prednum'
			sample_length = associations.length
		elsif freedom_degree == 'phennum'
			sample_length = associations.count{|s| s != 0}
		elsif freedom_degree == 'maxnum'
			sample_length = max_cluster_length
		else
			abort("Invalid freedom degree calculation method: #{freedom_degree}")
		end

		if scoring_system == 'mean'
			mean_association = associations.inject(0){|s,x| s + x } / sample_length
			regionAttributes_array[i] << mean_association
		elsif scoring_system == 'fisher'
			#hyper must be ln not log10 from net analyzer
			#https://en.wikipedia.org/wiki/Fisher%27s_method
			# STDERR.puts associations.inspect
			lns = associations.map{|a| Math.log(10 ** -a)} #hyper values come as log10 values
			sum = lns.inject(0){|s, a| s + a} 
			combined_pvalue = Statistics2.chi2_x(sample_length *2, -2*sum)	
			regionAttributes_array[i] << combined_pvalue
		elsif scoring_system == 'stouffer'
			sum = associations.inject(0){|s,x| s + x}
			combined_z_score = sum/Math.sqrt(sample_length)
			regionAttributes_array[i] << combined_z_score
		elsif scoring_system == 'geommean'
			#NOTE: if troubles with number size, use BigDecimal
			geommean_mult = associations.inject(1){|s,x| s * x}
			geommean_association = geommean_mult.to_f ** ( sample_length ** -1 )
			regionAttributes_array[i] << geommean_association
		elsif scoring_system == 'median'
			median_value = associations.length / 2
			if median_value % 2 == 0
				median_up = associations.sort[median_value]
				median_down = associations.sort[median_value - 1]
				pair_median = ( median_up + median_down ) / 2
				median_association = associations.sort[pair_median.ceil]
			else
				median_association = associations.sort[median_value.ceil]
			end
			regionAttributes_array[i] << median_association
		elsif scoring_system == 'maxnum'
			max_association = associations.max
			regionAttributes_array[i] << max_association
		elsif scoring_system == 'minnum'
			min_association = associations.min
			regionAttributes_array[i] << min_association
		else 
			abort("Invalid ranking method: #{scoring_system}")
		end
	end
	if scoring_system == 'mean' || 
		scoring_system == 'geommean' ||
		scoring_system == 'maxnum' ||
		scoring_system == 'minnum'
		regionAttributes.select!{|regionID, attributes| attributes.last >= pvalue_cutoff}
	elsif scoring_system == 'fisher'
		regionAttributes.select!{|regionID, attributes| attributes.last <= pvalue_cutoff}
	end
	#Combined p-value: less value equals better association -> not due randomly.
end

def join_regions(regions)
	#[chr, start, stop, association_values.keys, association_values.values, score]
	merged_regions = []
	sorted_regions = regions.sort_by{|reg | [reg[0], reg[1]]}
	ref_reg = sorted_regions.shift
	while !sorted_regions.empty?
		next_region = sorted_regions.shift
		if ref_reg[0] == next_region[0] &&
			(ref_reg[2] - next_region[1]).abs <= 1 &&
			(ref_reg[5] - next_region[5]).abs.fdiv([ref_reg[5], next_region[5]].max) <= 0.05 &&
			ref_reg[3] == next_region[3]

			ref_reg[2] = next_region[2]
			ref_assoc_values = ref_reg[4] 
			next_assoc_values = next_region[4]
			assoc_values = []
			ref_assoc_values.each_with_index do |ref_val, i|
				#assoc_values << (ref_val + next_assoc_values[i]).fdiv(2)
				assoc_values << [ref_val, next_assoc_values[i]].max
			end
			ref_reg[4] = assoc_values
			#ref_reg[5] = (ref_reg[5] + next_region[5]).fdiv(2)
			ref_reg[5] = [ref_reg[5], next_region[5]].max
		else
			merged_regions << ref_reg
			ref_reg = next_region
		end
	end
	merged_regions << ref_reg
	return merged_regions
end

# def hpo_quality_control(prediction_data, hpo_metadata_file, information_coefficient_file)
def hpo_quality_control(prediction_data, hpos_ci_values, hpo)
	characterised_hpos = []
	##information_coef_file= hpo_code, ci
	##prediction_data = [hpo1, hpo2, hpo3...]
	# hpos_ci_values = load_hpo_ci_values(information_coefficient_file)
	prediction_data.each do |hpo_code|
		# names, rejected = hpo.translate_codes2names([hpo_code])
		names, rejected = hpo.translate_ids([hpo_code])
		tmp = [names.first, hpo_code] # col hpo name, col hpo code
		ci = hpos_ci_values[hpo_code]
		unless ci.nil? # col exists? and ci values
			tmp.concat(["yes", ci])
		else
			tmp.concat(["no", "-"])
		end
		# parent = prediction_data & hpo.get_parents(hpo_code)
		parent = prediction_data & hpo.get_ancestors(hpo_code, true)
		if parent.empty?
			parent << "-"
		else
			# n, r = hpo.translate_codes2names(parent)
			n, r = hpo.translate_ids(parent)
			parent = parent.zip(n) # Combine code ids with hpo names
		end
		tmp << parent # col parents
		specific_childs = hpo.get_childs_table([hpo_code], true)
		tmp << specific_childs.first.last
		characterised_hpos << tmp
	end
	return characterised_hpos
end

def calculate_hpo_recovery_and_filter(adjacent_regions_joined, patient_original_phenotypes, predicted_hpo_percentage, min_hpo_recovery_percentage, patient_number)     
  records_to_delete = []
  counter = 0
  adjacent_regions_joined.each do |chr, start, stop, hpo_list, association_values, score|
    hpo_coincidences = patient_original_phenotypes & hpo_list
    original_hpo_recovery_percentage = hpo_coincidences.length / patient_original_phenotypes.length.to_f * 100
    records_to_delete << counter if original_hpo_recovery_percentage < min_hpo_recovery_percentage
    query = predicted_hpo_percentage[patient_number] 
    if query.nil?
     predicted_hpo_percentage[patient_number] = [original_hpo_recovery_percentage]
    else
     query << original_hpo_recovery_percentage
    end   
    counter += 1 
  end
  records_to_delete.reverse_each do |record_number|
    adjacent_regions_joined.delete_at(record_number)
  end
end

def report_data(characterised_hpos, hpo_associated_regions, html_file, hpo, genes_with_kegg_data, pathway_stats)
	container = {:characterised_hpos => characterised_hpos,
	 	:merged_regions => hpo_associated_regions,
	 	:hpo => hpo,
	 	:genes_with_kegg_data => genes_with_kegg_data,
	 	:pathway_stats => pathway_stats
	 }
	template = File.open(File.join(REPORT_FOLDER, 'patient_report.erb')).read
	report = Report_html.new(container, 'Patient HPO profile summary')
	report.build(template)
	report.write(html_file)
end

def save_patient_matrix(output, patient_hpo_profile, regionAttributes, hpo_region_matrix)
	File.open(output, "w") do |f|
		f.puts "Region\t#{patient_hpo_profile.join("\t")}"
		regionAttributes_array = regionAttributes.values
		hpo_region_matrix.each_with_index do |association_values, i|
		  chr, start, stop = regionAttributes_array[i]
		  f.puts "#{chr}:#{start}-#{stop}\t#{association_values.join("\t")}"
		end
	end
end