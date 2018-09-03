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


def generate_hpo_region_matrix(region2hpo, association_scores, info2predict)
	# #method for creating the hpo to region matrix for plotting
	# #info2predict = hpo list from user
	# #hpo_associated_regions = [[chr, start, stop, [hpos_list], [weighted_association_scores]]]
	hpo_region_matrix = []
	region2hpo.each do |regionID, hpos_list|
		row = []
		info2predict.each do |user_hpo|
			value = association_scores[regionID][user_hpo]
			if value.nil?
				row << 0
			else
				row << value
			end
		end
		hpo_region_matrix << row
	end
	return hpo_region_matrix
end

def scoring_regions(regionAttributes, hpo_region_matrix, scoring_system, pvalue_cutoff, freedom_degree)
	#hpo_associated_regions = [[chr, start, stop, [hpos_list], [weighted_association_scores]]]
	#hpo_region_matrix = [[0, 0.4, 0, 0.4], [0, 0, 0.5, 0.4]]
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
			lns = associations.map{|a| Math.log(10 ** -a)} #hyper values come as log10 values
			#https://en.wikipedia.org/wiki/Fisher%27s_method
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
			abort("Invalid ranking method: #{ranking_style}")
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

# def hpo_quality_control(prediction_data, hpo_metadata_file, information_coefficient_file)
def hpo_quality_control(prediction_data, hpo_metadata, hpo_child_metadata, hpos_ci_values)
	characterised_hpos = []
	##information_coef_file= hpo_code, ci
	##prediction_data = [hpo1, hpo2, hpo3...]
	##hpo_metadata = {hpo_code => [phenotype, relations]}, relations = [hpo_code_relation, name_relation]
	# hpo_metadata = load_hpo_metadata(hpo_metadata_file)
	# hpo_child_metadata = inverse_hpo_metadata(hpo_metadata)
	# hpos_ci_values = load_hpo_ci_values(information_coefficient_file)
	prediction_data.each do |hpo_code|
		tmp = []
		ci = hpos_ci_values[hpo_code]
		#STDERR.puts hpo_metadata[hpo_code]
		hpo_name, relations = hpo_metadata[hpo_code]
		tmp << hpo_name # col hpo name
		tmp << hpo_code # col hpo code
		unless ci.nil? # col exists? and ci values
			tmp << "yes"
			tmp << ci
		else
			tmp << "no"
			tmp << "-"
		end
		parent = check_parents(relations, prediction_data, hpo_metadata)		
		parent << "-" if parent.empty?
		tmp << parent # col parents
		childs = hpo_child_metadata[hpo_code]
		if childs.nil?
			childs = []
		else
			childs = childs.last
		end
		tmp << childs
		characterised_hpos << tmp
	end
	return characterised_hpos, hpo_metadata
end

def check_parents(relations, prediction_data, hpo_metadata)
	parent = []
	relations.each do |par_hpo_code, par_hpo_name|
		if prediction_data.include?(par_hpo_code)
			parent << [par_hpo_code, par_hpo_name]
		end
		grand_par_hpo = hpo_metadata[par_hpo_code]
		if !grand_par_hpo.nil?
			parent.concat(check_parents(grand_par_hpo.last, prediction_data, hpo_metadata))
		end
	end
	return parent
end

def report_data(characterised_hpos, hpo_associated_regions, html_file, hpo_metadata)
	container = {:characterised_hpos => characterised_hpos,
	 	:merged_regions => hpo_associated_regions,
	 	:hpo_metadata => hpo_metadata}
	template = File.open(File.join(REPORT_FOLDER, 'patient_report.erb')).read
	report = Report_html.new(container, 'Patient HPO profile summary')
	report.build(template)
	report.write(html_file)
end

##############################################################################
##############################################################################
##### OLD CODE FOR JOIN REGIONS BY BORDERS
##############################################################################
##############################################################################

# def sorting_regions_by_shared_hpos(region2hpo)
# 	#if regions share the same hpos, sort regions from lowest to highest
# 	#this method returns an array for its use in cluster_regions_by_common_hpos method
# 	arr_region2hpo = []
# 	region2hpo.each do |region, hpos|
# 		arr_region2hpo << [region, hpos.sort]
# 	end
# 	arr_region2hpo.sort!{|r1, r2| r1.last <=> r2.last}
# 	# # arr_region2hpo = [[1.1.A.1, [hpo1, hpo2, hpo3]], [1.2.A.1, [hpo1, hpo2, hpo3]]...]
# 	return arr_region2hpo
# end

# def cluster_regions_by_common_hpos(arr_region2hpo)
# 	#method for grouping hpos within different locations
# 	regions_by_hpos = {}
# 	last_hpos = []
# 	regions = []
# 	all_regions = []
# 	arr_region2hpo.each do |region, hpos|
# 		all_regions << region
# 		if last_hpos == hpos
# 			regions << region
# 		else
# 			regions_by_hpos[last_hpos] = regions if !last_hpos.empty?
# 			regions = [region]
# 		end
# 		last_hpos = hpos 
# 	end
# 	regions_by_hpos[last_hpos] = regions
# 	#puts regions_by_hpos.inspect
# 	# #regions_by_hpos = {[hpo1, hpo2, hpo3] => [1.1.A.1, 1.2.A.4, 1.3.A.12]...}
# 	return regions_by_hpos
# end

# def prepare_regions_for_profile_analysis(region2hpo, regionAttributes, association_scores, weight_style)
# 	# region2hpo = {region => [hpo1, hpo2...]}
# 	# regionAttributes = {region => [chr, start, stop, patients_number, region_length, region]}
# 	hpo_associated_regions = []
# 	arr_region2hpo = sorting_regions_by_shared_hpos(region2hpo)
# 	regions_by_hpos = cluster_regions_by_common_hpos(arr_region2hpo)
# 	regions_by_hpos.each do |hpos_list, regions|
# 		regionIDs = []
# 		regions_lengths = []
# 		patients_numbers = []
# 		region_attributes = regions.map { |region| regionAttributes[region] }
# 		region_attributes.each do |attributes|
# 			cur_chr, cur_start, cur_stop, cur_patients_number, cur_region_length, cur_regionID = attributes
# 			add_region(hpo_associated_regions, cur_chr, cur_start, cur_stop, hpos_list, [cur_regionID], association_scores, [cur_region_length], [cur_patients_number], weight_style)
# 		end
# 	end
# 	#puts hpo_associated_regions.inspect
# 	return hpo_associated_regions
# end

# def join_regions_by_borders(region2hpo, regionAttributes, association_scores, weight_style)
# 	# region2hpo = {region => [hpo1, hpo2...]}
# 	# regionAttributes = {region => [chr, start, stop, patients_number, region_length, region]}
# 	joined_regions_by_borders = []
# 	arr_region2hpo = sorting_regions_by_shared_hpos(region2hpo)
# 	regions_by_hpos = cluster_regions_by_common_hpos(arr_region2hpo)
# 	regions_by_hpos.each do |hpos_list, regions|
# 		regionIDs = []
# 		regions_lengths = []
# 		patients_numbers = []
# 		region_attributes = regions.map { |region| regionAttributes[region] }
# 		region_attributes.sort! { |r1, r2| [r1[0], r1[1]] <=> [r2[0], r2[1]] }
# 		tmp_chr = nil
# 		tmp_start = nil
# 		tmp_stop = nil
# 		region_attributes.each_with_index do |attributes, counter|
# 			break if counter + 1 == region_attributes.length
# 			cur_chr, cur_start, cur_stop, cur_patients_number, cur_region_length, cur_regionID = attributes
# 			next_chr, next_start, next_stop, next_patients_number, next_region_length, next_regionID = region_attributes[counter + 1]
# 			if cur_chr == next_chr
# 				if cur_stop == next_start || cur_stop == next_start + 1
# 					tmp_chr = cur_chr
# 					tmp_start = cur_start if tmp_start.nil?
# 					tmp_stop = cur_stop
# 				else
# 					add_region(joined_regions_by_borders, tmp_chr, tmp_start, tmp_stop, hpos_list, regionIDs, association_scores, regions_lengths, patients_numbers, weight_style)
# 					tmp_chr = nil
# 					tmp_start = nil
# 					tmp_stop = nil
# 				end
# 			else
# 				add_region(joined_regions_by_borders, tmp_chr, tmp_start, tmp_stop, hpos_list, regionIDs, association_scores, regions_lengths, patients_numbers, weight_style)
# 				tmp_chr = nil
# 				tmp_start = nil
# 				tmp_stop = nil
# 			end
# 			regionIDs << cur_regionID if regionIDs.empty?
# 			regionIDs << next_regionID
# 			regions_lengths << cur_region_length if regions_lengths.empty?
# 			regions_lengths << next_region_length
# 			patients_numbers << cur_patients_number if patients_numbers.empty?
# 			patients_numbers << next_patients_number
# 		end
# 		add_region(joined_regions_by_borders, tmp_chr, tmp_start, tmp_stop, hpos_list, regionIDs, association_scores, regions_lengths, patients_numbers, weight_style)
# 	end
# 	#puts joined_regions_by_borders.inspect	
# 	return joined_regions_by_borders
# end

# def add_region(hpo_associated_regions, tmp_chr, tmp_start, tmp_stop, hpos_list, regionIDs, association_scores, region_lengths, patients_numbers, weight_style)	
# 	#region_lengths = number of regions that have the same HPOs
# 	unless tmp_chr.nil? && tmp_start.nil? && tmp_stop.nil?
# 		association_values_by_region = regionIDs.map {|r| association_scores[r]}
# 		weighted_association_scores = []
# 		hpos_list.each do |hpo|
# 			scores = association_values_by_region.map{|hpo_scores| hpo_scores[hpo] }
# 			weighted_score = 0
# 			weight = 0
# 			if scores.length == 1
# 				weighted_score = scores.first
# 				weight = 1
# 			else
# 				scores.each_with_index do |s, i|
# 					if weight_style == 'double'
# 						weighted_score += s * region_lengths[i] * patients_numbers[i]
# 						weight += region_lengths[i] * patients_numbers[i]
# 					elsif weight_style == 'simple'
# 						weighted_score += s * region_lengths[i]
# 						weight += region_lengths[i]
# 					else
# 						abort("Invalid weight method: #{weight_style}")
# 					end
# 				end
# 			end
# 			weighted_association_scores << weighted_score/weight			
# 		end
# 		hpo_associated_regions << [tmp_chr, tmp_start, tmp_stop, hpos_list, weighted_association_scores]
# 	end
# end
