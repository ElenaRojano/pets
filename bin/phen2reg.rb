#! /usr/bin/env ruby
# E. Rojano, September 2016
# Program to predict the position from given HPO codes, sorted by their association values.

#data2predict = file to predict
#training_file.txt = file with training data (association values and hpo codes).
require 'optparse'
require 'statistics2'
require 'terminal-table'
require 'report_html'
require File.join(File.dirname(__FILE__), '..', 'lib', 'gephepred', 'generalMethods.rb')

##########################
#REPORTING
##########################
REPORT_FOLDER=File.expand_path(File.join(File.dirname(__FILE__), '..', 'templates'))

##########################
#METHODS
##########################

def search4HPO(info2predict, trainingData)
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
	#puts hpo_regions.inspect # hpo=> [[chr, start, stop, nodo, assoc], [...]]
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
				regionAttributes[regionID] = [chr, start, stop, total_patients_in_region, region_length.to_i, regionID]
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

def sorting_regions_by_shared_hpos(region2hpo)
	arr_region2hpo = []
	region2hpo.each do |region, hpos|
		arr_region2hpo << [region, hpos.sort]
	end
	arr_region2hpo.sort!{|r1, r2| r1.last <=> r2.last}
	# # arr_region2hpo = [[1.1.A.1, [hpo1, hpo2, hpo3]], [1.2.A.1, [hpo1, hpo2, hpo3]]...]
	return arr_region2hpo
end

def cluster_regions_by_common_hpos(arr_region2hpo)
	regions_by_hpos = {}
	last_hpos = []
	regions = []
	arr_region2hpo.each do |region, hpos|
		if last_hpos == hpos
			regions << region
		else
			if !last_hpos.empty?
				regions_by_hpos[last_hpos] = regions
				regions = []
			end
		end
		last_hpos = hpos 
	end
	regions_by_hpos[last_hpos] = regions
	#puts regions_by_hpos.inspect
	# #regions_by_hpos = {[hpo1, hpo2, hpo3] => [1.1.A.1, 1.2.A.4, 1.3.A.12]...}
	return regions_by_hpos
end


def merge_regions(region2hpo, regionAttributes, association_scores, weight_style)
	# region2hpo = {region => [hpo1, hpo2...]}
	# regionAttributes = {region => [chr, start, stop, patients_number, region_length, region]}
	merged_regions = []
	arr_region2hpo = sorting_regions_by_shared_hpos(region2hpo)
	regions_by_hpos = cluster_regions_by_common_hpos(arr_region2hpo)
	regions_by_hpos.each do |hpos_list, regions|
		regionIDs = []
		regions_lengths = []
		patients_numbers = []
		region_attributes = regions.map { |region| regionAttributes[region] }
		region_attributes.sort! { |r1, r2| [r1[0], r1[1]] <=> [r2[0], r2[1]] }
		tmp_chr = nil
		tmp_start = nil
		tmp_stop = nil
		region_attributes.each_with_index do |attributes, counter|
			break if counter + 1 == region_attributes.length
			cur_chr, cur_start, cur_stop, cur_patients_number, cur_region_length, cur_regionID = attributes
			next_chr, next_start, next_stop, next_patients_number, next_region_length, next_regionID = region_attributes[counter + 1]
			if cur_chr == next_chr
				if cur_stop == next_start || cur_stop == next_start + 1
					tmp_chr = cur_chr
					tmp_start = cur_start if tmp_start.nil?
					tmp_stop = cur_stop
				else
					add_region(merged_regions, tmp_chr, tmp_start, tmp_stop, hpos_list, regionIDs, association_scores, regions_lengths, patients_numbers, weight_style)
					tmp_chr = nil
					tmp_start = nil
					tmp_stop = nil
				end
			else
				add_region(merged_regions, tmp_chr, tmp_start, tmp_stop, hpos_list, regionIDs, association_scores, regions_lengths, patients_numbers, weight_style)
				tmp_chr = nil
				tmp_start = nil
				tmp_stop = nil
			end
			regionIDs << cur_regionID if regionIDs.empty?
			regionIDs << next_regionID
			regions_lengths << cur_region_length if regions_lengths.empty?
			regions_lengths << next_region_length
			patients_numbers << cur_patients_number if patients_numbers.empty?
			patients_numbers << next_patients_number
		end
		add_region(merged_regions, tmp_chr, tmp_start, tmp_stop, hpos_list, regionIDs, association_scores, regions_lengths, patients_numbers, weight_style)
	end	
	return merged_regions
end

def add_region(merged_regions, tmp_chr, tmp_start, tmp_stop, hpos_list, regionIDs, association_scores, region_lengths, patients_numbers, weight_style)	
	#region_lengths = number of regions that have the same HPOs
	unless tmp_chr.nil? && tmp_start.nil? && tmp_stop.nil?
		association_values_by_region = regionIDs.map {|r| association_scores[r]}
		weighted_association_scores = []
		hpos_list.each do |hpo|
			scores = association_values_by_region.map{|hpo_scores| hpo_scores[hpo] }
			weighted_score = 0
			weight = 0
			scores.each_with_index do |s, i|
				if weight_style == 'double'
					weighted_score += s * region_lengths[i] * patients_numbers[i]
					weight += region_lengths[i] * patients_numbers[i]
				elsif weight_style == 'simple'
					weighted_score += s * region_lengths[i]
					weight += region_lengths[i]
				else
					abort("Invalid weight method: #{weight_style}")
				end
			end
			weighted_association_scores << weighted_score/weight			
		end
		merged_regions << [tmp_chr, tmp_start, tmp_stop, hpos_list, weighted_association_scores]
	end
end

def generate_hpo_region_matrix(merged_regions, info2predict)
	# #info2predict = hpo list from user
	# #merged_regions = [[chr, start, stop, [hpos_list], [weighted_association_scores]]]
	hpo_region_matrix = []
	merged_regions.each do |chr, start, stop, hpos_list, association_scores|
		row = []
		info2predict.each do |user_hpo|
			pos = hpos_list.index(user_hpo)
			if pos.nil?
				row << 0
			else
				row << association_scores[pos]
			end
		end
		hpo_region_matrix << row
	end
	return hpo_region_matrix
end

def ranking_regions(merged_regions, hpo_region_matrix, ranking_style, pvalue_cutoff)
	#merged_regions = [[chr, start, stop, [hpos_list], [weighted_association_scores]]]
	#hpo_region_matrix = [[0, 0.4, 0, 0.4], [0, 0, 0.5, 0.4]]
	hpo_region_matrix.each_with_index do |associations, i|
		if ranking_style == 'mean'
			mean_association = associations.inject(0){|s,x| s + x } / associations.length
			merged_regions[i] << mean_association
		elsif ranking_style == 'fisher'
			#https://en.wikipedia.org/wiki/Fisher%27s_method
			lns = associations.map{|a| Math.log(Math::E ** -a)}	
			sum = lns.inject(0){|s, a| s + a} 
			combined_pvalue = Statistics2.chi2_x(lns.length * 2, -2 * sum)
			merged_regions[i] << combined_pvalue
		else 
			abort("Invalid ranking method: #{ranking_style}")
		end
	end
	if ranking_style == 'mean'
		merged_regions.sort!{|r1, r2| r2.last <=> r1.last}
		merged_regions.select!{|a| a.last >= pvalue_cutoff}
	elsif ranking_style == 'fisher'
		merged_regions.sort!{|r1, r2| r1.last <=> r2.last}
		merged_regions.select!{|a| a.last <= pvalue_cutoff}
	end
	#Combined p-value: less value equals better association, not due randomly.
end

def hpo_quality_control(prediction_data, hpo_metadata_file, information_coefficient_file)
	characterised_hpos = []
	#information_coef_file= hpo_code, ci
	#prediction_data = [hpo1, hpo2, hpo3...]
	hpo_metadata = load_hpo_metadata(hpo_metadata_file)
	hpo_child_metadata = inverse_hpo_metadata(hpo_metadata)
	#hpo_metadata = {hpo_code => [phenotype, relations]}, relations = [hpo_code_relation, name_relation]
	hpos_ci_values = {}
	File.open(information_coefficient_file).each do |line|
		line.chomp!
		hpo_code, ci = line.split("\t")
		hpos_ci_values[hpo_code] = ci.to_f
	end
	prediction_data.each do |hpo_code|
		tmp = []
		ci = hpos_ci_values[hpo_code]
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

def report_data(characterised_hpos, merged_regions, html_file, hpo_metadata)
	container = {:characterised_hpos => characterised_hpos,
	 	:merged_regions => merged_regions,
	 	:hpo_metadata => hpo_metadata}
	template = File.open(File.join(REPORT_FOLDER, 'patient_report.erb')).read
	report = Report_html.new(container, 'Patient HPO profile summary')
	report.build(template)
	report.write(html_file)
end
##########################
#OPT-PARSER
##########################

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"

  options[:training_file] = nil
  #chr\tstart\tstop\tphenotype\tassociation_value
  opts.on("-t", "--training_file PATH", "Input training file, with association values") do |training_path|
    options[:training_file] = training_path
  end

  options[:HPO_file] = nil
  opts.on("-H", "--training_file PATH", "Input HPO file, used as dictionary") do |hpo_file|
    options[:HPO_file] = hpo2_file
  end

  options[:prediction_data] = []
  #chr\tstart\tstop
  opts.on("-p", "--prediction_file PATH", "Input data with HPO codes for predicting their location. It can be either, a file path or string with HPO separated by commas") do |input_path|
    if File.exist?(input_path)
    	options[:prediction_data] = File.open(input_path).readlines.map!{|line| line.chomp}
    else
    	options[:prediction_data] = input_path.split(',')
    end
  end

  options[:group_by_region] = FALSE
  opts.on("-s", "--group_by_region", "Predict which HPOs are located in the same region") do
  	options[:group_by_region] = TRUE
  end

  options[:best_thresold] = 1.5
  opts.on("-b", "--best_thresold FLOAT", "Association value thresold") do |best_thresold|
    options[:best_thresold] = best_thresold.to_f
  end

  options[:pvalue_cutoff] = 0.1
  opts.on("-P", "--pvalue_cutoff FLOAT", "P-value cutoff") do |pvalue_cutoff|
    options[:pvalue_cutoff] = pvalue_cutoff.to_f
  end

  options[:merge_regions] = FALSE
  	opts.on("-j", "--merge_regions", "Join coordinates by borders") do
  options[:merge_regions] = TRUE
  end

  options[:output_matrix] = 'output_matrix.txt'
  opts.on("-o", "--output_matrix PATH", "Output matrix file, with associations for each input HPO") do |output_matrix|
    options[:output_matrix] = output_matrix
  end

  options[:weight_style] = ''
  opts.on("-w", "--weight_style STRING", "Weight style: simple (regions length) or double weight (regions length and patients number") do |weight_style|
    options[:weight_style] = weight_style
  end

  options[:ranking_style] = ''
  opts.on("-r", "--ranking_style STRING", "Ranking style: mean or fisher") do |ranking_style|
    options[:ranking_style] = ranking_style
  end

  options[:hpo_is_name] = TRUE
  	opts.on("-n", "--hpo_is_name", "Set this flag if phenotypes are given as names instead of codes") do
  options[:hpo_is_name] = FALSE
  end  

  options[:hpo2name_file] = nil
  opts.on("-f", "--hpo2name_file PATH", "Input hpo2name file") do |hpo2name_file|
    options[:hpo2name_file] = hpo2name_file
  end

  options[:information_coefficient] = nil
  opts.on("-i", "--hpo2name_file PATH", "Input file with information coefficients") do |information_coefficient|
    options[:information_coefficient] = information_coefficient
  end

  options[:output_quality_control] = "output_quality_control.txt"
  opts.on("-O", "--output_quality_control PATH", "Output file with quality control of all input HPOs") do |output_quality_control|
    options[:output_quality_control] = output_quality_control
  end

  options[:html_file] = "patient_profile_report.html"
  opts.on("-F", "--html_file PATH", "HTML file with patient information HPO profile summary") do |html_file|
    options[:html_file] = html_file
  end

end.parse!


##########################
#MAIN
##########################
if options[:hpo_is_name]
	hpo_dictionary = load_hpo_dictionary_name2code(options[:HPO_file])
	options[:prediction_data].each_with_index do |name, i|
		hpo_code = hpo_dictionary[name]
		if hpo_code.nil?
			abort("Invalid HPO name: #{name}.")
		end
		options[:prediction_data][i] = hpo_code
	end
end

#HPO quality control
#---------------------------
characterised_hpos, hpo_metadata = hpo_quality_control(options[:prediction_data], options[:hpo2name_file], options[:information_coefficient])
output_quality_control = File.open(options[:output_quality_control], "w")
header = ["HPO name", "HPO code", "Exists?", "CI value", "Is child of", "Childs"]
output_quality_control.puts Terminal::Table.new :headings => header, :rows => characterised_hpos
output_quality_control.close

#Prediction steps
#---------------------------
trainingData = load_training_file4HPO(options[:training_file], options[:best_thresold])
hpo_regions = search4HPO(options[:prediction_data], trainingData)
adjacent_regions_joined = []
if hpo_regions.empty?
	puts "Results not found"
elsif options[:group_by_region] == FALSE
	hpo_regions.each do |hpo, regions|
		regions.each do |region|
			puts "#{hpo}\t#{region.join("\t")}"
		end
	end
elsif options[:group_by_region] == TRUE
	region2hpo, regionAttributes, association_scores = group_by_region(hpo_regions)
	if options[:merge_regions] == FALSE
		# region2regions = {hpo => [[chr, start, stop, assocval],[...]]}
		region2hpo.each do |region, hpoCodes|
			puts "#{region}\t#{hpoCodes.join(",")}"
		end
	elsif options[:merge_regions] == TRUE
		adjacent_regions_joined = merge_regions(region2hpo, regionAttributes, association_scores, options[:weight_style])
		hpo_region_matrix = generate_hpo_region_matrix(adjacent_regions_joined, options[:prediction_data])
		output_matrix = File.open(options[:output_matrix], "w")
		headers = options[:prediction_data]
		output_matrix.puts "Region\t#{headers.join("\t")}"
		hpo_region_matrix.each_with_index do |association_values, i|
			chr, start, stop, hpo_list, association_scores = adjacent_regions_joined[i]
			output_matrix.puts "#{chr}:#{start}-#{stop}\t#{association_values.join("\t")}"
		end
		output_matrix.close
		ranking_regions(adjacent_regions_joined, hpo_region_matrix, options[:ranking_style], options[:pvalue_cutoff])
		adjacent_regions_joined.each do |chr, start, stop, hpo_list, association_score, association_mean|
			puts "#{chr}\t#{start}\t#{stop}\t#{hpo_list.join(',')}\t#{association_score.join(',')}\t#{association_mean}"
		end
	end
end

#Creating html report
#-------------------
report_data(characterised_hpos, adjacent_regions_joined, options[:html_file], hpo_metadata)