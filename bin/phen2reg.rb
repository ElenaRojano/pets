#! /usr/bin/env ruby
# E. Rojano, September 2016
# Program to predict the position from given HPO codes, sorted by their association values.
# TODO: return HPO located in the same region


#data2predict = file to predict
#training_file.txt = file with training data (association values and hpo codes).
require 'optparse'
require File.join(File.dirname(__FILE__), 'generalMethods.rb')

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


def merge_regions(region2hpo, regionAttributes, association_scores)
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
					add_region(merged_regions, tmp_chr, tmp_start, tmp_stop, hpos_list, regionIDs, association_scores, regions_lengths, patients_numbers)
					tmp_chr = nil
					tmp_start = nil
					tmp_stop = nil
				end
			else
				add_region(merged_regions, tmp_chr, tmp_start, tmp_stop, hpos_list, regionIDs, association_scores, regions_lengths, patients_numbers)
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
		add_region(merged_regions, tmp_chr, tmp_start, tmp_stop, hpos_list, regionIDs, association_scores, regions_lengths, patients_numbers)
	end	
	return merged_regions
end

def add_region(merged_regions, tmp_chr, tmp_start, tmp_stop, hpos_list, regionIDs, association_scores, region_lengths, patients_numbers)	
	#region_lengths = number of regions that have the same HPOs
	unless tmp_chr.nil? && tmp_start.nil? && tmp_stop.nil?
		association_values_by_region = regionIDs.map {|r| association_scores[r]}
		weighted_association_scores = []
		hpos_list.each do |hpo|
			scores = association_values_by_region.map{|hpo_scores| hpo_scores[hpo] }
			weighted_score = 0
			weight = 0
			scores.each_with_index do |s, i|
				weighted_score += s * region_lengths[i] * patients_numbers[i]
				weight += region_lengths[i] * patients_numbers[i]
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

def ranking_regions(merged_regions, hpo_region_matrix)
	#merged_regions = [[chr, start, stop, [hpos_list], [weighted_association_scores]]]
	#hpo_region_matrix = [[0, 0.4, 0, 0.4], [0, 0, 0.5, 0.4]]
	hpo_region_matrix.each_with_index do |associations, i|
		mean_association = associations.inject(0){|s,x| s + x } / associations.length
		merged_regions[i] << mean_association
	end
	merged_regions.sort!{|r1, r2| r2.last <=> r1.last}
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
  opts.on("-b", "--best_thresold INTEGER", "Association value thresold") do |best_thresold|
    options[:best_thresold] = best_thresold.to_f
  end

  options[:merge_regions] = FALSE
  	opts.on("-j", "--merge_regions", "Join coordinates by borders") do
  options[:merge_regions] = TRUE
  end

  options[:output_matrix] = 'output_matrix.txt'
  opts.on("-o", "--output_matrix PATH", "Output matrix file, with associations for each input HPO") do |output_matrix|
    options[:output_matrix] = output_matrix
  end

end.parse!


##########################
#MAIN
##########################
trainingData = load_training_file4HPO(options[:training_file], options[:best_thresold])
hpo_regions = search4HPO(options[:prediction_data], trainingData)
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
		adjacent_regions_joined = merge_regions(region2hpo, regionAttributes, association_scores)
		#adjacent_regions_joined.each do |chr, start, stop, hpo_list, association_scores|
		#	puts "#{chr}\t#{start}\t#{stop}\t#{hpo_list.join(',')}\t#{association_scores.join(',')}"
		#end
		hpo_region_matrix = generate_hpo_region_matrix(adjacent_regions_joined, options[:prediction_data])
		output_matrix = File.open(options[:output_matrix], "w")
		headers = options[:prediction_data]
		output_matrix.puts "Region\t#{headers.join("\t")}"
		hpo_region_matrix.each_with_index do |association_values, i|
			chr, start, stop, hpo_list, association_scores = adjacent_regions_joined[i]
			output_matrix.puts "#{chr}:#{start}-#{stop}\t#{association_values.join("\t")}"
		end
		output_matrix.close
		ranking_regions(adjacent_regions_joined, hpo_region_matrix)
		adjacent_regions_joined.each do |chr, start, stop, hpo_list, association_score, association_mean|
			puts "#{chr}\t#{start}\t#{stop}\t#{hpo_list.join(',')}\t#{association_score.join(',')}\t#{association_mean}"
		end
	end
end