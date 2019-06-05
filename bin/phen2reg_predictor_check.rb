#! /usr/bin/env ruby

##########################
#LIBRARIES
##########################
require 'optparse'

##########################
#METHODS
##########################
@success_percentage_distribution = []
@prediction_vector = []
@rankings = { :in => [], :out => []}
@genome_fraction_predicted = 0 #All positive cases
@good_predicted_subregions = 0 #True positive cases

def compute_rankings
	if !@prediction_vector.empty?
		n_preds = @prediction_vector.length.fdiv(100)
		@prediction_vector.each_with_index do |in_control, i|
			ranking = (i + 1).fdiv(n_preds)
			if in_control
				@rankings[:in] << ranking
			else
				@rankings[:out] << ranking
			end
		end
		@prediction_vector = []
	end
end

def load_prediction_file(input_file)
	predicted_regions = []
	File.open(input_file).each do |line|
		line.chomp!
		predicted_info = line.split("\t")
		profile_index = predicted_info[0].gsub('ProfID:','').to_i
		if predicted_info[1] != 'Results not found'
			predicted_hpos_number = predicted_info[4].split(',').length
			predicted_regions << [profile_index, predicted_info[1], predicted_info[2].to_i, predicted_info[3].to_i, predicted_info[6].to_f, predicted_hpos_number] 
		else
			predicted_regions << [profile_index]
		end
	end
	return predicted_regions # profile_index, pred_chr, pred_start, pred_stop, score
end

def get_imputation_scores(predicted_regions, integration_method)
	# predicted_regions.map{|pred_reg| pred_reg.pop}
	# STDERR.puts predicted_regions.inspect
	selected_regions = predicted_regions.select{|r| r[4].class == Float }
	score_regionLength_pairs = selected_regions.map{|r| [r[4], r[3] - r[2]]} #Get combined score and region length
	score_regionLength_pairs.sort!{|p1, p2| p1.first <=> p2.first}
	score_regionLength_pairs.reverse! if integration_method == 'fisher'
	total_region_length = score_regionLength_pairs.map{|p| p.last }.inject{|sum, l| sum + l }
	
	length2inspect = total_region_length/1000		 
	acumulated_score = 0
	inspected_length = 0
	score_regionLength_pairs.each do |score, length|
		acumulated_score += score * length
		inspected_length += length
		break if inspected_length >= length2inspect
	end
	return acumulated_score.fdiv(inspected_length)
end

def generate_random_imp_score(imputation_score, desv)
	range = imputation_score * desv * rand()
	if [true, false].sample
		final_score = imputation_score - range
	else
		final_score = imputation_score + range
	end
	return final_score
end

def load_patient_data(input_data_file)
	patient_data = []
	File.open(input_data_file).each do |line|
		line.chomp!
		mutation_coords, hpo_profile = line.split("\t")
		number_of_phenotypes = hpo_profile.split('|').length
		chr, start_pos, stop_pos = mutation_coords.split(':')
		patient_data << [chr, start_pos.to_i, stop_pos.to_i, number_of_phenotypes]	
	end
	return patient_data #ctrl_chr, ctrl_start, ctrl_stop, #number_of_phens
end

def get_perfomance_table(ctrl_regions, predicted_regions, scale, imputation_score, hpo_min_recovery, apply_imputation)
	table = []
	last_profile_id = ctrl_chr = ctrl_start = ctrl_stop = predicted_hpos_number = number_of_phenotypes = nil
	in_out_regions = []
	predicted_regions.each do |profile_index, pred_chr, pred_start, pred_stop, score, predicted_hpos_number|
		if last_profile_id != profile_index && !last_profile_id.nil?
			table.concat(process_in_out_regions(ctrl_start, ctrl_stop, scale, in_out_regions, imputation_score, apply_imputation))
			@genome_fraction_predicted = 0 
			@good_predicted_subregions = 0
			in_out_regions = []
			compute_rankings
		end
		ctrl_chr, ctrl_start, ctrl_stop, number_of_phenotypes = ctrl_regions[profile_index] #get position in array, for each prediction
		unless predicted_hpos_number.nil? || number_of_phenotypes.nil?
			hpo_recovery_percentage = ( predicted_hpos_number / number_of_phenotypes.to_f ) * 100 
			#STDERR.puts "#{predicted_hpos_number}\t#{number_of_phenotypes}"
			if hpo_recovery_percentage > hpo_min_recovery
				in_out_regions.concat(get_in_out_regions(ctrl_chr, ctrl_start, ctrl_stop, pred_chr, pred_start, pred_stop, score))
				last_profile_id = profile_index	
			end
		end
	end
	table.concat(process_in_out_regions(ctrl_start, ctrl_stop, scale, in_out_regions, imputation_score, apply_imputation))
	return table
end

def process_in_out_regions(ctrl_start, ctrl_stop, scale, in_out_regions, imputation_score, apply_imputation)
	@success_percentage_distribution << get_sucess_percentage(in_out_regions)
	table = []
	ctrl_length = ctrl_stop - ctrl_start
	non_predicted_regions = ctrl_length - @good_predicted_subregions
	if non_predicted_regions > 0
		#total_predicted_region_length = in_out_regions.map{|s| s.last}.inject(0){|i, sum| i + sum}
		#imputation_score = in_out_regions.map{|s| s[1] * s.last}.inject(0){|i, sum| i + sum}.fdiv(total_predicted_region_length)
		#imputation_score += 0.25 * imputation_score
		#imputation_score = in_out_regions.map{|s| s[1]}.max
		#index = (9 * in_out_regions.length).fdiv(10).ceil
		#imputation_score = in_out_regions.map{|s| s[1]}.sort[index]


		#in_out_regions << ["in", generate_random_imp_score(0.764, 0.35), non_predicted_regions] if apply_imputation
		in_out_regions << ["in", generate_random_imp_score(imputation_score, 0.35), non_predicted_regions] if apply_imputation
	end
	evaluated_genome_fraction = @genome_fraction_predicted + ctrl_length
	in_out_regions.each do |group, score, region_length|
		list_entries = (region_length.fdiv(evaluated_genome_fraction) * scale).ceil	
		table.concat(Array.new(list_entries, [group, score])) 
	end
	return table
end

def get_in_out_regions(ctrl_chr, ctrl_start, ctrl_stop, pred_chr, pred_start, pred_stop, score)
	in_out_regions = []
	if ctrl_chr == pred_chr
		if pred_start < ctrl_start && pred_stop > ctrl_stop # predicted region larger than ctrl region	
			region_length = ctrl_stop - ctrl_start
			in_out_regions << ["in", score, region_length]
			@good_predicted_subregions += region_length
			region_length = ctrl_start - pred_start
			in_out_regions << ["out", score, region_length]
			@genome_fraction_predicted += region_length
			region_length = pred_stop - ctrl_stop
			in_out_regions << ["out", score, region_length]
			@genome_fraction_predicted += region_length
		elsif pred_start >= ctrl_start && pred_stop <= ctrl_stop #within ctrl region		
			region_length = pred_stop - pred_start
			in_out_regions << ["in", score, region_length]
			@good_predicted_subregions += region_length
		elsif ctrl_start < pred_stop && ctrl_stop >= pred_stop #upstream region out of ctrl region
			region_length = pred_stop - ctrl_start
			in_out_regions << ["in", score, region_length]
			@good_predicted_subregions += region_length
			region_length = ctrl_start - pred_start
			in_out_regions << ["out", score, region_length]
			@genome_fraction_predicted += region_length
		elsif ctrl_start <= pred_start && ctrl_stop > pred_start #downstream region out of ctrl region
			region_length = ctrl_stop - pred_start
			in_out_regions << ["in", score, region_length]
			@good_predicted_subregions += region_length
			region_length = pred_stop - ctrl_stop
			in_out_regions << ["out", score, region_length]
			@genome_fraction_predicted += region_length
		else #in same chr but not in ctrl region
			region_length = pred_stop - pred_start
			in_out_regions << ["out", score, region_length]
			@genome_fraction_predicted += region_length
		end	
	elsif !pred_chr.nil? #in different chr
		region_length = pred_stop - pred_start
		in_out_regions << ["out", score, region_length]
		@genome_fraction_predicted += region_length
	end
	if in_out_regions.map{|reg| reg.first }.include?('in')
		@prediction_vector << true
	else
		@prediction_vector << false
	end
	return in_out_regions
end

def get_sucess_percentage(in_out_regions) 
	percentage = 0
	if !in_out_regions.empty?
		count = 0
		total = 0
		in_out_regions.each do |group, score, reg_length|
			count += reg_length if group == 'in'
			total += reg_length
		end
		percentage = count.fdiv(total)
	end
	return percentage
end


##########################
#OPT-PARSE
##########################

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"
  
  options[:input_prediction] = nil
  opts.on("-i", "--input_prediction PATH", "Input prediction file for checking") do |input_prediction|
    options[:input_prediction] = input_prediction
  end

  options[:meth] = nil
  opts.on("-m", "--meth STRING", "Method used in score integration calculation, affects to the imputation algorithm (if used)") do |meth|
    options[:meth] = meth
  end

  options[:output_file] = 'final_values_for_pr_curve.txt'
  opts.on("-o", "--output_file PATH", "Output results for PR curve") do |output_file|
    options[:output_file] = output_file
  end

  options[:hpo_recovery] = 0
  opts.on("-p", "--hpo_recovery INTEGER", "Minimum percentage of HPO terms to consider predictions") do |hpo_recovery|
    options[:hpo_recovery] = hpo_recovery.to_f
    abort("Please, choose a recovery value higher than 0") if options[:hpo_recovery] <= 0
  end

  options[:input_regions] = nil
  opts.on("-r", "--input_regions PATH", "Input patients true affected regions (ctrl file)") do |input_regions|
  	options[:input_regions] = input_regions
  end

  options[:success_percentage] = 'success_percentage'
  opts.on("-s", "--success_percentage PATH", "Output results with success percentage for each prediction") do |success_percentage|
    options[:success_percentage] = success_percentage
  end

  options[:apply_imputation] = false
  opts.on("-y", "--apply_imputation", "Activates imputation") do
    options[:apply_imputation] = true
  end

  options[:scale_size] = 100
  opts.on("-z", "--scale_size INTEGER", "Scale region size to avoid long range regions") do |scale_size|
    options[:scale_size] = scale_size.to_i
    abort("Please, choose a scale value higher than 0") if options[:scale_size] <= 0
  end


end.parse!

##########################
#MAIN
##########################

predicted_regions = load_prediction_file(options[:input_prediction])
imputation_score = get_imputation_scores(predicted_regions, options[:meth])
patient_data = load_patient_data(options[:input_regions])

table = get_perfomance_table(patient_data, predicted_regions, options[:scale_size], imputation_score, options[:hpo_recovery], options[:apply_imputation])
File.open(options[:output_file], 'w') do |f|
	f.puts "group\tscore"
	table.each do |output, score|
		if options[:apply_imputation]
			#score = generate_random_imp_score(0.764, 0.35) if score.nil?
			score = generate_random_imp_score(imputation_score, 0.35) if score.nil?
		else
			next if score.nil? #when no imputation
		end
		f.puts "#{output}\t#{score}"
	end
end
File.open(options[:success_percentage], 'w') do |f|
	f.puts 'perc'
	@success_percentage_distribution.each do |pg|
		f.puts pg
	end
end
File.open('ranking', 'w') do |f|
	f.puts "reg\tranking"
	@rankings.each do |reg, ranks|
		ranks.each do |rank|
			f.puts "#{reg.to_s}\t#{rank}"
		end
	end
end
	
