#! /usr/bin/env ruby
#DECIPHER predictor system, using data from cross validation
#data2predict = file to predict
#training_file.txt = file with training data (association values and hpo codes).
require 'optparse'
require File.join(File.dirname(__FILE__), 'generalMethods.rb')

##########################
#METHODS
##########################


def predict_patient(prediction_file, training_set, thresold)
	File.open(prediction_file).each do |line|
		line.chomp!
		fields = line.split("\t")
		chr = fields.shift
		pt_start = fields.shift.to_i
		pt_stop = fields.shift.to_i
    #associationValue = fields[1]
		#next if associationValue < thresold
    query = training_set[chr]
		puts training_set.inspect
    if !query.nil?
			query.each do |hpo_start, hpo_stop, hpo_code, association_score|
				if (hpo_stop > pt_start && hpo_stop <= pt_stop) ||
					(hpo_start >= pt_start && hpo_start < pt_stop) ||
					(hpo_start <= pt_start && hpo_stop >= pt_stop) ||
					(hpo_start > pt_start && hpo_stop < pt_stop) 
					puts [chr, pt_start, pt_stop].concat(fields).concat([hpo_code, association_score, hpo_start, hpo_stop]).join("\t") if association_score >= thresold
				end
			end
		end
	end
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

  options[:prediction_file] = nil
  #chr\tstart\tstop
  opts.on("-p", "--prediction_file PATH", "Input file with regions to predict") do |input_path|
    options[:prediction_file] = input_path
  end

  options[:association_limit] = 0
  opts.on("-l", "--association_limit FLOAT", "Thresold for association values") do |association_limit|
  	options[:association_limit] = association_limit.to_f
  end
  # options[:output_path] = "patient_file_overlapping.txt"
  # opts.on("-o", '--output_path PATH', 'Output path for overlapping patient file') do |output_path|
  # 	options[:output_path] = output_path
  # end

end.parse!


##########################
#MAIN
##########################

training_set = load_training_file4regions(options[:training_file])
predict_patient(options[:prediction_file], training_set, options[:association_limit])

