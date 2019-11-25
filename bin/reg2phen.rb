#! /usr/bin/env ruby
#DECIPHER predictor system, using data from cross validation
#data2predict = file to predict
#training_file.txt = file with training data (association values and hpo codes).
require 'optparse'
require File.join(File.dirname(__FILE__), '..', 'lib', 'gephepred', 'generalMethods.rb')

##########################
#METHODS
##########################


# def predict_patient(prediction_file, training_set, threshold)
#   predicted_hpos = []
#   File.open(prediction_file).each do |line|
#     line.chomp!
#     chr, pt_start, pt_stop = line.split("\t")
#     pt_start = pt_start.to_i
#     pt_stop = pt_stop.to_i
#     query = training_set[chr]
#     association_value = query.last.last
#     next if association_value.to_f < threshold
#     if !query.nil?
#       query.each do |hpo_start, hpo_stop, nodeID, hpo_code, association_value|
#         if (hpo_stop > pt_start && hpo_stop <= pt_stop) ||
#           (hpo_start >= pt_start && hpo_start < pt_stop) ||
#           (hpo_start <= pt_start && hpo_stop >= pt_stop) ||
#           (hpo_start > pt_start && hpo_stop < pt_stop) 
#           if association_value >= threshold
#             predicted_hpos << [chr, pt_start, pt_stop].concat([hpo_code, association_value])
#           end
#         end
#       end
#     end
#   end
#   # STDERR.puts predicted_hpos.inspect
#   return predicted_hpos
# end

def predict_patient(prediction_file, training_set, threshold)
  File.open(prediction_file).each do |line|
    line.chomp!
    # fields = line.split("\t")
    chr, pt_start, pt_stop = line.split("\t")
    pt_start = pt_start.to_i
    pt_stop = pt_stop.to_i
    #associationValue = fields[1]
    # next if associationValue < threshold
    query = training_set[chr]
    # STDERR.puts query.inspect
    if !query.nil?
      query.each do |hpo_start, hpo_stop, nodeID, hpo_code, association_score|
        # Process.exit
        if (hpo_stop > pt_start && hpo_stop <= pt_stop) ||
          (hpo_start >= pt_start && hpo_start < pt_stop) ||
          (hpo_start <= pt_start && hpo_stop >= pt_stop) ||
          (hpo_start > pt_start && hpo_stop < pt_stop) 
          puts [chr, pt_start, pt_stop].concat([hpo_code, association_score, hpo_start, hpo_stop]).join("\t") if association_score.to_f >= threshold
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
  opts.on("-l", "--association_limit FLOAT", "Threshold for association values") do |association_limit|
  	options[:association_limit] = association_limit.to_f
  end

  options[:output_path] = "output.txt"
  opts.on("-o", '--output_path PATH', 'Output path for overlapping patient file') do |output_path|
  	options[:output_path] = output_path
  end

  opts.on_tail("-h", "--help", "Show this message") do
    puts opts
    exit
  end

end.parse!


##########################
#MAIN
##########################

training_set = load_training_file4regions(options[:training_file])
predict_patient(options[:prediction_file], training_set, options[:association_limit])
# predicted_hpos = predict_patient(options[:prediction_file], training_set, options[:association_limit])
# handler = File.open(options[:output_path], 'w')
# predicted_hpos.each do |predictions|
#   handler.puts predictions.join("\t") 
# end
# handler.close
