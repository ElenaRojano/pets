#! /usr/bin/env ruby
# E. Rojano, 2017
# Code for generate table with CI and association values.
# Input files: dictionary HPO to CI and association values files (NetAnalyzer output)

##########################
#LIBRARIES
##########################
require 'optparse'

##########################
#METHODS
##########################

def load_dictionary(file)
	hpo_to_ci_dictionary = {}
	File.open(file).each do |line|
		line.chomp!
		patient, hpo, ci = line.split("\t")
		hpo_to_ci_dictionary[hpo] = ci
	end
	return hpo_to_ci_dictionary
end

def load_association_file_with_ci_values(file, hpo_to_ci_dictionary, association_cutoff)
	association_values_2_ci = []
	last_ci = nil
	File.open(file).each do |line|
		line.chomp!
		hpo, region, association_value = line.split("\t")
		association_value = association_value.to_f
		next if association_value <= association_cutoff
		current_ci = hpo_to_ci_dictionary[hpo]
		association_values_2_ci << [current_ci, association_value]
	end
	return association_values_2_ci
end

##############################
#OPTPARSE
##############################

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"

  options[:input_association] = nil
  opts.on("-a", "--input_association PATH", "Input association file") do |input_association|
    options[:input_association] = input_association
  end

  options[:input_ci] = nil
  opts.on("-c", "--input_ci PATH", "Input CI file") do |input_ci|
    options[:input_ci] = input_ci
  end

  options[:output_file] = "association2ci.txt"
  opts.on("-f", "--output_file PATH", "Output file with associations to ci") do |output_file|
    options[:output_file] = output_file
  end

  options[:association_cutoff] = 0
  opts.on("-o", "--association_cutoff FLOAT", "Association cutoff for filtering") do |association_cutoff|
    options[:association_cutoff] = association_cutoff.to_f
  end

  options[:add_header] = FALSE
  opts.on("-H", "--add_header", "Add headers") do
    options[:add_header] = TRUE
  end


end.parse!

##############################
#MAIN
##############################

hpo_to_ci_dictionary = load_dictionary(options[:input_ci])
association_values_2_ci = load_association_file_with_ci_values(options[:input_association], hpo_to_ci_dictionary, options[:association_cutoff])
handler = File.open(options[:output_file], 'w')
handler.puts "coef\tassoc" if options[:add_header]
association_values_2_ci.each do |association2ci|
	handler.puts association2ci.join("\t")
end
handler.close