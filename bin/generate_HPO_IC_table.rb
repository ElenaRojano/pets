#! /usr/bin/env ruby

ROOT_PATH = File.dirname(__FILE__)
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'pets'))

require 'optparse'
require 'generalMethods.rb'
require 'semtools'

##############################
#METHODS
##############################

def load_paco_files(path, hpo)
	patient_profiles = {}
  counter = 0
	Dir.glob(path).each do |input_file|
		File.open(input_file).each do |line|
			next if line.include?('patient_id')
			line.chomp!
      data = line.split("\t")
      phens = data.last.split('|').map{|a| a.to_sym}
      phens = hpo.clean_profile_hard(phens)
      patient_profiles[(data.first + counter.to_s).to_sym] = phens
      counter += 1
		end
	end
  return patient_profiles
end

##############################
#OPTPARSE
##############################

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"

  options[:hpo_file] = nil
  opts.on("-b", "--hpo_file PATH", "Input HPO obo file") do |hpo_file|
    options[:hpo_file] = hpo_file
  end

  options[:input_path] = nil
  opts.on("-i", "--input_path PATH", "Input path with files to combine") do |data|
    options[:input_path] = data
  end

  options[:output_file] = 'output.txt'
  opts.on("-o", "--output_file PATH", "Output file") do |data|
    options[:output_file] = data
  end

  opts.on_tail("-h", "--help", "Show this message") do
    puts opts
    exit
  end

end.parse!

##############################
#MAIN
##############################
hpo = Ontology.new(file: options[:hpo_file], load_file: true)

patient_profiles = load_paco_files(options[:input_path], hpo)
hpo.load_profiles(patient_profiles)
onto_ic, freq_ic = hpo.get_observed_ics_by_onto_and_freq
File.open(options[:output_file], 'w') do |f|
  freq_ic.each do |term, freq|
    f.puts "#{term.to_s}\t#{freq}"
  end
end