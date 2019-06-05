#! /usr/bin/env ruby
#Tool for calculating averages between different association values file.
#File structure: prec rec cut meth
#Load all files (7) stored in the same directory and calculate average;
#of lines for each method. Return a file with the same structure;
#giving name as "average" to the last column

require 'optparse'

##########################
#METHODS
##########################

def load_association_file(filename)
	fileInfo = []
	header = ''
	line_number = 0
	File.open(filename).each do |line|
		line.chomp!
		if line_number == 0
			header = line
		else
			cut, precision, recall, meth = line.split("\t")
			fileInfo << [cut.to_f, precision.to_f, recall.to_f, meth]
		end
		line_number += 1
	end
	return fileInfo, header
end

def calculate_average(all_files, cols_for_average)
	average = []
	n_files = all_files.length.to_f
	ref_file = all_files.shift
	summatory_file = []
	ref_file.each_with_index do |line, i|
		all_files.each do |file|
			line2 = file[i]
			cols_for_average.each do |col|
				line[col] = line[col] + line2[col]
			end
		end
		summatory_file << line
	end
	summatory_file.each do |line|
		cols_for_average.each do |col|
			line[col] = line[col]/n_files
		end
		average << line
	end
	return average
end


##########################
#OPT-PARSER
##########################

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"
 
  options[:file_names] = nil
  opts.on("-f", "--file_names STRING", "Input file names to calculate averages. Please separate names by commas") do |file_names|
    options[:file_names] = file_names.split(',')
  end

  options[:which_cols] = nil
  opts.on("-c", "--which_cols STRING", "Cols for performing average analysis") do |which_cols|
    options[:which_cols] = which_cols.split(',').map{|i| i.to_i - 1}
  end

end.parse!

##########################
#MAIN
##########################

all_files = []
header = nil
options[:file_names].each do |filename|
	file, header = load_association_file(filename)
	all_files << file
end

average = calculate_average(all_files, options[:which_cols])

puts header
average.each do |line|
	puts line.join("\t")
end
	


