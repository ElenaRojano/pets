#! /usr/bin/env ruby


require 'optparse'

##########################
#METHODS
##########################

def load_file(input_file)
	nodeID_regions = {}
	File.open(input_file).each do |line|
		line.chomp!
		reg_start, reg_stop, reg_chr, nodeID = line.split("\t")
		node_properties = nodeID.split('.')
		patients_number = node_properties.pop.to_i
		query = nodeID_regions[reg_chr]
		metadata = [reg_start.to_i, reg_stop.to_i, patients_number]
		if query.nil?
			nodeID_regions[reg_chr] = [metadata]
		else
			query << metadata
		end
	end
	return nodeID_regions
end

def get_patients_per_bin(chromosomes, bin_size)
	bins = []
	chromosomes.each do |chr, regions|
		regions.sort!{|reg1, reg2| reg1[0] <=> reg2[0]}
		chr_length = regions.last[1]
		packages = (chr_length.to_f/bin_size).ceil

		package_array = Array.new(packages){|n| []}
		regions.each do |reg_start, reg_stop, patients_number|
			low_bin = (reg_start/bin_size).floor
			upper_bin = (reg_stop/bin_size).floor
			(low_bin..upper_bin).each do |bin_index|
				package_array[bin_index] << patients_number
			end
		end
		package_array.each_with_index do |patients_number_package, pkg_index|
			if !patients_number_package.empty?
				patients_average = patients_number_package.inject{ |s, i| s + i } / patients_number_package.length
			else
				patients_average = 0
			end	
			bins << [chr, pkg_index * bin_size, patients_average]
		end
	end
	return bins
end


##########################
#OPT-PARSER
##########################

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"
 
  options[:input_file] = nil
  opts.on("-f", "--input_file PATH", "Input file with regions and nodeIDs") do |input_file|
    options[:input_file] = input_file
  end

end.parse!

##########################
#MAIN
##########################

nodeID_regions = load_file(options[:input_file])
sorted_regions = get_patients_per_bin(nodeID_regions, 50000)
sorted_regions.each do |data|
	puts data.join("\t")
end
