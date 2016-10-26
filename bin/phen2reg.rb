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

def search4HPO(prediction_file, trainingData)
	hpoInfo = {}
	File.open(prediction_file).each do |line|
		line.chomp!
		hpo_coords = trainingData[line]
		#always check if user gives an HPO code:
		if !hpo_coords.nil?
			hpoInfo[line] = hpo_coords
		end
	end
	return hpoInfo
end

def sortByRegion(hpoInfo)
	storageHPO = {}
	hpoInfo.each do |referenceHPO, refRegions|
		hpoInfo.each do |currentHPO, currentRegions|
			next if referenceHPO == currentHPO
			refRegions.each do |refChr, refStart, refStop, refAsocVal|
				currentRegions.each do |curChr, curStart, curStop, curAsocVal|
					if refChr == curChr && 
						refStart == curStart &&
						refStop == curStop
						region = "#{curChr}:#{curStart}-#{curStop}"
						query = storageHPO[region]
						if !query.nil?
							query << currentHPO if !query.include?(currentHPO)
						else
							storageHPO[region] = [currentHPO, referenceHPO]
						end
					end
				end	
			end
		end
	end
	return storageHPO
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
  opts.on("-p", "--prediction_file PATH", "Input file with HPO codes for predicting their location") do |input_path|
    options[:prediction_file] = input_path
  end

  options[:sort_by_region] = FALSE
  opts.on("-s", "--sort_by_region", "Predict which HPOs are located in the same region") do
  options[:sort_by_region] = TRUE
  end

end.parse!


##########################
#MAIN
##########################

trainingData = load_training_file4HPO(options[:training_file])
hpoInfo = search4HPO(options[:prediction_file], trainingData)
if hpoInfo.empty?
	puts "Results not found"
elsif options[:sort_by_region] != TRUE
	hpoInfo.each do |hpo, regions|
		regions.each do |region|	
			puts "#{hpo}\t#{region.join("\t")}"
		end
	end
elsif options[:sort_by_region] == TRUE	
	fullInfo = sortByRegion(hpoInfo)
	fullInfo.each do |position, hpoCodes|
		puts "#{position}\t#{hpoCodes.join(",")}"		
	end
end
