#! /usr/bin/env ruby
#Tool to create the training file, taking as input the cluster_coords.txt file and phenotype_mutations_relations.txt

##########################
#RUBY GEMS
##########################
require 'optparse'

##########################
#METHODS
##########################

def load_cluster_file(cluster_file)
	clusters_info = {}
	File.open(cluster_file).each do |line|
		line.chomp!
		start, stop, chr, node = line.split("\t")
		clusters_info[node] = [chr, start, stop]
	end
	return clusters_info
end

def obtain_training(relations_file, clusters, filter)
	File.open(relations_file).each do |line|
		line.chomp!
		hpo, node, score = line.split("\t")
		next if score.to_f.abs <= filter
		clustersFileInfo = clusters[node]
		puts "#{clustersFileInfo.join("\t")}\t#{hpo}\t#{score}\t#{node}"
	end
end

##########################
#OPT-PARSE
##########################
options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"

  options[:cluster_file] = nil
  opts.on("-c", "--cluster_file PATH", "Input file with patient clusters") do |cluster_path|
    options[:cluster_file] = cluster_path
  end

  options[:relations_file] = nil
  opts.on("-n", "--relations_file PATH", "Input relations file from tripartite network") do |relations_file|
    options[:relations_file] = relations_file
  end

  options[:filter_association] = 0
  opts.on("-f", "--filter_minimun INTEGER", "Filter for association values") do |filter_association|
    options[:filter_association] = filter_association.to_f
  end


end.parse!

##########################
#MAIN
##########################
clusters = load_cluster_file(options[:cluster_file])
obtain_training(options[:relations_file], clusters, options[:filter_association])