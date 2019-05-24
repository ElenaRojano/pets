#! /usr/bin/env ruby

require 'optparse'

#################################
## METHODS
#################################
def load_pairs(file, key)
	pairsA = {}
	pairsB = {}
	File.open(file).each do |line|
		line.chomp!
		fields = line.split("\t")
		if fields.first =~ /#{key}/#.include?(key)
			save_record(pairsA, fields.last, fields.first )
		else
			save_record(pairsB, fields.last, fields.first )
		end
	end
	return pairsA, pairsB
end

def save_record(hash, key, val)
  query = hash[key]
  if query.nil?
    hash[key] = [val]
  else
    query << val
  end
end

def generate_files(n_files, output)
	files = []
	n_files.times do |n|
		files << File.open("#{output}#{n+1}.txt", 'w')
	end
	return files
end

def connect_pairs_write(pairsA, pairsB, n_files, files)
	pairsA.each do |keyA, valA|
	  valB = pairsB[keyA]
	  if !valB.nil?
	    valA.each do |vA|
	      valB.each do |vB|
	        files[rand(n_files)].puts "#{vA}\t#{keyA}\t#{vB}"
	      end
	    end
	  end
	end
end

def get_relations(pairsA, pairsB)
	relations = {}
	pairsA.each do |keyA, valA|
		valB = pairsB[keyA]
		if !valB.nil?
			valA.each do |vA|
		        	valB.each do |vB|
      					rel_key = vA + '_' + vB 
					query = relations[rel_key]
					if query.nil?
						relations[rel_key] = [keyA]
		      			else
      						query << keyA
			      		end
			        end
			end
	    	end
	end
	return relations
end

##############################
#OPTPARSE
##############################

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"

  options[:input_file] = nil
  opts.on("-i", "--input_file PATH", "Input file for create adjacency matrix") do |input_file|
    options[:input_file] = input_file
  end

  options[:key] = ''
  opts.on("-k", "--key STRING", "String to split th two groups") do |key|
    options[:key] = key
  end

  options[:output] = 'tri_'
  opts.on("-o", "--output PATH", "Output network pairs") do |output|
    options[:output] = output
  end

  options[:n_files] = 10
  opts.on("-n", "--n_files INTEGER", "Split data onto n files") do |n|
    options[:n_files] = n.to_i
  end

  options[:min_connections] = 1
  opts.on("-m", "--min_connections INTEGER", "Minimun connections to take into account a relation") do |n|
    options[:min_connections] = n.to_i
  end

end.parse!

################################
## MAIN
################################
files = generate_files(options[:n_files], options[:output])

pairsA, pairsB = load_pairs(options[:input_file], options[:key])
if options[:min_connections] == 1
	connect_pairs_write(pairsA, pairsB, options[:n_files], files)
else
	STDERR.puts "MIN. NUMBER OF CONNECTIONS = #{options[:min_connections]}"
	relations = get_relations(pairsA, pairsB)
	count = 0
	discarded = 0
	relations.each do |rel, connections|
		if connections.length >= options[:min_connections]
			fields = rel.split('_')
			connections.each do |con|
				files[rand(options[:n_files])].puts "#{fields.first}\t#{con}\t#{fields.last}"
			end
		else
			discarded += connections.length
		end
		count += connections.length
	end
	STDERR.puts "Relations: #{count}"
	STDERR.puts "Discarded: #{discarded}"
end
files.each do |f|
  f.close
end
