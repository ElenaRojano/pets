#! /usr/bin/env ruby
#Code for calculating F measure for precision-recall curves 

require 'optparse'

##########################
#METHODS
##########################

def load_pr_data(file)
	counter = 0
	pr_data = {}
	File.open(file).each do |line|
		if counter > 0
			line.chomp!
			cutoff, prec, rec, meth = line.split("\t")
			query = pr_data[meth]
			pr_info = [cutoff.to_f, prec.to_f, rec.to_f]
			if query.nil?
				pr_data[meth] = [pr_info]
			else
				query << pr_info
			end
		end
		counter += 1
	end
	return pr_data
end

def calculate_youden(pr_data)
	#the max f_measure is the best cutoff
	best_cutoffs = []
	pr_data.each do |meth, pr_values|
		max_f_measure = 0
		best_cutoff = 0
		#next if meth != 'cosine'
		pr_values.each do |cutoff, prec, rec|
			f_measure = 2 * prec * rec / (prec + rec)
			if max_f_measure < f_measure
				max_f_measure = f_measure
				best_cutoff = cutoff
			end
		end
		best_cutoffs << [meth, best_cutoff]
	end
	return best_cutoffs
end

##########################
#OPT-PARSE
##########################
options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"

  options[:input_file] = nil
  opts.on("-f", "--input_file PATH", "Input file with precision-recall values") do |input_file|
    options[:input_file] = input_file
  end

end.parse!

##########################
#MAIN
##########################

pr_data = load_pr_data(options[:input_file])
best_cutoffs = calculate_youden(pr_data)
#puts best_cutoffs.inspect
best_cutoffs.each do |cutoffs|
	puts cutoffs.join("\t")
end