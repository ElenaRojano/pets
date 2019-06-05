#! /usr/bin/env ruby

#Tool for calculate the AUC on PR curves.

##########################
#LIBRARIES
##########################

require 'optparse'

##########################
#METHODS
##########################

def load_file(input_file, x_val_col, y_val_col)
	coordinates = []
	File.open(input_file).each do |line|
		line.chomp!
		next if line.include?('prec') || line.include?('rec')
		info = line.split("\t")
		x_value = info[x_val_col - 1].to_f
		y_value = info[y_val_col - 1].to_f
		#STDERR.puts y_value
		coordinates << [x_value, y_value]
	end
	return coordinates.sort{|r1, r2| r1[0] <=> r2[0]}
end


def calculate_auc(pr_values)
	#pr_values = [[x, y], [x', y']...]
	x_val = 0
	y_val = 0
	total_area = 0
	pr_values.each_with_index do |xy_pair, counter|
		if counter != 0
			current_x = xy_pair[0]
			current_y = xy_pair[1]
			#puts x_val
			total_area += (x_val - current_x).abs * current_y
			#STDERR.puts total_square_area
			total_area += (x_val - current_x).abs * (y_val - current_y).abs / 2
			x_val = current_x
			y_val = current_y
		else
			x_val = xy_pair[0]
			y_val = xy_pair[1]
			total_area += x_val * y_val
			#STDERR.puts total_area
		end
		#STDERR.puts total_area
	end
	#STDERR.puts total_area
	return total_area
end

# def calculate_auc(pr_values)
# 	#pr_values = [[x, y], [x', y']...]
# 	counter = 0
# 	x_val = 0
# 	y_val = 0
# 	total_square_area = 0
# 	total_triangle_area = 0
# 	pr_values.each do |xy_pair|
# 		if counter != 0
# 			x_prime = xy_pair[0]
# 			y_prime = xy_pair[1]
# 			#puts "#{x_prime}\t#{x_val}"
# 			#puts "#{y_prime}\t#{y_val}"
# 			total_square_area += (x_val - x_prime) * y_prime
# 			total_triangle_area += (x_val - x_prime) * (y_prime - y_val) / 2
# 			x_val = x_prime
# 			y_val = y_prime
# 		else
# 			x_val = xy_pair[0]
# 			y_val = xy_pair[1]
# 			counter += 1
# 		end
# 	end
# 	total_area = total_square_area + total_triangle_area
# 	STDERR.puts total_area
# 	return total_area
# end


##########################
#OPT-PARSE
##########################

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"
  
  options[:input_file] = nil
  opts.on("-f", "--input_file PATH", "Precision-recall values file") do |input_file|
    options[:input_file] = input_file
  end

  options[:x_values] = nil
  opts.on("-x", "--:x_values INTEGER", "Set column for extracting x values") do |x_values|
    options[:x_values] = x_values.to_i
  end

  options[:y_values] = nil
  opts.on("-y", "--:y_values INTEGER", "Set column for extracting y values") do |y_values|
    options[:y_values] = y_values.to_i
  end

end.parse!

##########################
#MAIN
##########################

pr_values = load_file(options[:input_file], options[:x_values], options[:y_values])
#puts pr_values
final_area = calculate_auc(pr_values)
puts final_area