#! /usr/bin/env ruby

ROOT_PATH = File.dirname(__FILE__)
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'pets'))

require 'optparse'
require 'pets'

##################################
## METHODS
##################################

def get_data(options)
	fields2extract = get_fields2extract(options)
	field_numbers = fields2extract.values
	records = read_records(options, fields2extract, field_numbers)
end

def read_records(options, fields2extract, field_numbers) # Modified from cohort_parset
	records = []
	count = 0
	File.open(options[:input_file]).each do |line|
		line.chomp!
		if options[:header] && count == 0
			line.gsub!(/#\s*/,'') # correct comment like	headers
			field_names = line.split("\t")
			get_field_numbers2extract(field_names, fields2extract)
			field_numbers = fields2extract.values
		else
			fields = line.split("\t")
			record = field_numbers.map{|n| fields[n]}
			if fields2extract[:id_col].nil?
				id = "rec_#{count}" #generate ids
			else
				id = record.shift
			end
			record[1] = record[1].to_i 
			record[2] = record[2].to_i
			record << id
			records << record
		end
		count +=1
	end
	return records
end

def get_fields2extract(options)
	fields2extract = {}
	[:id_col, :chromosome_col, :start_col, :end_col].each do |field|
		col = options[field]
		if !col.nil?
			col = col.to_i if !options[:header]
			fields2extract[field] = col
		end
	end
	return fields2extract
end

def get_field_numbers2extract(field_names, fields2extract)
	fields2extract.each do |field, name|
		fields2extract[field] = field_names.index(name)
	end
end

##########################
#OPT-PARSER
##########################

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"

  options[:chromosome_col] = nil
  opts.on("-c", "--chromosome_col INTEGER/STRING", "Column name if header is true, otherwise 0-based position of the column with the chromosome") do |data|
    options[:chromosome_col] = data
  end

  options[:id_col] = nil
  opts.on("-d", "--id_col INTEGER/STRING", "Column name if header is true, otherwise 0-based position of the column with the id") do |data|
    options[:id_col] = data
  end

  options[:end_col] = nil
  opts.on("-e", "--end_col INTEGER/STRING", "Column name if header is true, otherwise 0-based position of the column with the end mutation coordinate") do |data|
    options[:end_col] = data
  end

  options[:header] = true
  #chr\tstart\tstop
  opts.on("-H", "--header", "Set if the file has a line header. Default true") do 
    options[:header] = false
  end

  options[:input_file] = nil
  opts.on("-i", "--input_file PATH", "Input file path") do |data|
    options[:input_file] = data
  end

  options[:reference_file] = nil
  opts.on("-r", "--reference_file PATH", "Reference file with genome annotation") do |data|
    options[:reference_file] = data
  end

  options[:output_file] = nil
  opts.on("-o", "--output_file PATH", "Output file with patient data") do |data|
    options[:output_file] = data
  end

  options[:start_col] = nil
  opts.on("-s", "--start_col INTEGER/STRING", "Column name if header is true, otherwise 0-based position of the column with the start mutation coordinate") do |data|
  	options[:start_col] = data
  end

  options[:feature_type] = nil
  opts.on("-t", "--feature_type STRING", "Keep features from reference whose are tagged with this feature type") do |data|
  	options[:feature_type] = data
  end

  options[:feature_name] = nil
  opts.on("-n", "--feature_name STRING", "Use this feature id that is present in attributes/annotation field of reference") do |data|
  	options[:feature_name] = data
  end

  opts.on_tail("-h", "--help", "Show this message") do
    puts opts
    exit
  end

end.parse!

regions = Genomic_Feature.new(get_data(options))
Genomic_Feature.add_reference(
	Reference_parser.load(
		options[:reference_file], 
		feature_type: options[:feature_type]
	)
)
gene_features = regions.get_features(attr_type: options[:feature_name])

File.open(options[:output_file], 'w') do |f|
	gene_features.each do |id, feat_ids|
		feat_ids.each do |ft_id|
			f.puts "#{id}\t#{ft_id}"
		end
	end 
end