#! /usr/bin/env ruby

ROOT_PATH = File.dirname(__FILE__)
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'pets'))

require 'optparse'
require 'report_html'
require 'semtools'
require 'pets'


#############################################################################################
## METHODS
############################################################################################
def procces_patient_data(patient_data, hpo)
	clean_profiles = {}
	patient_data.each do |pat_id, data|
		profile = hpo.clean_profile_hard(data.first.map{|c| c.to_sym})
		clean_profiles[pat_id] = profile if !profile.empty?
	end
	return clean_profiles
end

#############################################################################################
## OPTPARSE
############################################################################################

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"

  options[:chromosome_col] = nil
  opts.on("-c", "--chromosome_col INTEGER/STRING", "Column name if header is true, otherwise 0-based position of the column with the chromosome") do |data|
    options[:chromosome_col] = data
  end

  options[:pat_id_col] = nil
  opts.on("-d", "--pat_id_col INTEGER/STRING", "Column name if header is true, otherwise 0-based position of the column with the patient id") do |data|
    options[:pat_id_col] = data
  end

  options[:end_col] = nil
  opts.on("-e", "--end_col INTEGER/STRING", "Column name if header is true, otherwise 0-based position of the column with the end mutation coordinate") do |data|
    options[:end_col] = data
  end
  
  options[:header] = true
  opts.on("-H", "--header", "File has a line header. Default true") do 
    options[:header] = false
  end

  options[:output_file] = 'report.html'
  opts.on("-o", "--output_file PATH", "Output paco file with HPO names") do |data|
    options[:output_file] = data
  end

  options[:input_file] = nil
  opts.on("-P", "--input_file PATH", "Input file with PACO extension") do |value|
    options[:input_file] = value
  end

  options[:hpo_col] = nil
  opts.on("-p", "--hpo_term_col INTEGER/STRING", "Column name if header true or 0-based position of the column with the HPO terms") do |data|
    options[:hpo_col] = data
  end

  options[:start_col] = nil
  opts.on("-s", "--start_col INTEGER/STRING", "Column name if header is true, otherwise 0-based position of the column with the start mutation coordinate") do |data|
    options[:start_col] = data
  end

  options[:hpo_separator] = '|'
  opts.on("-S", "--hpo_separator STRING", "Set which character must be used to split the HPO profile. Default '|'") do |data|
    options[:hpo_separator] = data
  end

  options[:matrix_limits] = [20, 40]
  opts.on("-L", "--matrix_limits STRING", "Number of rows and columns to show in heatmap defined as 'Nrows,Ncols'. Default 20,40") do |data|
    options[:matrix_limits] = data.split(",").map{|i| i.to_i}
  end

  options[:ref_prof] = nil
  opts.on("-r", "--ref_profile PATH", "Path to reference profile. One term code per line") do |data|
    options[:ref_prof] = File.open(data).readlines.map{|c| c.chomp.to_sym}
  end

	 opts.on_tail("-h", "--help", "Show this message") do
	    puts opts
	    exit
	  end
end.parse!

#############################################################################################
## MAIN
############################################################################################
patient_data = load_patient_cohort(options)

hpo_file = !ENV['hpo_file'].nil? ? ENV['hpo_file'] : HPO_FILE
hpo = Ontology.new
hpo.read(hpo_file)

clean_profiles = procces_patient_data(patient_data, hpo)
if !options[:ref_prof].nil?
  ref_profile = hpo.clean_profile_hard(options[:ref_prof])
else
  ref_profile = hpo.clean_profile_hard(clean_profiles.flatten.uniq)
end
hpo.load_profiles({ref: ref_profile})

similarities = hpo.compare_profiles(external_profiles: clean_profiles, sim_type: :lin, bidirectional: false)

candidate_sim_matrix, candidates, candidates_ids = get_similarity_matrix(ref_profile, similarities[:ref], clean_profiles, hpo, options[:matrix_limits].first,options[:matrix_limits].last)
candidate_sim_matrix.unshift(['HP'] + candidates_ids)

template = File.open(File.join(REPORT_FOLDER, 'similarity_matrix.erb')).read
container = {
	similarity_matrix: candidate_sim_matrix
}
report = Report_html.new(container, 'Similarity matrix')
report.build(template)
report.write(options[:output_file])

File.open(options[:output_file].gsub('.html','') +'.txt', 'w') do |f|
  similarities[:ref].each do |candidate, value|
    f.puts [candidate.to_s, value].join("\t")
  end
end
