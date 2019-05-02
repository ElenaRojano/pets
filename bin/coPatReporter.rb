#! /usr/bin/env ruby

ROOT_PATH = File.dirname(__FILE__)
REPORT_FOLDER = File.expand_path(File.join(ROOT_PATH, '..', 'templates'))
EXTERNAL_DATA = File.expand_path(File.join(ROOT_PATH, '..', 'external_data'))
HPO2NAME_DICTIONARY = File.join(EXTERNAL_DATA, 'hpo2name.txt')
IC_FILE = File.join(EXTERNAL_DATA, 'uniq_hpo_with_CI.txt')
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'gephepred'))

require 'optparse'
require 'csv'
require 'generalMethods.rb'
require 'report_html'

##########################
#METHODS
##########################
HPOS = 0
CHR = 1
START = 2
STOP = 3

def load_patient_cohort(options)
	patient_data = {}
	count = 0
	fields2extract = get_fields2extract(options)
	field_numbers = fields2extract.values
  File.open(options[:input_file]).each do |line|
    line.chomp!
    if options[:header] && count == 0
      line.gsub!(/#\s*/,'') # correct comment like  headers
      field_names = line.split("\t")
      get_field_numbers2extract(field_names, fields2extract)
      field_numbers = fields2extract.values
    else
      fields = line.split("\t")
      pat_record = field_numbers.map{|n| fields[n]}
      if fields2extract[:pat_id_col].nil?
        pat_id = "pat_#{count}" #generate ids
      else
        pat_id = pat_record.shift + "_i#{count}" # make sure that ids are uniq
      end
      patient_data[pat_id] = pat_record
    end
    count +=1
  end
  options[:pat_id_col] = 'generated' if fields2extract[:pat_id_col].nil?
	return patient_data
end 

def get_fields2extract(options)
	fields2extract = {}
	[:pat_id_col, :hpo_col, :chromosome_col, :start_col, :end_col].each do |field|
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

def format_patient_data(patient_data, options, name2code_dictionary, hpo_storage, hpo_parent_child_relations)
  all_hpo = []
  suggested_childs = {}
  patient_data.each do |pat_id, patient_record|
    string_hpos, chr, start, stop = patient_record
    hpos = string_hpos.split(options[:hpo_separator])
    translate_hpo_names2codes(hpos, name2code_dictionary, pat_id) if options[:hpo_names]
    suggested_childs[pat_id] = check_hpo_codes(hpos, hpo_storage, hpo_parent_child_relations, pat_id)
    all_hpo.concat(hpos)
    patient_record[HPOS] = hpos
    patient_record[START] = start.to_i if !start.nil?
    patient_record[STOP] = stop.to_i if !stop.nil?
  end
  return all_hpo.uniq, suggested_childs
end

def translate_hpo_names2codes(hpos, hpo_dictionary, pat_id)
  hpo_codes = []
  hpos.each_with_index do |hpo_name, i|
    hpo_code = hpo_dictionary[hpo_name]
    if hpo_code.nil?
      STDERR.puts "WARNING: patient #{pat_id} has the unknown hpo NAME '#{hpo_name}'. Rejected."
    else
      hpo_codes << hpo_code
    end
  end
  hpos.clear
  hpos.concat(hpo_codes)
end

def check_hpo_codes(hpos, hpo_storage, hpo_parent_child_relations, pat_id)
  more_specific_hpo = []
  hpos.each_with_index do |hpo_code, i|
    hpo_data = hpo_storage[hpo_code]
    if hpo_data.nil?
      hpos[i] = nil
      STDERR.puts "WARNING: patient #{pat_id} has the unknown hpo CODE '#{hpo_code}'. Rejected."
    else
      main_hpo_code, name = hpo_data
      hpos[i] = main_hpo_code # change from alternate hpo codes to the main ones
      childs = hpo_parent_child_relations[main_hpo_code]
      if childs.nil?
        specific_childs = []
      else
        specific_childs = childs.last
      end
      more_specific_hpo << [[main_hpo_code, name], specific_childs]
    end
  end
  hpos.compact!
  return more_specific_hpo
end

def generate_patient_hpo_matrix(patient_data, cohort_hpos)
  matrix = []
  n = cohort_hpos.length
  patient_data.each do |pat_id, patient_record|
    pat_hpos = patient_record[HPOS]
    vector = Array.new(n, 0)
    pat_hpos.each do |hpo|
      vector[cohort_hpos.index(hpo)] = 1
    end
    matrix << vector
  end
  return matrix
end

def write_matrix_for_R(matrix, x_names, y_names, file)
  File.open(file, 'w') do |f|
    f.puts x_names.join("\t")
    matrix.each_with_index do |row, i|
      f.puts [y_names[i]].concat(row).join("\t")
    end
  end
end

def process_clustered_patients(options, clustered_patients, patient_data) # get ic and chromosomes
  ic_file = ENV['ic_file']
  ic_file = IC_FILE if ic_file.nil?
  phenotype_ic = load_hpo_ci_values(ic_file)
  all_ics = []
  top_cluster_phenotypes = []
  cluster_data_by_chromosomes = []
  multi_chromosome_patients = 0
  processed_clusters = 0
  clustered_patients.sort_by{|cl_id, pat_ids| pat_ids.length }.reverse.each do |cluster_id, patient_ids|
    patient_number = patient_ids.length
    next if patient_number == 1
    chrs = Hash.new(0)
    all_phens = []
    profile_ics = []
    patient_ids.each do |pat_id|
      patient = patient_data[pat_id]
      phenotypes = patient[HPOS]       
      profile_ics << get_profile_ic(phenotypes, phenotype_ic)
      #optional
      all_phens << phenotypes if processed_clusters < options[:clusters2show_detailed_phen_data]      
      chrs[patient[CHR]] += 1 if !options[:chromosome_col].nil?
    end
    top_cluster_phenotypes << all_phens if processed_clusters < options[:clusters2show_detailed_phen_data]
    all_ics << profile_ics
    if !options[:chromosome_col].nil?
      multi_chromosome_patients += patient_number if chrs.length > 1
      chrs.each do |chr, count|
        cluster_data_by_chromosomes << [cluster_id, patient_number, chr, count]
      end
    end
    processed_clusters += 1
  end
  return all_ics, cluster_data_by_chromosomes, top_cluster_phenotypes, multi_chromosome_patients
end

def get_profile_ic(hpo_names, phenotype_ic)
  ic = 0
  profile_length = 0
  hpo_names.each do |hpo_id|
    hpo_ic = phenotype_ic[hpo_id] 
    ic +=  hpo_ic if !hpo_ic.nil?
    profile_length += 1
  end
  profile_length = 1 if profile_length == 0
  return ic/profile_length
end

def write_cluster_ic_data(all_ics, cluster_ic_data_file, limit)
  File.open(cluster_ic_data_file, 'w') do |f|
    f.puts %w[cluster_id ic].join("\t")
    all_ics.each_with_index do |cluster_ics, i|
      break if i == limit
      cluster_length = cluster_ics.length
      cluster_ics.each do |clust_ic|
        f.puts "#{cluster_length}_#{i}\t#{clust_ic}"
      end
    end
  end
end

def write_cluster_chromosome_data(cluster_data, cluster_chromosome_data_file, limit)
  File.open(cluster_chromosome_data_file, 'w') do |f|
    f.puts %w[cluster_id chr count].join("\t")
    index = 0
    last_id = cluster_data.first.first
    cluster_data.each do |cluster_id, patient_number, chr, count|
      index += 1 if cluster_id != last_id 
      break if index == limit
      f.puts ["#{patient_number}_#{index}", chr, count].join("\t")
      last_id = cluster_id
    end
  end
end

def get_hpo_profile(patient_data)
  hpo_profiles = []
  patient_data.each do |pat_id, pat_data|
      hpo_profiles << pat_data[HPOS]
  end
  return hpo_profiles
end

def get_summary_stats(patient_data, cohort_hpos, all_hpo_profiles)
  stats = []
  stats << ['Unique HPOs', cohort_hpos.length]
  n_pat = patient_data.length
  stats << ['Patient number', n_pat]
  all_hpo_prof_lengths = all_hpo_profiles.map{|p| p.length}.sort
  stats << ['HPOs per patient (average)', all_hpo_prof_lengths.inject(0){|sum, n| sum + n}.fdiv(n_pat)]
  hpo_pat90 = nil
  rate = 0
  count = 0
  while rate <= 0.1
    hpo_pat90 = all_hpo_prof_lengths[count+1]
    rate = count.fdiv(n_pat)
    count += 1
  end
  stats << ['HPOs for patient in percentile 90', hpo_pat90]
  return stats
end

def hpo_stats(all_hpo_profiles)
  stats = Hash.new(0)
  all_hpo_profiles.each do |profile|
    profile.each do |hpo|
      stats[hpo] += 1
    end
  end
  n_profiles = all_hpo_profiles.length
  hpo_stats = []
  stats.each do |hpo, count|
    hpo_stats << [hpo, count.fdiv(n_profiles)*100]
  end
  hpo_stats.sort!{|h1, h2| h2[1] <=> h1[1]}
  return hpo_stats[0..20]
end

def translate_hpo_codes2names(all_hpo_profiles, hpo_storage)
  all_hpo_profiles.each do |profile|
    profile.each_with_index do |hpo, i|
      hpo_data = hpo_storage[hpo]
      if hpo_data.nil?
        STDERR.puts "WARNING: hpo code '#{hpo}' not exists."
      else
        profile[i] = hpo_data[1]
      end
    end
  end
end

def write_detailed_hpo_profile_evaluation(suggested_childs, detailed_profile_evaluation_file, summary_stats)
  hpo_count = 0
  parent_hpo_count = 0
  CSV.open(detailed_profile_evaluation_file, "wb") do |csv|
    suggested_childs.each do |pat_id, suggestions|
      warning = nil
      warning = 'WARNING: Very few phenotypes' if suggestions.length < 4
      csv << ["PATIENT #{pat_id}", "#{warning}"]
      csv << ["CURRENT PHENOTYPES", "PUTATIVE MORE SPECIFIC PHENOTYPES"]
      suggestions.each do |parent, childs|
        hpo_count += 1
        parent_code, parent_name = parent
        if childs.empty?
          csv << ["#{parent_name} (#{parent_code})", '-']
        else
          parent_hpo_count += 1
          parent_writed = false
          childs.each do |child_code, child_name|
            if !parent_writed
              parent_field = "#{parent_name} (#{parent_code})"
              parent_writed = true
            else
              parent_field = ""
            end
            csv << [parent_field, "#{child_name} (#{child_code})"]
          end
        end
      end
      csv << ["", ""]
    end
  end
  summary_stats << ['Percentage of defined HPOs that have more specific childs', parent_hpo_count.fdiv(hpo_count) * 100]
end

##########################
#OPT-PARSER
##########################

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"

  options[:clusters2show_detailed_phen_data] = 3
  opts.on("-C", "--clusters2show INTEGER", "How many patient clusters are show in detailed phenotype cluster data section. Default 3") do |data|
    options[:clusters2show_detailed_phen_data] = data.to_i
  end

  options[:chromosome_col] = nil
  opts.on("-c", "--chromosome_col INTEGER/STRING", "Column name if header is true, otherwise 0-based position of the column with the chromosome") do |data|
    options[:chromosome_col] = data
  end

  options[:pat_id_col] = nil
  opts.on("-d", "--pat_id_col INTEGER/STRING", "Column name if header is true, otherwise 0-based position of the column with the patient id") do |data|
    options[:pat_id_col] = data
  end

  options[:excluded_hpo] = nil
  opts.on("-E", "--excluded_hpo PATH", "List of HPO phenotypes to exclude (low informative)") do |excluded_hpo|
    options[:excluded_hpo] = excluded_hpo
  end

  options[:end_col] = nil
  opts.on("-e", "--end_col INTEGER/STRING", "Column name if header is true, otherwise 0-based position of the column with the end mutation coordinate") do |data|
    options[:end_col] = data
  end

  options[:clusters2graph] = 30
  opts.on("-g", "--clusters2graph INTEGER", "How may patient clusters are plotted in cluster plots. Default 30") do |data|
    options[:clusters2graph] = data.to_i
  end

  options[:header] = true
  #chr\tstart\tstop
  opts.on("-H", "--header", "Set if the file has a line header. Default true") do 
    options[:header] = false
  end

  options[:input_file] = nil
  opts.on("-i", "--input_file PATH", "Input file with patient data") do |data|
    options[:input_file] = data
  end

  options[:hpo_names] = false
  opts.on("-n", "--hpo_names", "Define if the input HPO are human readable names. Default false") do
    options[:hpo_names] = true
  end

  options[:output_file] = nil
  opts.on("-o", "--output_file PATH", "Output file with patient data") do |data|
    options[:output_file] = data
  end

  options[:hpo_file] = nil
  opts.on("-P", "--hpo_file PATH", "Input HPO file for extracting HPO codes") do |value|
    options[:hpo_file] = value
  end

  options[:hpo_col] = nil
  opts.on("-p", "--hpo_term_col INTEGER/STRING", "Column name if header true or 0-based position of the column with the HPO terms") do |data|
  	options[:hpo_col] = data
  end

  options[:hpo_separator] = '|'
  opts.on("-S", "--hpo_separator STRING", "Set which character must be used to split the HPO profile. Default '|'") do |data|
  	options[:hpo_separator] = data
  end

  options[:start_col] = nil
  opts.on("-s", "--start_col INTEGER/STRING", "Column name if header is true, otherwise 0-based position of the column with the start mutation coordinate") do |data|
  	options[:start_col] = data
  end

end.parse!


##########################
#MAIN
##########################
output_folder = File.dirname(options[:output_file])
detailed_profile_evaluation_file = File.join(output_folder, 'detailed_hpo_profile_evaluation.csv')
temp_folder = File.join(output_folder, 'temp')
matrix_file = File.join(temp_folder, 'pat_hpo_matrix.txt')
clustered_patients_file = File.join(temp_folder, 'cluster_asignation')
cluster_ic_data_file = File.join(temp_folder, 'cluster_ic_data.txt')
cluster_chromosome_data_file = File.join(temp_folder, 'cluster_chromosome_data.txt')
Dir.mkdir(temp_folder) if !File.exists?(temp_folder)

# LOAD HPO DATA
#-------------------------

# #load hpo dictionaries
hpo_black_list = load_hpo_black_list(options[:excluded_hpo])
hpo_storage = load_hpo_file(options[:hpo_file], hpo_black_list)
hpo_parent_child_relations = get_child_parent_relations(hpo_storage)
name2code_dictionary = create_hpo_dictionary(hpo_storage) if options[:hpo_names]


# hpo_dictionary_file = ENV['hpo2name_file']
# hpo_dictionary_file = HPO2NAME_DICTIONARY if hpo_dictionary_file.nil?
# hpo_storage = load_hpo_metadata(hpo_dictionary_file)
# hpo_parent_child_relations = inverse_hpo_metadata(hpo_storage)
# name2code_dictionary = {}
# name2code_dictionary = load_hpo_dictionary_name2code(hpo_dictionary_file) if options[:hpo_names]

patient_data = load_patient_cohort(options)
cohort_hpos, suggested_childs = format_patient_data(patient_data, options, name2code_dictionary, hpo_storage, hpo_parent_child_relations)
pat_hpo_matrix = generate_patient_hpo_matrix(patient_data, cohort_hpos)
write_matrix_for_R(pat_hpo_matrix, cohort_hpos, patient_data.keys, matrix_file)

system("get_clusters.R #{matrix_file} #{temp_folder}") if !File.exists?(clustered_patients_file)
clustered_patients = load_clustered_patients(clustered_patients_file)
all_ics, cluster_data_by_chromosomes, top_cluster_phenotypes, multi_chromosome_patients = process_clustered_patients(options, clustered_patients, patient_data)
write_cluster_ic_data(all_ics, cluster_ic_data_file, options[:clusters2graph])
system("plot_boxplot.R #{cluster_ic_data_file} #{temp_folder} cluster_id ic 'Cluster size/id' 'Information coefficient'")
# Process.exit

if !options[:chromosome_col].nil?
  write_cluster_chromosome_data(cluster_data_by_chromosomes, cluster_chromosome_data_file, options[:clusters2graph])
  system("plot_scatterplot.R #{cluster_chromosome_data_file} #{temp_folder} cluster_id chr count 'Cluster size/id' 'Chromosome' 'Patients'")
end
all_hpo_profiles = get_hpo_profile(patient_data)
translate_hpo_codes2names(all_hpo_profiles, hpo_storage)
summary_stats = get_summary_stats(patient_data, cohort_hpos, all_hpo_profiles)
write_detailed_hpo_profile_evaluation(suggested_childs, detailed_profile_evaluation_file, summary_stats)
summary_stats << ['Number of clusters with mutations accross > 1 chromosomes', multi_chromosome_patients] if !options[:chromosome_col].nil?
hpo_stats = hpo_stats(all_hpo_profiles)

#report
container = {
  :temp_folder => temp_folder,
  :top_cluster_phenotypes => top_cluster_phenotypes.length,
  :summary_stats => summary_stats,
  :hpo_stats => hpo_stats
 }
top_cluster_phenotypes.each_with_index do |cluster, i|
  clust_pr = cluster.map{|pr| [pr.join(', ')] }
  container["clust_#{i}"] = clust_pr
end
template = File.open(File.join(REPORT_FOLDER, 'cohort_report.erb')).read
report = Report_html.new(container, 'Cohort quality report')
report.build(template)
report.write(options[:output_file]+'.html')

