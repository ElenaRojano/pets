#! /usr/bin/env ruby

ROOT_PATH = File.dirname(__FILE__)
REPORT_FOLDER = File.expand_path(File.join(ROOT_PATH, '..', 'templates'))
EXTERNAL_DATA = File.expand_path(File.join(ROOT_PATH, '..', 'external_data'))
EXTERNAL_CODE = File.expand_path(File.join(ROOT_PATH, '..', 'external_code'))
HPO_FILE = File.join(EXTERNAL_DATA, 'hp.obo')
IC_FILE = File.join(EXTERNAL_DATA, 'uniq_hpo_with_CI.txt')
CHR_SIZE = File.join(EXTERNAL_DATA, 'chromosome_sizes_hg19.txt')
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'pets'))
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'pets', 'simnetome', 'lib'))

require 'optparse'
require 'csv'
require 'generalMethods.rb'
require 'coPatReporterMethods.rb'
require 'report_html'
require 'ontology'

##########################
#METHODS
##########################
HPOS = 0
CHR = 1
START = 2
STOP = 3

def format_patient_data(patient_data, options, hpo)
  all_hpo = []
  rejected_hpos = []
  suggested_childs = {}
  total_terms = 0
  terms_with_more_specific_childs = 0
  patient_data.each do |pat_id, patient_record|
    string_hpos, chr, start, stop = patient_record
    hpos = string_hpos.split(options[:hpo_separator])

    if options[:hpo_names]
      hpos, pat_rejected_hpos = hpo.translate_names2codes(hpos)
      if !pat_rejected_hpos.empty?
        STDERR.puts "WARNING: patient #{pat_id} has the unknown hpo NAMES '#{pat_rejected_hpos.join(',')}'. Rejected."
        rejected_hpos.concat(pat_rejected_hpos)
      end
    end

    hpos, pat_rejected_hpos = hpo.check_codes(hpos)
    if !pat_rejected_hpos.empty?
      STDERR.puts "WARNING: patient #{pat_id} has the unknown hpo CODES '#{pat_rejected_hpos.join(',')}'. Rejected."
      rejected_hpos.concat(pat_rejected_hpos)
    end
    total_terms += hpos.length
    more_specific_childs = hpo.get_more_specific_childs_table(hpos)
    terms_with_more_specific_childs += more_specific_childs.select{|hpo_record| !hpo_record.last.empty?}.length #Exclude phenotypes with no childs
    suggested_childs[pat_id] = more_specific_childs  
    all_hpo.concat(hpos)
    patient_record[HPOS] = hpos
    patient_record[START] = start.to_i if !start.nil?
    patient_record[STOP] = stop.to_i if !stop.nil?
  end
  return all_hpo.uniq, suggested_childs, rejected_hpos.uniq, terms_with_more_specific_childs.fdiv(total_terms)
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

def process_clustered_patients(options, clustered_patients, patient_data, hpo) # get ic and chromosomes
  if options[:ic_stats]
    ic_file = ENV['ic_file']
    ic_file = IC_FILE if ic_file.nil?
    phenotype_ic = load_hpo_ci_values(ic_file)
  else
    phenotype_ic = compute_IC_values(patient_data, patient_data.length)
  end
  all_ics = []
  top_cluster_phenotypes = []
  cluster_data_by_chromosomes = []
  multi_chromosome_patients = 0
  processed_clusters = 0
  clustered_patients.sort_by{|cl_id, pat_ids| pat_ids.length }.reverse.each do |cluster_id, patient_ids|
    num_of_patients = patient_ids.length
    next if num_of_patients == 1
    chrs = Hash.new(0)
    all_phens = []
    profile_ics = []
    patient_ids.each do |pat_id|
      patient = patient_data[pat_id]
      phenotypes = patient[HPOS]
      profile_ics << get_profile_ic(phenotypes, phenotype_ic)
      #optional
      if processed_clusters < options[:clusters2show_detailed_phen_data]
        phen_names, rejected_codes = hpo.translate_codes2names(phenotypes)
        all_phens << phen_names
      end
      chrs[patient[CHR]] += 1 if !options[:chromosome_col].nil?
    end
    top_cluster_phenotypes << all_phens if processed_clusters < options[:clusters2show_detailed_phen_data]
    all_ics << profile_ics
    # STDERR.puts [cluster_id, num_of_patients, chr, count].inspect
    if !options[:chromosome_col].nil?
      multi_chromosome_patients += num_of_patients if chrs.length > 1
      chrs.each do |chr, count|
        cluster_data_by_chromosomes << [cluster_id, num_of_patients, chr, count]
      end
    end
    processed_clusters += 1
  end
  # STDERR.puts cluster_data_by_chromosomes.inspect
  return all_ics, cluster_data_by_chromosomes, top_cluster_phenotypes, multi_chromosome_patients
end

def get_profile_ic(hpo_names, phenotype_ic)
  ic = 0
  profile_length = 0
  hpo_names.each do |hpo_id|
    hpo_ic = phenotype_ic[hpo_id]
    # STDERR.puts phenotype_ic.inspect
    ic += hpo_ic if !hpo_ic.nil?
    profile_length += 1
  end
  profile_length = 1 if profile_length == 0
  return ic.fdiv(profile_length)
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
    last_id = cluster_data.first.first unless cluster_data.empty?
    cluster_data.each do |cluster_id, patient_number, chr, count|
      index += 1 if cluster_id != last_id 
      break if index == limit
      f.puts ["#{patient_number}_#{index}", chr, count].join("\t")
      last_id = cluster_id
    end
  end
end

def write_coverage_data(coverage_to_plot, coverage_to_plot_file)
  File.open(coverage_to_plot_file, 'w') do |f|
    coverage_to_plot.each do |chr, position, freq|
     f.puts "#{chr}\t#{position}\t#{freq}"
   end
  end
end

def get_uniq_hpo_profiles(patient_data) # To aviod duplications due to more one mutation in the same patient
  hpo_profiles = []
  parsed_pats = []
  patient_data.each do |variant_id, patient_rec|
    pat_id, count = variant_id.split('_i')
    if parsed_pats.include?(pat_id)
      next
    else
      parsed_pats << pat_id
      hpo_profiles << patient_rec[HPOS]
    end
  end
  return hpo_profiles
end

def get_patient_ids(patient_data) # To aviod duplications due to more one mutation in the same patient
  ids = []
  patient_ids = patient_data.keys
  patient_ids.each do |pat_id|
    id, count = pat_id.split('_i')
    ids << id
  end
  return  ids.uniq
end

def get_summary_stats(patient_data, cohort_hpos, hpo)
  stats = []
  stats << ['Unique HPOs', cohort_hpos.length]
  stats << ['Number of patients in the cohort', get_patient_ids(patient_data).length]
  stats << ['HPOs per patient (average)', hpo.get_profile_mean_length]
  stats << ['HPOs for patient in percentile 90', hpo.get_profile_length_at_percentile(perc=90)]
  return stats
end

def cluster_patients(patient_data, cohort_hpos, matrix_file, clustered_patients_file)
  pat_hpo_matrix = generate_patient_hpo_matrix(patient_data, cohort_hpos)
  write_matrix_for_R(pat_hpo_matrix, cohort_hpos, patient_data.keys, matrix_file)
  system("#{File.join(EXTERNAL_CODE, 'get_clusters.R')} #{matrix_file} #{clustered_patients_file}") if !File.exists?(clustered_patients_file)
  clustered_patients = load_clustered_patients(clustered_patients_file)
  return(clustered_patients)
end

def get_profile_ontology_distribution_tables(hpo)
  cohort_ontology_levels = hpo.get_ontology_levels_from_profiles(uniq=false)
  uniq_cohort_ontology_levels = hpo.get_ontology_levels_from_profiles
  hpo_ontology_levels = hpo.term_level
  total_ontology_terms = hpo_ontology_levels.values.flatten.length
  total_cohort_terms = cohort_ontology_levels.values.flatten.length
  total_uniq_cohort_terms = uniq_cohort_ontology_levels.values.flatten.length

  ontology_levels = []
  distribution_percentage = []
  hpo_ontology_levels.each do |level, terms|
    cohort_terms = cohort_ontology_levels[level]
    uniq_cohort_terms = uniq_cohort_ontology_levels[level]
    if cohort_terms.nil? || uniq_cohort_terms.nil?
      num = 0
      u_num = 0
    else
      num = cohort_terms.length
      u_num = uniq_cohort_terms.length
    end
    ontology_levels << [level, terms.length, num]
    distribution_percentage << [
      level,
      (terms.length.fdiv(total_ontology_terms)*100).round(3),
      (num.fdiv(total_cohort_terms)*100).round(3),
      (u_num.fdiv(total_uniq_cohort_terms)*100).round(3)
    ]
  end
  ontology_levels.sort! { |x,y| x.first <=> y.first }
  distribution_percentage.sort! { |x,y| x.first <=> y.first }
  ontology_levels.unshift ["level", "ontology", "cohort"]
  distribution_percentage.unshift ["level", "ontology", "weighted cohort", "uniq terms cohort"]
  return ontology_levels, distribution_percentage
end


def write_detailed_hpo_profile_evaluation(suggested_childs, detailed_profile_evaluation_file, summary_stats)
  CSV.open(detailed_profile_evaluation_file, "wb") do |csv|
    suggested_childs.each do |pat_id, suggestions|
      warning = nil
      warning = 'WARNING: Very few phenotypes' if suggestions.length < 4
      csv << ["PATIENT #{pat_id}", "#{warning}"]
      csv << ["CURRENT PHENOTYPES", "PUTATIVE MORE SPECIFIC PHENOTYPES"]
      suggestions.each do |parent, childs|
        parent_code, parent_name = parent
        if childs.empty?
          csv << ["#{parent_name} (#{parent_code})", '-']
        else
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
end



def write_arrays4scatterplot(onto_ic, freq_ic, hpo_ic_file)
  File.open(hpo_ic_file, 'w') do |f|
    f.puts "OntoIC\tFreqIC"
    onto_ic.each_with_index do |value,i|
        f.puts [value, freq_ic[i]].join("\t")
    end
  end
end

##########################
#OPT-PARSER
##########################

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"

  options[:coverage_analysis] = true
  opts.on("-a", "--coverage_analysis", "Deactivate genome coverage analysis. Default true") do 
    options[:coverage_analysis] = false
  end

  options[:bin_size] = 50000
  opts.on("-b", "--bin_size INTEGER", "Maximum number of bins to plot the coverage") do |data|
    options[:bin_size] = data.to_i
  end

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

  options[:patients_filter] = 2
  opts.on("-f", "--patients_filter INTEGER", "Minimum number of patients sharing SORs. Default 0") do |data|
    options[:patients_filter] = data.to_i
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

  options[:ic_stats] = false
  opts.on("-t", "--ic_stats", "Use internal IC stats. Default false") do
    options[:ic_stats] = true
  end

  opts.on_tail("-h", "--help", "Show this message") do
    puts opts
    exit
  end

end.parse!


##########################
#MAIN
##########################
output_folder = File.dirname(options[:output_file])
detailed_profile_evaluation_file = File.join(output_folder, 'detailed_hpo_profile_evaluation.csv')
temp_folder = File.join(output_folder, 'temp')
matrix_file = File.join(temp_folder, 'pat_hpo_matrix.txt')
hpo_ic_file = File.join(temp_folder, 'hpo_ic.txt')
hpo_profile_ic_file = File.join(temp_folder, 'hpo_ic.txt')
clustered_patients_file = File.join(temp_folder, 'cluster_asignation')
cluster_ic_data_file = File.join(temp_folder, 'cluster_ic_data.txt')
cluster_chromosome_data_file = File.join(temp_folder, 'cluster_chromosome_data.txt')
coverage_to_plot_file = File.join(temp_folder, 'coverage_data.txt')
sor_coverage_to_plot_file = File.join(temp_folder, 'sor_coverage_data.txt')
# cnvs_lenght_to_plot_file = File.join(temp_folder, 'cnvs_lenght.txt')
Dir.mkdir(temp_folder) if !File.exists?(temp_folder)

hpo_black_list = []
hpo_file = ENV['hpo_file']
hpo_file = HPO_FILE if hpo_file.nil?

hpo = Ontology.new
hpo.load_black_list(options[:excluded_hpo]) if !options[:excluded_hpo].nil?
hpo.load_data(hpo_file)

patient_data, patient_number = load_patient_cohort(options)
cohort_hpos, suggested_childs, rejected_hpos, fraction_terms_specific_childs = format_patient_data(patient_data, options, hpo)
hpo.load_profiles(get_uniq_hpo_profiles(patient_data))

ontology_levels, distribution_percentage = get_profile_ontology_distribution_tables(hpo)

onto_ic, freq_ic = hpo.get_ic_by_onto_and_freq(hpo_file)
onto_ic_profile, freq_ic_profile = hpo.get_ic_profile_by_onto_and_freq

clustered_patients = cluster_patients(patient_data, cohort_hpos, matrix_file, clustered_patients_file) 
all_ics, cluster_data_by_chromosomes, top_cluster_phenotypes, multi_chromosome_patients = process_clustered_patients(options, clustered_patients, patient_data, hpo)

summary_stats = get_summary_stats(patient_data, cohort_hpos, hpo)
summary_stats << ['Percentage of defined HPOs that have more specific childs', (fraction_terms_specific_childs * 100).round(4)]

# Move code 'Percentage of defined HPOs that have more specific childs' outside the next function
hpo_stats = hpo.get_term_frequency_from_profiles(names=true)[0..20]
summary_stats << ['Number of unknown phenotypes', rejected_hpos.length]

all_cnvs_length = []
if !options[:chromosome_col].nil?
  summary_stats << ['Number of clusters with mutations accross > 1 chromosomes', multi_chromosome_patients]
  
  #----------------------------------
  # Prepare data to plot coverage
  #----------------------------------
  if options[:coverage_analysis]
    processed_patient_data = process_patient_data(patient_data)
    cnv_sizes = []
    processed_patient_data.each do |chr, metadata|
      metadata.each do |patientID, start, stop|
        cnv_sizes << stop - start
      end
    end
    cnv_size_average = cnv_sizes.inject{ |sum, el| sum + el }.fdiv(cnv_sizes.length.to_f)
    patients_by_cluster, sors = generate_cluster_regions(processed_patient_data, 'A', 0)
    total_patients_sharing_sors = []
    all_patients = patients_by_cluster.keys
    all_patients.each do |identifier|
      total_patients_sharing_sors << identifier.split('_i').first
    end
    all_cnvs_length = get_cnvs_length(patient_data)
    
    ###1. Process CNVs
    raw_coverage, n_cnv, nt, pats_per_region = calculate_coverage(sors)
    summary_stats << ['Number of genome windows', n_cnv]
    summary_stats << ['Nucleotides affected by mutations', nt]
    summary_stats << ['Patient average per region', pats_per_region.round(4)]
    summary_stats << ['CNV size average', cnv_size_average.round(4)]
    coverage_to_plot = get_final_coverage(raw_coverage, options[:bin_size])

    ###2. Process SORs
    raw_sor_coverage, n_sor, nt, pats_per_region = calculate_coverage(sors, options[:patients_filter] - 1)
    summary_stats << ["Number of patients with at least 1 SOR", total_patients_sharing_sors.uniq.length]
    summary_stats << ["Number of SORs with >= #{options[:patients_filter]} patients", n_sor]
    summary_stats << ['Nucleotides affected by mutations', nt]
    # summary_stats << ['Patient average per region', pats_per_region]
    sor_coverage_to_plot = get_final_coverage(raw_sor_coverage, options[:bin_size])

    all_sor_length = get_sor_length_distribution(raw_sor_coverage)  
  end
end

#----------------------------------
# Write files for report
#----------------------------------
write_detailed_hpo_profile_evaluation(suggested_childs, detailed_profile_evaluation_file, summary_stats)
write_arrays4scatterplot(onto_ic.values, freq_ic.values, hpo_ic_file)
write_arrays4scatterplot(onto_ic_profile, freq_ic_profile, hpo_profile_ic_file)
system("#{File.join(EXTERNAL_CODE, 'plot_scatterplot_simple.R')} #{hpo_ic_file} #{File.join(temp_folder, 'hpo_ics.pdf')} 'OntoIC' 'FreqIC' 'HP Ontology IC' 'HP Frequency based IC'")
system("#{File.join(EXTERNAL_CODE, 'plot_scatterplot_simple.R')} #{hpo_profile_ic_file} #{File.join(temp_folder, 'hpo_profile_ics.pdf')} 'OntoIC' 'FreqIC' 'HP Ontology Profile IC' 'HP Frequency based Profile IC'")

write_cluster_ic_data(all_ics, cluster_ic_data_file, options[:clusters2graph])
system("#{File.join(EXTERNAL_CODE, 'plot_boxplot.R')} #{cluster_ic_data_file} #{temp_folder} cluster_id ic 'Cluster size/id' 'Information coefficient'")

if !options[:chromosome_col].nil?
  write_cluster_chromosome_data(cluster_data_by_chromosomes, cluster_chromosome_data_file, options[:clusters2graph])
  system("#{File.join(EXTERNAL_CODE, 'plot_scatterplot.R')} #{cluster_chromosome_data_file} #{temp_folder} cluster_id chr count 'Cluster size/id' 'Chromosome' 'Patients'")
  if options[:coverage_analysis]
    ###1. Process CNVs
    write_coverage_data(coverage_to_plot, coverage_to_plot_file)
    cmd = "#{File.join(EXTERNAL_CODE, 'plot_area.R')} -d #{coverage_to_plot_file} -o #{temp_folder}/coverage_plot -x V2 -y V3 -f V1 -H -m #{CHR_SIZE} -t CNV"
    system(cmd)
    ###2. Process SORs
    write_coverage_data(sor_coverage_to_plot, sor_coverage_to_plot_file)
    system("#{File.join(EXTERNAL_CODE, 'plot_area.R')} -d #{sor_coverage_to_plot_file} -o #{temp_folder}/sor_coverage_plot -x V2 -y V3 -f V1 -H -m #{CHR_SIZE} -t SOR")    
  end
end
#----------------------------------
# Report
#----------------------------------
total_patients = 0
new_cluster_phenotypes = {}
phenotypes_frequency = Hash.new(0)
top_cluster_phenotypes.each_with_index do |cluster, clusterID|
  total_patients = cluster.length
  cluster.each do |phenotypes|
    phenotypes.each do |p|
      phenotypes_frequency[p] += 1
    end
  end
  new_cluster_phenotypes[clusterID] = [total_patients, phenotypes_frequency.keys, phenotypes_frequency.values.map{|v| v.fdiv(total_patients) * 100}]
  phenotypes_frequency = Hash.new(0)
end

container = {
  :temp_folder => temp_folder,
  # :top_cluster_phenotypes => top_cluster_phenotypes.length,
  :summary_stats => summary_stats,
  :hpo_stats => hpo_stats,
  :all_cnvs_length => all_cnvs_length,
  :all_sor_length => all_sor_length,
  :new_cluster_phenotypes => new_cluster_phenotypes.keys.length,
  :ontology_levels => ontology_levels,
  :distribution_percentage => distribution_percentage,
 }
# top_cluster_phenotypes.each_with_index do |cluster, i|
#   clust_pr = cluster.map{|pr| [pr.join(', ')] }
#   container["clust_#{i}"] = clust_pr
# end

clust_info = []
new_cluster_phenotypes.each do |clusterID, info|
    phens = info[1].join(', ')
    freqs = info[2].map{|a| a.round(4)}.join(', ')
    clust_info << [info[0], phens, freqs]
    container["clust_#{clusterID}"] = clust_info
    clust_info = []
end

template = File.open(File.join(REPORT_FOLDER, 'cohort_report.erb')).read
report = Report_html.new(container, 'Cohort quality report')
report.build(template)
report.write(options[:output_file]+'.html')