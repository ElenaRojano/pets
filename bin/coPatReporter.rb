#! /usr/bin/env ruby

ROOT_PATH = File.dirname(__FILE__)
REPORT_FOLDER = File.expand_path(File.join(ROOT_PATH, '..', 'templates'))
EXTERNAL_DATA = File.expand_path(File.join(ROOT_PATH, '..', 'external_data'))
EXTERNAL_CODE = File.expand_path(File.join(ROOT_PATH, '..', 'external_code'))
HPO_FILE = File.join(EXTERNAL_DATA, 'hp.obo')
IC_FILE = File.join(EXTERNAL_DATA, 'uniq_hpo_with_CI.txt')
CHR_SIZE = File.join(EXTERNAL_DATA, 'chromosome_sizes_hg19.txt')
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'pets'))

require 'optparse'
require 'csv'
require 'generalMethods.rb'
require 'coPatReporterMethods.rb'
require 'report_html'
require 'semtools'

#Expand class (semtools modifications if necessary):
class Ontology

end

##########################
# FUNCTIONS
##########################
def read_excluded_hpo_file(file)
  excluded_hpo = []
  File.open(file).each do |line|
    excluded_hpo << line.chomp
  end
  return excluded_hpo
end

# def translate_codes(clusters, hpo)
#   clusters.each do |clusterID, num_of_pats, patientIDs_ary, patient_hpos_ary|
#     patient_hpos_ary.map!{|patient_hpos| patient_hpos.map{|hpo_code| hpo.translate_id(hpo_code)}}
#   end
#   return clusters
# end


def translate_codes(clusters, hpo)
  translated_clusters = []
  clusters.each do |clusterID, num_of_pats, patientIDs_ary, patient_hpos_ary|
    #patient_hpos_ary.each do |patient_hpos| 
      #patient_hpos.each do |hpo_code| 
        translate_codes = patient_hpos_ary.map{|patient_hpos| patient_hpos.map{|hpo_code| hpo.translate_id(hpo_code)}}
        translated_clusters << [clusterID, 
          num_of_pats, 
          patientIDs_ary, 
          patient_hpos_ary, 
          translate_codes
        ]
      #end
    #end
  end
  #STDERR.puts translated_clusters.inspect
  #Process.exit
  return translated_clusters
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

  options[:ic_stats] = 'freq'
  opts.on("-t", "--ic_stats STRING", "'freq' to compute IC based en hpo frequency in the input cohort. 'freq_internal' to use precomputed internal IC stats. 'onto' to compute ic based on hpo ontology structure.. Default freq") do |data|
    options[:ic_stats] = data
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
hpo_profile_ic_file = File.join(temp_folder, 'hpo_ic_profile.txt')
hpo_frequency_file = File.join(temp_folder, 'hpo_cohort_frequency.txt')
parents_per_term_file = File.join(temp_folder, 'parents_per_term.txt')
clustered_patients_file = File.join(temp_folder, 'cluster_asignation')
cluster_ic_data_file = File.join(temp_folder, 'cluster_ic_data.txt')
cluster_chromosome_data_file = File.join(temp_folder, 'cluster_chromosome_data.txt')
coverage_to_plot_file = File.join(temp_folder, 'coverage_data.txt')
sor_coverage_to_plot_file = File.join(temp_folder, 'sor_coverage_data.txt')
similarity_matrix_resnik_file = File.join(temp_folder, 'similarity_matrix_resnik.txt')
similarity_matrix_lin_file = File.join(temp_folder, 'similarity_matrix_lin.txt')
similarity_matrix_jiang_file = File.join(temp_folder, 'similarity_matrix_jiang.txt')
profiles_similarity_resnik_file = File.join(temp_folder, 'profiles_similarity_resnik.txt')
jiang_clusters_distribution_file = File.join(temp_folder, 'jiang_clusters_distribution.txt')
resnik_clusters_distribution_file = File.join(temp_folder, 'resnik_clusters_distribution.txt')
lin_clusters_distribution_file = File.join(temp_folder, 'lin_clusters_distribution.txt')

# cnvs_lenght_to_plot_file = File.join(temp_folder, 'cnvs_lenght.txt')
Dir.mkdir(temp_folder) if !File.exists?(temp_folder)

hpo_file = ENV['hpo_file']
hpo_file = HPO_FILE if hpo_file.nil?

#############
# >> UPDATED TO OBO_HANDLER
# # add_removable_terms_from_file
# hpo = Ontology.new
# hpo.load_black_list(options[:excluded_hpo]) if !options[:excluded_hpo].nil?
# hpo.load_data(hpo_file)
if !options[:excluded_hpo].nil?
  hpo = Ontology.new(file: hpo_file, load_file: true, removable_terms: read_excluded_hpo_file(options[:excluded_hpo]))
else
  hpo = Ontology.new(file: hpo_file, load_file: true)
end

#hash_keys
patient_data = load_patient_cohort(options)

cohort_hpos, suggested_childs, rejected_hpos, fraction_terms_specific_childs = format_patient_data(patient_data, options, hpo)
patient_uniq_profiles = get_uniq_hpo_profiles(patient_data)
hpo.load_profiles(patient_uniq_profiles)

profiles_similarity_resnik = hpo.compare_profiles
profiles_similarity_lin = hpo.compare_profiles(sim_type: :lin)
profiles_similarity_jiang = hpo.compare_profiles(sim_type: :jiang_conrath)

resnik_profile_pairs = write_profile_pairs(profiles_similarity_resnik, profiles_similarity_resnik_file)


similarity_matrix_resnik = format_profiles_similarity_data(profiles_similarity_resnik)
similarity_matrix_lin = format_profiles_similarity_data(profiles_similarity_lin)
similarity_matrix_jiang = format_profiles_similarity_data(profiles_similarity_jiang)


profile_sizes, parental_hpos_per_profile = get_profile_redundancy(hpo)
ontology_levels, distribution_percentage = get_profile_ontology_distribution_tables(hpo)

# onto_ic, freq_ic = hpo.get_ic_by_onto_and_freq(hpo_file)
onto_ic, freq_ic = hpo.get_observed_ics_by_onto_and_freq
observed_keys = onto_ic.keys
observed_keys.each do |k| # duplitacte to symbols
  # onto_ic[k.to_s] = onto_ic.delete(k)
  # freq_ic[k.to_s] = freq_ic.delete(k)
  onto_ic[k.to_s] = onto_ic[k]
  freq_ic[k.to_s] = freq_ic[k]
end
# onto_ic_profile, freq_ic_profile = hpo.get_ic_profile_by_onto_and_freq
onto_ic_profile, freq_ic_profile = hpo.get_profiles_resnik_dual_ICs
onto_ic_profile = onto_ic_profile.values
freq_ic_profile = freq_ic_profile.values
clustered_patients = cluster_patients(patient_data, cohort_hpos, matrix_file, clustered_patients_file) 
all_ics, cluster_data_by_chromosomes, top_cluster_phenotypes, multi_chromosome_patients = process_clustered_patients(options, clustered_patients, patient_data, hpo, onto_ic, freq_ic, options[:pat_id_col])
get_patient_hpo_frequency(patient_uniq_profiles, hpo_frequency_file)



summary_stats = get_summary_stats(patient_data, cohort_hpos, hpo)
summary_stats << ['Percentage of defined HPOs that have more specific childs', (fraction_terms_specific_childs * 100).round(4)]

# Move code 'Percentage of defined HPOs that have more specific childs' outside the next function
# hpo_stats = hpo.get_term_frequency_from_profiles(names=true)[0..20]
hpo_stats = hpo.get_profiles_terms_frequency()
hpo_stats.map{ |stat| stat[1] = stat[1]*100}
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
write_arrays4scatterplot(onto_ic.values, freq_ic.values, hpo_ic_file, 'OntoIC', 'FreqIC')
write_arrays4scatterplot(onto_ic_profile, freq_ic_profile, hpo_profile_ic_file, 'OntoIC', 'FreqIC')
write_arrays4scatterplot(profile_sizes, parental_hpos_per_profile, parents_per_term_file, 'ProfileSize', 'ParentTerms')
write_similarity_matrix(similarity_matrix_resnik, similarity_matrix_resnik_file)
write_similarity_matrix(similarity_matrix_lin, similarity_matrix_lin_file)
write_similarity_matrix(similarity_matrix_jiang, similarity_matrix_jiang_file)
system("#{File.join(EXTERNAL_CODE, 'plot_scatterplot_simple.R')} #{hpo_ic_file} #{File.join(temp_folder, 'hpo_ics.pdf')} 'OntoIC' 'FreqIC' 'HP Ontology IC' 'HP Frequency based IC'")
system("#{File.join(EXTERNAL_CODE, 'plot_scatterplot_simple.R')} #{hpo_profile_ic_file} #{File.join(temp_folder, 'hpo_profile_ics.pdf')} 'OntoIC' 'FreqIC' 'HP Ontology Profile IC' 'HP Frequency based Profile IC'")
system("#{File.join(EXTERNAL_CODE, 'plot_scatterplot_simple.R')} #{parents_per_term_file} #{File.join(temp_folder, 'parents_per_term.pdf')} 'ProfileSize' 'ParentTerms' 'Patient HPO profile size' 'Parent HPO terms within the profile'")

###Cohort frequency calculation
system("#{File.join(EXTERNAL_CODE, 'ronto_plotter.R')} -i #{hpo_frequency_file} -o #{File.join(temp_folder, 'hpo_freq_colour')} -O #{hpo_file}") 

###PLOTTING HEATMAPS
system("#{File.join(EXTERNAL_CODE, 'plot_heatmap.R')} -d #{similarity_matrix_resnik_file} -o #{File.join(temp_folder, 'resnik')} -m max -H") 
system("#{File.join(EXTERNAL_CODE, 'plot_heatmap.R')} -d #{similarity_matrix_lin_file} -o #{File.join(temp_folder, 'lin')} -m comp1 -H")    
system("#{File.join(EXTERNAL_CODE, 'plot_heatmap.R')} -d #{similarity_matrix_jiang_file} -o #{File.join(temp_folder, 'jiang')} -H")

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
# GENERAL COHORT ANALYZER REPORT
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
  :similarity_matrix_resnik => similarity_matrix_resnik,
  :similarity_matrix_lin => similarity_matrix_lin,
  :similarity_matrix_jiang => similarity_matrix_jiang
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

#----------------------------------
# CLUSTER COHORT ANALYZER REPORT
#----------------------------------

#jiang_clusters = parse_clusters_file(File.join(temp_folder, 'jiang_clusters.txt'), patient_uniq_profiles)
jiang_clusters_codes, jiang_clusters_info = parse_clusters_file(File.join(temp_folder, 'jiang_clusters.txt'), patient_uniq_profiles)
get_cluster_metadata(jiang_clusters_info, jiang_clusters_distribution_file)
system("#{File.join(EXTERNAL_CODE, 'xyplot_graph.R')} -d #{jiang_clusters_distribution_file} -o #{File.join(temp_folder, 'jiang_clusters_distribution')} -x PatientsNumber -y HPOAverage")
jiang_clusters = translate_codes(jiang_clusters_codes, hpo)


#resnik_clusters = parse_clusters_file(File.join(temp_folder, 'resnik_clusters.txt'), patient_uniq_profiles)
resnik_clusters_codes, resnik_clusters_info = parse_clusters_file(File.join(temp_folder, 'resnik_clusters.txt'), patient_uniq_profiles)
get_cluster_metadata(resnik_clusters_info, resnik_clusters_distribution_file)
system("#{File.join(EXTERNAL_CODE, 'xyplot_graph.R')} -d #{resnik_clusters_distribution_file} -o #{File.join(temp_folder, 'resnik_clusters_distribution')} -x PatientsNumber -y HPOAverage")
resnik_clusters = translate_codes(resnik_clusters_codes, hpo)

#lin_clusters = parse_clusters_file(File.join(temp_folder, 'lin_clusters.txt'), patient_uniq_profiles)
lin_clusters_codes, lin_clusters_info = parse_clusters_file(File.join(temp_folder, 'lin_clusters.txt'), patient_uniq_profiles)
get_cluster_metadata(lin_clusters_info, lin_clusters_distribution_file)
system("#{File.join(EXTERNAL_CODE, 'xyplot_graph.R')} -d #{lin_clusters_distribution_file} -o #{File.join(temp_folder, 'lin_clusters_distribution')} -x PatientsNumber -y HPOAverage")
lin_clusters = translate_codes(lin_clusters_codes, hpo)

container = {
  :temp_folder => temp_folder,
  :jiang_clusters => jiang_clusters,
  :resnik_clusters => resnik_clusters,
  :lin_clusters => lin_clusters,
  :hpo => hpo
 }

template = File.open(File.join(REPORT_FOLDER, 'cluster_report.erb')).read
report = Report_html.new(container, 'Patient clusters report')
report.build(template)
report.write(options[:output_file]+'_clusters.html')