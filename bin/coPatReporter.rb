#! /usr/bin/env ruby

ROOT_PATH = File.dirname(__FILE__)
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'pets'))

require 'benchmark'
require 'parallel'
require 'optparse'
require 'csv'
require 'npy'
require 'report_html'
require 'semtools'
require 'pets'

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

  options[:id_col] = nil
  opts.on("-d", "--pat_id_col INTEGER/STRING", "Column name if header is true, otherwise 0-based position of the column with the patient id") do |data|
    options[:id_col] = data
  end

  options[:excluded_hpo] = nil
  opts.on("-E", "--excluded_hpo PATH", "List of HPO phenotypes to exclude (low informative)") do |excluded_hpo|
    options[:excluded_hpo] = excluded_hpo
  end

  options[:end_col] = nil
  opts.on("-e", "--end_col INTEGER/STRING", "Column name if header is true, otherwise 0-based position of the column with the end mutation coordinate") do |data|
    options[:end_col] = data
  end

  options[:minClusterProportion] = 0.01
  opts.on("-M", "--minClusterProportion FLOAT", "Minimum percentage of patients per cluster") do |data|
    options[:minClusterProportion] = data.to_f
  end

  options[:patients_filter] = 2
  opts.on("-f", "--patients_filter INTEGER", "Minimum number of patients sharing SORs. Default 0") do |data|
    options[:patients_filter] = data.to_i
  end

  options[:clusters2graph] = 30
  opts.on("-g", "--clusters2graph INTEGER", "How may patient clusters are plotted in cluster plots. Default 30") do |data|
    options[:clusters2graph] = data.to_i
  end

  options[:genome_assembly] = 'hg38'
  opts.on("-G", "--genome_assembly STRING", "Genome assembly version. Please choose between hg18, hg19 and hg38. Default hg38") do |data|
    options[:genome_assembly] = data
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

  options[:clustering_methods] = ['resnik', 'jiang_conrath', 'lin']
  opts.on("-m", "--clustering_methods ARRAY", "Clustering methods") do |data|
    options[:clustering_methods] = data.split(',')
  end

  options[:names] = false
  opts.on("-n", "--hpo_names", "Define if the input HPO are human readable names. Default false") do
    options[:names] = true
  end

  options[:output_file] = nil
  opts.on("-o", "--output_file PATH", "Output file with patient data") do |data|
    options[:output_file] = data
  end

  options[:hpo_file] = nil
  opts.on("-P", "--hpo_file PATH", "Input HPO file for extracting HPO codes") do |value|
    options[:hpo_file] = value
  end

  options[:ont_col] = nil
  opts.on("-p", "--hpo_term_col INTEGER/STRING", "Column name if header true or 0-based position of the column with the HPO terms") do |data|
  	options[:ont_col] = data
  end

  options[:separator] = '|'
  opts.on("-S", "--hpo_separator STRING", "Set which character must be used to split the HPO profile. Default '|'") do |data|
  	options[:separator] = data
  end

  options[:start_col] = nil
  opts.on("-s", "--start_col INTEGER/STRING", "Column name if header is true, otherwise 0-based position of the column with the start mutation coordinate") do |data|
  	options[:start_col] = data
  end

  options[:root_node] = "HP:0000118"
  opts.on("-r", "--root_node", "Root node from which the ontology will be represented ") do |root_node|
    options[:root_node] = root_node
  end

  options[:ic_stats] = 'freq'
  opts.on("-t", "--ic_stats STRING", "'freq' to compute IC based en hpo frequency in the input cohort. 'freq_internal' to use precomputed internal IC stats. 'onto' to compute ic based on hpo ontology structure.. Default freq") do |data|
    options[:ic_stats] = data
  end

  options[:threads] = 1
  opts.on("-T", "--threads INTEGER", "Number of threads to be used in calculations. Default 1") do |data|
    options[:threads] = data.to_i
  end



  opts.on_tail("-h", "--help", "Show this message") do
    puts opts
    exit
  end

end.parse!


##########################
#MAIN
##########################

if options[:genome_assembly] == 'hg19' || options[:genome_assembly] == 'hg37'
  CHR_SIZE = File.join(EXTERNAL_DATA, 'chromosome_sizes_hg19.txt')
elsif options[:genome_assembly] == 'hg38'
  CHR_SIZE = File.join(EXTERNAL_DATA, 'chromosome_sizes_hg38.txt')
elsif options[:genome_assembly] == 'hg18'
  CHR_SIZE = File.join(EXTERNAL_DATA, 'chromosome_sizes_hg18.txt')
else
  abort('Wrong human genome assembly. Please choose between hg19, hg18 or hg38.')
end

output_folder = File.dirname(options[:output_file])
detailed_profile_evaluation_file = File.join(output_folder, 'detailed_hpo_profile_evaluation.csv')
rejected_file = File.join(output_folder, 'rejected_records.txt')
temp_folder = File.join(output_folder, 'temp')
matrix_file = File.join(temp_folder, 'pat_hpo_matrix.npy')
hpo_ic_file = File.join(temp_folder, 'hpo_ic.txt')
hpo_profile_ic_file = File.join(temp_folder, 'hpo_ic_profile.txt')
hpo_frequency_file = File.join(temp_folder, 'hpo_cohort_frequency.txt')
parents_per_term_file = File.join(temp_folder, 'parents_per_term.txt')
clustered_patients_file = File.join(temp_folder, 'cluster_asignation')
cluster_ic_data_file = File.join(temp_folder, 'cluster_ic_data.txt')
cluster_chromosome_data_file = File.join(temp_folder, 'cluster_chromosome_data.txt')
coverage_to_plot_file = File.join(temp_folder, 'coverage_data.txt')
sor_coverage_to_plot_file = File.join(temp_folder, 'sor_coverage_data.txt')
ronto_file = File.join(temp_folder, 'hpo_freq_colour')


Dir.mkdir(temp_folder) if !File.exists?(temp_folder)

hpo_file = !ENV['hpo_file'].nil? ? ENV['hpo_file'] : HPO_FILE
Cohort.load_ontology(:hpo, hpo_file, options[:excluded_hpo])
Cohort.act_ont = :hpo

patient_data, rejected_hpos_L, rejected_patients_L = Cohort_Parser.load(options)
rejected_hpos_C, rejected_patients_C = patient_data.check
rejected_hpos = rejected_hpos_L | rejected_hpos_C
rejected_patients = rejected_patients_L + rejected_patients_C
File.open(rejected_file, 'w'){|f| f.puts (rejected_patients).join("\n")}

patient_data.link2ont(Cohort.act_ont) # TODO: check if method load should call to this and use the semtools checking methods (take care to only remove invalid terms)

profile_sizes, parental_hpos_per_profile = patient_data.get_profile_redundancy
_, _ = patient_data.check(hard=true)
hpo_stats = patient_data.get_profiles_terms_frequency() # hpo NAME, freq
hpo_stats.each{ |stat| stat[1] = stat[1]*100}
File.open(hpo_frequency_file, 'w') do |f|
  patient_data.get_profiles_terms_frequency(translate: false).each do |hpo_code, freq| # hpo CODE, freq
    f.puts "#{hpo_code.to_s}\t#{freq}"
  end
end
suggested_childs, fraction_terms_specific_childs = patient_data.compute_term_list_and_childs()
ontology_levels, distribution_percentage = patient_data.get_profile_ontology_distribution_tables()
onto_ic, freq_ic, onto_ic_profile, freq_ic_profile = patient_data.get_ic_analysis()

if options[:ic_stats] == 'freq_internal'
  ic_file = !ENV['ic_file'].nil? ? ENV['ic_file'] : IC_FILE
  freq_ic = load_hpo_ci_values(ic_file)
  phenotype_ic = freq_ic
  freq_ic_profile = {}
  patient_data.each_profile do |pat_id, phenotypes|
    freq_ic_profile[pat_id] = get_profile_ic(phenotypes, phenotype_ic)
  end
elsif options[:ic_stats] == 'freq'
  phenotype_ic = freq_ic
elsif options[:ic_stats] == 'onto'
  phenotype_ic = onto_ic
end

clustered_patients = dummy_cluster_patients(patient_data.profiles, matrix_file, clustered_patients_file) 
all_ics, prof_lengths, clust_by_chr, top_clust_phen, multi_chr_pats = process_dummy_clustered_patients(options, clustered_patients, patient_data, phenotype_ic)

summary_stats = get_summary_stats(patient_data, rejected_patients, hpo_stats, fraction_terms_specific_childs, rejected_hpos)

all_cnvs_length = []
if !options[:chromosome_col].nil?
  summary_stats << ['Number of clusters with mutations accross > 1 chromosomes', multi_chr_pats]
  
  #----------------------------------
  # Prepare data to plot coverage
  #----------------------------------
  if options[:coverage_analysis]
    patient_data.index_vars
    all_cnvs_length = patient_data.get_vars_sizes(true)
    cnv_size_average = get_mean_size(all_cnvs_length)
    patients_by_cluster, sors = patient_data.generate_cluster_regions(:reg_overlap, 'A', 0)

    ###1. Process CNVs
    raw_coverage, n_cnv, nt, pats_per_region = calculate_coverage(sors)
    summary_stats << ['Average variant size', cnv_size_average.round(4)]
    summary_stats << ['Nucleotides affected by mutations', nt]
    summary_stats << ['Number of genome windows', n_cnv]
    summary_stats << ['Mean patients per genome window', pats_per_region.round(4)]
    coverage_to_plot = get_final_coverage(raw_coverage, options[:bin_size])

    ###2. Process SORs
    raw_sor_coverage, n_sor, nt, pats_per_region = calculate_coverage(sors, options[:patients_filter] - 1)
    summary_stats << ["Number of genome window shared by >= #{options[:patients_filter]} patients", n_sor]
    summary_stats << ["Number of patients with at least 1 SOR", patients_by_cluster.length]
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
write_arrays4scatterplot(onto_ic.values, freq_ic.values, hpo_ic_file, 'OntoIC', 'FreqIC') # hP terms
write_arrays4scatterplot(onto_ic_profile.values, freq_ic_profile.values, hpo_profile_ic_file, 'OntoIC', 'FreqIC') #HP profiles
write_arrays4scatterplot(profile_sizes, parental_hpos_per_profile, parents_per_term_file, 'ProfileSize', 'ParentTerms')
write_cluster_ic_data(all_ics, prof_lengths, cluster_ic_data_file, options[:clusters2graph])

system_call(EXTERNAL_CODE, 'plot_scatterplot_simple.R', "-i #{hpo_ic_file} -o #{File.join(temp_folder, 'hpo_ics.pdf')} -x 'OntoIC' -y 'FreqIC' --x_tag 'HP Ontology IC' --y_tag 'HP Frequency based IC' --x_lim '0,4.5' --y_lim '0,4.5'") if !File.exists?(File.join(temp_folder, 'hpo_ics.pdf'))
system_call(EXTERNAL_CODE, 'plot_scatterplot_simple.R', "-i #{hpo_profile_ic_file} -o #{File.join(temp_folder, 'hpo_profile_ics.pdf')} -x 'OntoIC' -y 'FreqIC' --x_tag 'HP Ontology Profile IC' --y_tag 'HP Frequency based Profile IC' --x_lim '0,4.5' --y_lim '0,4.5'") if !File.exists?(File.join(temp_folder, 'hpo_profile_ics.pdf'))
system_call(EXTERNAL_CODE, 'plot_scatterplot_simple.R', "-i #{parents_per_term_file} -o #{File.join(temp_folder, 'parents_per_term.pdf')} -x 'ProfileSize' -y 'ParentTerms' --x_tag 'Patient HPO profile size' --y_tag 'Parent HPO terms within the profile'")
system_call(EXTERNAL_CODE, 'ronto_plotter.R', "-i #{hpo_frequency_file} -o #{ronto_file} --root_node #{options[:root_node]} -O #{hpo_file.gsub('.json','.obo')}") if !File.exist?(ronto_file + '.png') ###Cohort frequency calculation
system_call(EXTERNAL_CODE, 'plot_boxplot.R', "#{cluster_ic_data_file} #{temp_folder} cluster_id ic 'Cluster size/id' 'Information coefficient' 'Plen' 'Profile size'")

if !options[:chromosome_col].nil?
  write_cluster_chromosome_data(clust_by_chr, cluster_chromosome_data_file, options[:clusters2graph])
  system_call(EXTERNAL_CODE, 'plot_scatterplot.R', "#{cluster_chromosome_data_file} #{temp_folder} cluster_id chr count 'Cluster size/id' 'Chromosome' 'Patients'")
  if options[:coverage_analysis]
    ###1. Process CNVs
    write_coverage_data(coverage_to_plot, coverage_to_plot_file)
    system_call(EXTERNAL_CODE, 'plot_area.R', "-d #{coverage_to_plot_file} -o #{temp_folder}/coverage_plot -x V2 -y V3 -f V1 -H -m #{CHR_SIZE} -t CNV")
    ###2. Process SORs
    write_coverage_data(sor_coverage_to_plot, sor_coverage_to_plot_file)
    system_call(EXTERNAL_CODE, 'plot_area.R', "-d #{sor_coverage_to_plot_file} -o #{temp_folder}/sor_coverage_plot -x V2 -y V3 -f V1 -H -m #{CHR_SIZE} -t SOR")
  end
end

#----------------------------------
# GENERAL COHORT ANALYZER REPORT
#----------------------------------
new_cluster_phenotypes = get_top_dummy_clusters_stats(top_clust_phen)

container = {
  :temp_folder => temp_folder,
  # :top_clust_phen => top_clust_phen.length,
  :summary_stats => summary_stats,
  :clustering_methods => options[:clustering_methods],
  :hpo_stats => hpo_stats,
  :all_cnvs_length => all_cnvs_length,
  :all_sor_length => all_sor_length,
  :new_cluster_phenotypes => new_cluster_phenotypes.keys.length,
  :ontology_levels => ontology_levels,
  :distribution_percentage => distribution_percentage
}

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
get_semantic_similarity_clustering(options, patient_data, temp_folder)