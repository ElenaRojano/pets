#! /usr/bin/env ruby
#DECIPHER predictor system, using data from cross validation
#data2predict = file to predict
#training_file.txt = file with training data (association values and hpo codes).

REPORT_FOLDER=File.expand_path(File.join(File.dirname(__FILE__), '..', 'templates'))
ROOT_PATH = File.dirname(__FILE__)
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'pets'))
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'pets', 'simnetome', 'lib'))
require 'optparse'
require 'generalMethods.rb'
require 'ontology'

##########################
#METHODS
##########################
def predict_patient(predictions, training_set, threshold, transform, genes, genes_dictionary)
  results = {}
  predictions.each do |info|
    info = genes_dictionary[info.shift] if genes
    chr, pt_start, pt_stop = info
    query = training_set[chr]
    next if query.nil?
    query.each do |hpo_start, hpo_stop, nodeID, hpo_code, association_score|
      next if !coor_overlap?(pt_start, pt_stop, hpo_start, hpo_stop) 
      if association_score >= threshold 
        association_score = 10**(-association_score) if transform
        record2save = [chr, pt_start, pt_stop, hpo_code, association_score, hpo_start, hpo_stop]
        key = info.join(":")
        save = results[key]
        if save.nil?
            results[key] = [record2save]
        else
            save << record2save
        end
      end
    end
  end
  return results
end

def translate_hpos_in_results(results, hpo)
  results.each do |coords, data|
    data.each do |info|
      hpo_name, rejected = hpo.translate_codes2names([info[3]])
      info[3] = hpo_name.first
    end
  end
end

def generate_gene_locations(gene_location)
  gene_locations = {}
  gene_location.each do |chr, genesInfo|
    genesInfo.each do |geneID, geneStart, geneStop|
      gene_locations[geneID] = [chr, geneStart, geneStop]
    end
  end
  return gene_locations
end

def generate_genes_dictionary(gene_location, genes_with_kegg)
  gene_locations = generate_gene_locations(gene_location)
  genes_dictionary = {}
  genes_with_kegg.each do |geneID, annotInfo|
    #STDERR.puts annotInfo.shift.inspect
    gene_codes = annotInfo.shift
    unless gene_codes.empty?
      gene_codes.each do |gene_code|
        genes_dictionary[gene_code] = gene_locations[geneID]
      end
    end
  end
  return genes_dictionary
end

##########################
#OPT-PARSER
##########################

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"

  options[:hpo_file] = nil
  opts.on("-b", "--hpo_file PATH", "Input HPO obo file") do |hpo_file|
    options[:hpo_file] = hpo_file
  end

  options[:association_limit] = 0
  opts.on("-l", "--association_limit FLOAT", "Threshold for association values") do |association_limit|
    options[:association_limit] = association_limit.to_f
  end

  options[:input_genes] = false
    opts.on("-g", "--input_genes", "Input list of genes instead of chromosomic coordinates") do
  options[:input_genes] = true
  end

  options[:output_path] = "output.txt"
  opts.on("-o", '--output_path PATH', 'Output path for overlapping patient file') do |output_path|
    options[:output_path] = output_path
  end

  options[:prediction_file] = nil
  #chr\tstart\tstop
  opts.on("-p", "--prediction_file PATH", "Input file with regions to predict") do |input_path|
    options[:prediction_file] = input_path
  end

  options[:transform_pvalues] = false
    opts.on("-P", "--transform_pvalues", "Set to transform association values into P-values (only for HyI associations)") do
  options[:transform_pvalues] = true
  end

  options[:top_results] = 10
  opts.on("-r", "--top_results FLOAT", "Threshold for association values") do |top_results|
    options[:top_results] = top_results.to_i
  end

  options[:training_file] = nil
  #chr\tstart\tstop\tphenotype\tassociation_value
  opts.on("-t", "--training_file PATH", "Input training file, with association values") do |training_path|
    options[:training_file] = training_path
  end

  opts.on_tail("-h", "--help", "Show this message") do
    puts opts
    exit
  end

end.parse!

##########################
#PATHS
##########################
all_paths = {code: File.join(File.dirname(__FILE__), '..')}
all_paths[:external_data] = File.join(all_paths[:code], 'external_data')
all_paths[:gene_data] = File.join(all_paths[:external_data], 'gene_data.gz')
all_paths[:biosystems_gene] = File.join(all_paths[:external_data], 'biosystems_gene.gz')
all_paths[:biosystems_info] = File.join(all_paths[:external_data], 'bsid2info.gz')
all_paths[:gene_data_with_pathways] = File.join(all_paths[:external_data], 'gene_data_with_pathways.gz')
all_paths[:gene_location] = File.join(all_paths[:external_data], 'gene_location.gz')


##########################
#MAIN
##########################
gene_location, genes_with_kegg = get_and_parse_external_data(all_paths)
hpo = Ontology.new
hpo.load_data(options[:hpo_file])
training_set = load_training_file4regions(options[:training_file])

genes_dictionary = {}
genes_dictionary = generate_genes_dictionary(gene_location, genes_with_kegg) if options[:input_genes]

multiple_regions = []
File.open(options[:prediction_file]).each do |line|
  line.chomp!
  fields = line.split("\t")
  if options[:input_genes]
    fields[2] = fields[2].to_i
    fields[3] = fields[3].to_i
  else
    fields[1] = fields[1].to_i
    fields[2] = fields[2].to_i
  end
  multiple_regions << fields
end

results = predict_patient(
  multiple_regions, 
  training_set, 
  options[:association_limit], 
  options[:transform_pvalues], 
  options[:input_genes], 
  genes_dictionary
)
translate_hpos_in_results(results, hpo)

results.each do |pred, values|
  values.sort! do |a, b|
    if !options[:transform_pvalues]
     b[4] <=> a[4]
   else
     a[4] <=> b[4] 
   end
  end
end

File.open(options[:output_path], 'w') do |f|
  results.each do |k, v|
    f.puts "Results for #{k}:"
    v.each_with_index do |i, c|      
        f.puts i.join("\t") if c < options[:top_results]
    end
  end
end