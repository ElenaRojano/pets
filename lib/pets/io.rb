require 'csv'

def load_hpo_ontology(hpo_file, excluded_hpo_file)
  hpo = nil
  if !hpo_file.include?('.json')
    if !excluded_hpo_file.nil?
      hpo = Ontology.new(file: hpo_file, load_file: true, removable_terms: read_excluded_hpo_file(excluded_hpo_file))
    else
      hpo = Ontology.new(file: hpo_file, load_file: true)
    end
  else
    hpo = Ontology.new
    hpo.read(hpo_file)
    if !excluded_hpo_file.nil?
      hpo.add_removable_terms(read_excluded_hpo_file(excluded_hpo_file))
      hpo.remove_removable()
      hpo.build_index()
    end
  end
  return hpo
end

def read_excluded_hpo_file(file)
  excluded_hpo = []
  File.open(file).each do |line|
    excluded_hpo << line.chomp
  end
  return excluded_hpo
end

def write_matrix_for_R(matrix, x_names, y_names, file)
  File.open(file, 'w') do |f|
    f.puts x_names.join("\t")
    matrix.each_with_index do |row, i|
      f.puts [y_names[i]].concat(row).join("\t")
    end
  end
end

def write_cluster_ic_data(all_ics, profile_lengths, cluster_ic_data_file, limit)
  File.open(cluster_ic_data_file, 'w') do |f|
    f.puts %w[cluster_id ic Plen].join("\t")
    all_ics.each_with_index do |cluster_ics, i|
      break if i == limit
      cluster_length = cluster_ics.length
      cluster_ics.each_with_index do |clust_ic, j|
        f.puts "#{cluster_length}_#{i}\t#{clust_ic}\t#{profile_lengths[i][j]}"
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

def write_arrays4scatterplot(x_axis_value, y_axis_value, filename, x_axis_name, y_axis_name)
  File.open(filename, 'w') do |f|
    f.puts "#{x_axis_name}\t#{y_axis_name}"
    x_axis_value.each_with_index do |value,i|
      y_value = y_axis_value[i]
      raise("The #{i} position is not presented in y_axis_value") if y_value.nil?
      f.puts [value, y_value].join("\t")
    end
  end
end


def write_similarity_matrix(similarity_matrix, similarity_matrix_file)  
  File.open(similarity_matrix_file, 'w') do |f|
    similarity_matrix.each do |row|
      f.puts row.join("\t")
    end
  end
end

def write_profile_pairs(similarity_pairs, filename)
  File.open(filename, 'w') do |f|
    similarity_pairs.each do |pairsA, pairsB_and_values|
      pairsB_and_values.each do |pairsB, values|
        f.puts "#{pairsA}\t#{pairsB}\t#{values}"
      end
    end
  end
end

def write_patient_hpo_stat(average_hp_per_pat_distribution, output_file)
  File.open(output_file, 'w') do |f|
    f.puts "#{'PatientsNumber'}\t#{'HPOAverage'}"
    average_hp_per_pat_distribution.each do |patient_num, ave|
      f.puts "#{patient_num}\t#{ave}"
    end
  end
end

def parse_clusters_file(clusters_file, patient_data)
  clusters_info = {}
  clusters_table = []
  File.open(clusters_file).each do |line|
    line.chomp!
    patientID, clusterID = line.split("\t")
    patientHPOProfile = patient_data.get_profile(patientID)
    query = clusters_info[clusterID]
    if query.nil? 
      clusters_info[clusterID] = {patientID => patientHPOProfile}
    else
      query[patientID] = patientHPOProfile
    end
  end
  clusters_info.each do |clusterID, patients_info|
    patients_per_cluster = patients_info.keys.length
    clusters_table << [clusterID, patients_per_cluster, patients_info.keys, patients_info.values]
  end
  return clusters_table, clusters_info
end

def load_profiles(file_path, hpo)
  profiles = {}
  #count = 0
  File.open(file_path).each do |line|
    id, profile = line.chomp.split("\t")
    hpos = profile.split(',').map{|a| a.to_sym}
    hpos, rejected_hpos = hpo.check_ids(hpos)
    if !hpos.empty?
      hpos = hpo.clean_profile(hpos)
      profiles[id] = hpos if !hpos.empty?
    end
  end
  return profiles
end

def load_variants(variant_folder)
  variants = {}
  coordinates = {}
  count = 0
  all_vars = {}
  Dir.glob(File.join(variant_folder, '*.tab')).each do |path|
    profile_id = File.basename(path, '.tab')
    vars = {}
    File.open(path).each do |line|
      fields = line.chomp.split("\t")
      chr = fields[0]
      start = fields[1].to_i
      query = coordinates[chr]
      if query.nil?
        coordinates[chr] = [start]
        count += 1
        id = "var_#{count}"
      else
        if !query.include?(start)
          query << start
          count += 1
          id = "var_#{count}"
        else
          id = all_vars.key([chr, start]) 
        end
      end
      vars[id] = [chr, start]
    end
    all_vars.merge!(vars)
    variants[profile_id] = vars
  end
  return variants
end

def load_evidences(evidences_path, hpo)
  genomic_coordinates = {}
  coord_files = Dir.glob(File.join(evidences_path, '*.coords'))
  coord_files.each do |cd_f|
    entity = File.basename(cd_f, '.coords')
    coordinates = load_coordinates(cd_f)
    genomic_coordinates[entity] = coordinates
  end
  evidences = {}
  evidence_files = Dir.glob(File.join(evidences_path, '*_HP.txt'))
  evidence_files.each do |e_f|
    pair = File.basename(e_f, '.txt')
    profiles, id2label = load_evidence_profiles(e_f, hpo)
    evidences[pair] = {prof: profiles, id2lab: id2label}
  end
  return evidences, genomic_coordinates
end

def load_coordinates(file_path)
  coordinates = {}
  header = true
  File.open(file_path).each do |line|
    fields = line.chomp.split("\t")
    if header
      header = false
    else
      entity, chr, strand, start, stop = fields
      coordinates[entity] = [chr, start.to_i, stop.to_i, strand]
    end
  end
  return coordinates
end

def load_evidence_profiles(file_path, hpo)
  profiles = {}
  id2label = {}
  #count = 0
  File.open(file_path).each do |line|
    id, label, profile = line.chomp.split("\t")
    hpos = profile.split(',').map{|a| a.to_sym}
    hpos, rejected_hpos = hpo.check_ids(hpos)
    if !hpos.empty?
      hpos = hpo.clean_profile(hpos)
      profiles[id] = hpos if !hpos.empty?
      id2label[id] = label
    end
  end
  return profiles, id2label
end

#Common methods for predictors
#Training file example = 9  131371492   131375954   HP:0010974  2.41161970596 9.3.A.5
#1. Indexing by chr (region)
def coor_overlap?(ref_start, ref_stop, start, stop)
  overlap = false
  if (stop > ref_start && stop <= ref_stop) ||
    (start >= ref_start && start < ref_stop) ||
    (start <= ref_start && stop >= ref_stop) ||
    (start > ref_start && stop < ref_stop)
    overlap = true
  end
  return overlap
end

def load_training_file4regions(training_file)
  training_set = {}
  posInfo = loadFile(training_file)
  posInfo.each do |info|
    chr = info.shift
    query = training_set[chr]
    if query.nil?
      training_set[chr] = [info]
    else
      query << info
    end
  end
  return training_set
end

#2. Indexing by hpo (code)
#prepare training file for analysis using phenotype2region prediction
def load_training_file4HPO(training_file, thresold=0)
  training_set = {}
  information = loadFile(training_file, thresold)
  information.each do |info|
    hpoCode = info.delete_at(4)
    query = training_set[hpoCode]
    if query.nil?
      training_set[hpoCode] = [info]  
    else
      query << info
    end
  end
  # STDERR.puts training_set.keys.inspect
  return training_set
end


#3. Load training info file:
#Chr;Start;Stop;HPO;Association;node
def loadFile(file, thresold=0)
  information = []
  File.open(file).each do |line|
    line.chomp!
    allInfo = line.split("\t")
    associationValue = allInfo[4].to_f
    if associationValue >= thresold
      chr = allInfo[0]
      startPos = allInfo[1].to_i
      stopPos = allInfo[2].to_i
      hpoCode = allInfo[3]
      nodeID = allInfo[5]
      information << [chr, startPos, stopPos, nodeID, hpoCode, associationValue]
    end
  end
  return information
end

def load_hpo_ci_values(information_coefficient_file)
  hpos_ci_values = {}
  File.open(information_coefficient_file).each do |line|
    line.chomp!
    hpo_code, ci = line.split("\t")
    hpos_ci_values[hpo_code.to_sym] = ci.to_f
  end
  return hpos_ci_values
end

def load_clustered_patients(file)
  clusters = {}
  File.open(file).each do |line|
    line.chomp!
    pat_id, cluster_id = line.split("\t")
    query = clusters[cluster_id]
    if query.nil?
      clusters[cluster_id] = [pat_id]
    else
      query << pat_id
    end
  end
  return clusters
end

def load_gene_data(gene_data_path)
  gene_list = {} #geneID => attr
  gene_location = {} # chr => gene
  infile = open(gene_data_path)
  gz = Zlib::GzipReader.new(infile)
  current_chr = nil
  genes = []
  gz.each_line do |line|
    line.chomp!
    next if line =~ /^#/ 
    fields = line.split("\t")
    if fields[8].include?('genome=chromosome')
      chr = fields[8].split(';')[1].split('=').last
      gene_location[current_chr] = genes
      genes = []
      current_chr = chr
    elsif fields[2] == 'gene'
      attributes = {}
      fields[8].split(';').each do |pair| 
        key, value = pair.split('=')
        attributes[key] = value
      end
      geneName = nil
      geneName = attributes['gene'] if !attributes['gene'].nil?
      geneSyns = []
      geneSyns = attributes['gene_synonym'].split(',') if !attributes['gene_synonym'].nil?
      description = attributes['description']
      description = URI.unescape(description) if !description.nil?
      attributes['Dbxref'] =~ /GeneID:(\d+)/
      gene_list[$1] = [geneName, geneSyns, description]
      genes << [$1, fields[3].to_i, fields[4].to_i]
    end
  end
  gene_location[current_chr] = genes
  return gene_list, gene_location
end

def parse_kegg_data(query_genes)
  kegg_data = {} #gene => attb
    while !query_genes.empty?
      gene_set = query_genes.shift(10)
      url = "http://rest.kegg.jp/get/#{gene_set.map{|qg| "hsa:#{qg}"}.join('+')}"
      uri = URI(url)
      response = Net::HTTP.get(uri) 
    geneID = nil
    gene_names = []
    definition = nil
    pathways = []
    parsing_pathway_field = false
    response.squeeze(' ').each_line do |line|
      line.chomp!
      if line =~ /^ENTRY/
        geneID = line.split(' ')[1]
      elsif line =~ /^NAME/
        gene_names = line.split(' ', 2).last.split(', ')
      elsif line =~ /^DEFINITION/
        definition = line.split(' ', 2)[1]
      elsif line =~ /^PATHWAY/
        pathways << line.split(' ', 3)[1..2]
        parsing_pathway_field = true
      elsif line =~ /^BRITE/ || line =~ /^POSITION/ || line =~ /^DISEASE/ || line =~ /^MODULE/ || line =~ /^DRUG_TARGET/ || line =~ /^NETWORK/
        parsing_pathway_field = false
      elsif parsing_pathway_field
        pathways << line.strip.split(' ', 2)
      elsif line == '///'
        parsing_pathway_field = false
        kegg_data[geneID] = [gene_names, definition, pathways]
        pathways = []
        gene_names = []
      end
    end
  end
  return kegg_data
end

def write_compressed_plain_file(data, path)
  File.open(path, 'w') do |f|
    gz = Zlib::GzipWriter.new(f)
    gz.write data.to_json
    gz.close
  end
end

def read_compressed_json(path)
  infile = open(path)
  gz = Zlib::GzipReader.new(infile)
  object = JSON.parse(gz.read)
  return object
end

def download(ftp_server, path, name)
  ftp = Net::FTP.new()
  ftp.connect(ftp_server)
  ftp.login
  ftp.getbinaryfile(path, name)
  ftp.close
end
