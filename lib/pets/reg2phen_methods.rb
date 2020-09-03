require 'report_html'

def predict_patient(predictions, training_set, threshold, transform, genes, genes_dictionary)
  results = {}
  predictions.each do |info|
    if genes
      geneID = info.shift
      info = genes_dictionary[geneID]
    end 
    chr, pt_start, pt_stop, is_GeneSynom = info
    sors = training_set[chr]
    next if sors.nil?
    sors.each do |sor_start, sor_stop, nodeID, hpo_code, association_score|
      next if !coor_overlap?(pt_start, pt_stop, sor_start, sor_stop) 
      patientsNum = nodeID.split('.').last
      if association_score >= threshold 
        association_score = 10**(-association_score) if transform
        genes_in_sor = get_genes_by_coordinates(genes_dictionary, chr, sor_start, sor_stop)
        record2save = [chr, pt_start, pt_stop, hpo_code, association_score, sor_start, sor_stop, patientsNum, genes_in_sor.join(',')]
        key = info[0...3].join(":")
        key.concat(" (#{geneID})") if genes
        query_res = results[key]
        if query_res.nil?
            results[key] = [record2save]
        else
            query_res << record2save
        end
      end
    end
  end
  return results
end

def get_genes_by_coordinates(genes_dictionary, chr, start, stop)
  genes_in_sor = genes_dictionary.select do |gene_sym, attributes|
    gene_chr, gene_start, gene_stop, is_GeneSynom = attributes
    if gene_chr == chr && !is_GeneSynom
      coor_overlap?(gene_start, gene_stop, start, stop)
    else
      false
    end
  end
  return genes_in_sor.keys
end

def translate_hpos_in_results(results, hpo)
  results.each do |coords, data|
    data.each do |info|
      # hpo_name, rejected = hpo.translate_codes2names([info[3]])
      info[3] = hpo.translate_id(info[3])
      # info[3] = hpo_name.first
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
  reestructured_gene_locations = generate_gene_locations(gene_location)
  genes_dictionary = {}
  genes_with_kegg.each do |geneID, annotInfo|
    #STDERR.puts annotInfo.shift.inspect
    gene_location_data = reestructured_gene_locations[geneID]
    unless gene_location_data.nil?
      geneName, description, pathways, geneSyns = annotInfo
      genes_dictionary[geneName] = gene_location_data.dup.concat([false])
      geneSyns.each do |gene_symbol|
        genes_dictionary[gene_symbol] = gene_location_data.dup.concat([true])
      end
    end
  end
  return genes_dictionary
end

def parse_patient_results(results)
  patient_results = []
  results.each do |k, val|
    val.each do |v|
      patient_results << v.unshift(k)
    end
  end
  return patient_results
end

def report_data(container, html_file)
  template = File.open(File.join(REPORT_FOLDER, 'reg2phen_report.erb')).read
  report = Report_html.new(container, 'Patient summary')
  report.build(template)
  report.write(html_file)
end
