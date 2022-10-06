options[:chromosome_col] = nil
opts.on("-c", "--chromosome_col INTEGER/STRING", "Column name if header is true, otherwise 0-based position of the column with the chromosome") do |data|
  options[:chromosome_col] = data
end

options[:id_col] = nil
opts.on("-d", "--pat_id_col INTEGER/STRING", "Column name if header is true, otherwise 0-based position of the column with the patient id") do |data|
  options[:id_col] = data
end

options[:end_col] = nil
opts.on("-e", "--end_col INTEGER/STRING", "Column name if header is true, otherwise 0-based position of the column with the end mutation coordinate") do |data|
  options[:end_col] = data
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

options[:sex_col] = nil
opts.on("-x", "--sex_col INTEGER/STRING", "Column name if header is true, otherwise 0-based position of the column with the patient sex") do |data|
  options[:sex_col] = data
end
