require 'genomic_features'
class Reference_parser

	def self.load(file_path, file_format: nil, feature_type: nil)
		file_format = file_path.split('.', 2).last if file_format.nil?
		if file_format == 'gtf'
			regions, all_attrs = parse_gtf(file_path, feature_type: feature_type)
		end

		return Genomic_Feature.new(regions, annotations: all_attrs)
	end

	def self.parse_gtf(file_path, feature_type: nil) # https://www.ensembl.org/info/website/upload/gff.html
		features = []
		all_attrs = {}
		File.open(file_path).each do |line|
			next if /^#/ =~ line
			seqname, source, feature, start, stop, score, strand, frame, attribute = line.chomp.split("\t")
			if feature_type.nil? || feature_type == feature
				attrs = process_attrs(attribute, ';', ' ')
				attrs['source'] = source
				attrs['feature'] = feature
				id = attrs['gene_id']
				features << [seqname.gsub('chr',''), start.to_i, stop.to_i, id]
				all_attrs[id] = attrs
			end
		end
		return features, all_attrs
	end

	private
	def self.process_attrs(attributes, tuple_sep, field_sep)
		return attributes.split(tuple_sep).map{|attr_pair| 
			tuple = attr_pair.strip.split(field_sep, 2)
			tuple.last.gsub!('"','')
			tuple
		}.to_h
	end
end
