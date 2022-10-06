class Cohort_Parser
	def self.load(options)
		fields2extract = get_fields2extract(options)
		field_numbers = fields2extract.values
		records = read_records(options, fields2extract, field_numbers)
		options[:extracted_fields] = fields2extract.keys
		cohort, rejected_terms, rejected_recs = create_cohort(records, options)
		return cohort, rejected_terms, rejected_recs
	end 

	def self.read_records(options, fields2extract, field_numbers)
		records = {}
		count = 0
		File.open(options[:input_file]).each do |line|
			line.chomp!
			if options[:header] && count == 0
				line.gsub!(/#\s*/,'') # correct comment like	headers
				field_names = line.split("\t")
				get_field_numbers2extract(field_names, fields2extract)
				field_numbers = fields2extract.values
			else
				fields = line.split("\t")
				record = field_numbers.map{|n| fields[n]}
				if fields2extract[:id_col].nil?
					id = "rec_#{count}" #generate ids
				else
					id = record.shift
				end
				if !record[0].nil?
					record[0] = record[0].split(options[:separator])
				else
					record[0] = []
				end
				record[2] = record[2].to_i if !options[:start_col].nil?
				record[3] = record[3].to_i if !options[:end_col].nil?
				query = records[id]
				if query.nil?
					records[id] = [record]
				else
					query << record
				end
			end
			count +=1
		end
		return records
	end

	def self.get_fields2extract(options)
		fields2extract = {}
		[:id_col, :ont_col, :chromosome_col, :start_col, :end_col, :sex_col].each do |field|
			col = options[field]
			if !col.nil?
				col = col.to_i if !options[:header]
				fields2extract[field] = col
			end
		end
		return fields2extract
	end

	def self.get_field_numbers2extract(field_names, fields2extract)
		fields2extract.each do |field, name|
			fields2extract[field] = field_names.index(name)
		end
	end

	def self.create_cohort(records, options)
		ont = Cohort.get_ontology(Cohort.act_ont)
		rejected_terms = []
		rejected_recs = []
		cohort = Cohort.new()
		records.each do |id, record|
			rec = record.first
			terms = rec.first
			if options[:names] # Translate hpo names 2 codes
				init_term_number = terms.length
				terms, rec_rejected_terms = ont.translate_names(terms)
				if !rec_rejected_terms.empty?
					STDERR.puts "WARNING: record #{id} has the unknown term NAMES '#{rec_rejected_terms.join(',')}'. Terms removed."
					rejected_terms.concat(rec_rejected_terms)
				end
				if terms.empty? && init_term_number > 0
					rejected_recs << id
					next
				end
			end
			if rec.length > 1 # there is genomic region attributes
				variants = record.map{|v| v[1..3] }
			else
				variants = [] # Not exists genomic region attributes so we create a empty array
			end
			other_attr = {}
			if options[:extracted_fields].include?(:sex_col) # Check for additional attributes. -1 is applied to ignore :id in extracted fields
				other_attr[:sex] = record.first[options[:extracted_fields].index(:sex_col) -1]
			end
			cohort.add_record([id, terms, check_variants(variants)], other_attr)
		end
		return cohort, rejected_terms.uniq, rejected_recs
	end

	def self.check_variants(vars)
		checked_vars = []
		vars.each do |var| #[chr, start, stop]
			if var.first == '-' # the chr must be defined
				STDERR.puts "WARNING: variant #{var.join(',')} has been removed"
			else
				checked_vars << var
			end
		end
		return vars
	end
end
