require 'obo_handler'

class Ontology
	attr_reader :term_level, :excluded_codes

	def initialize()
		@excluded_codes = []
		@ont_data = {}
		@term_level = {}
		@term_parent_child_relations = {}
		@name2code_dictionary = {}
		@profiles = []
		@onto_ic = {}
		@freq_ic = {}
		@onto_ic_profile = []
	  	@freq_ic_profile = []
	end

	def load_black_list(excluded_codes_file)
		File.open(excluded_codes_file).each do |line|
			line.chomp!
			@excluded_codes << line
		end
	end

	def load_data(obo_path)
		load_hpo_file(obo_path)
		get_child_parent_relations
		create_hpo_dictionary
		extract_ontology_levels_info
	end

	def load_hpo_file(hpo_file)
		hpo_obsolete = {}
		id = nil
		name = nil
		alt_ids = []
		syn = []
		is_a = []
		is_obsolete = false
		new_ids = []
		allowed_tags = %w[id name is_a synonym alt_id is_obsolete replaced_by consider]
		File.open(hpo_file).each do |line|
			line.chomp!
			tag, info = line.split(': ')
			if allowed_tags.include?(tag)
				if tag == 'id'
					if is_obsolete
						new_ids.each do |new_id|
							hpo_obsolete[id] = [new_id, name]
						end
					else
						add_record2storage(id, name, is_a, syn, alt_ids) if !name.nil?
					end
					id = info
					name = nil
					alt_id = []
					syn = []
					is_a = []
					is_obsolete = false
					new_ids = []
				end
				if tag == 'alt_id'
					alt_ids << info
				elsif tag == 'is_a'
					is_a << info.split(' ! ')[0]
				elsif tag == 'synonym'
					syn << info.split('"')[1] #to keep only the name of the synonym
				elsif tag == 'is_obsolete'
					is_obsolete = true
				elsif tag == 'replaced_by'
					new_ids << info
				elsif tag == 'consider' #there can be more than one "consider" tag for a given term
					new_ids << info
				else
					name = info
				end
			end
		end
		add_record2storage(id, name, is_a, syn, alt_ids)
		new_ids.each do |new_id|
			hpo_obsolete[id] = [new_id, name] if is_obsolete
		end
		hpo_obsolete.each do |obsoleteID, data|
			new_id = data[0]
			obsolete_name = data[1]
			info_for_obsolete = @ont_data[new_id]
			unless info_for_obsolete.nil?
				if info_for_obsolete[3].nil?
					info_for_obsolete[3] = [obsolete_name]
				else
					info_for_obsolete[3] << obsolete_name unless info_for_obsolete[3].include?(obsolete_name)
				end
			end
			info_for_obsolete = @ont_data[new_id]
			@ont_data[obsoleteID] = info_for_obsolete
		end
	end

	def add_record2storage(id, name, is_a, syn, alt_ids)
		if !@excluded_codes.include?(id)
			attributes = [id, name, is_a - @excluded_codes, syn]
			@ont_data[id] = attributes
			alt_ids.each do |altid|
				@ont_data[altid] = attributes
			end
		end 
	end

	def get_child_parent_relations
	# for getting hpo childs
		@ont_data.each do |hpo_code, hpo_data|
			id, name, is_a, syn = hpo_data
			next if is_a.nil?
			hpo_child = [id, name]
			is_a.each do |par_hpo_code|
				query = @term_parent_child_relations[par_hpo_code]
				if query.nil?
					@term_parent_child_relations[par_hpo_code] = [hpo_child]
				else
					query << hpo_child
				end
			end
		end
	end

	def extract_ontology_levels_info
	  @ont_data.each do |hpo_id, hpo_data|
	    parental_terms = hpo_data[2]
	    unless parental_terms.empty?
	      search_for_parentals(parental_terms.first, 0, hpo_id)
	    else
	      @term_level[0] = [hpo_id] #HP:0000001 or #HP:0000118, parental term in the HPO
	    end
	  end
	end

	def get_ontology_levels_from_profiles(uniq=true)
	  profile_terms = @profiles.flatten
	  profile_terms.uniq! if uniq

	  cohort_hpo_levels = {}
	  profile_terms.each do |term|
	    @term_level.each do |level, hpo_ids|
	      hpo_ids.each do |hpo|
	        if term == hpo
	          query = cohort_hpo_levels[level]
	          if query.nil?
	            cohort_hpo_levels[level] = [term]
	          else
	            query << term
	          end
	        end
	      end
	    end
	  end
	  return cohort_hpo_levels
	end

	def create_hpo_dictionary
		@ont_data.each do |hpo, metadata|
			hpo_code, hpo_name, hpo_parents, hpo_synonyms = metadata 
			@name2code_dictionary[hpo_name] = hpo_code
			next if hpo_synonyms.nil? # To remove hpos without parents (i.e: Obsolete HPO)
			hpo_synonyms.each do |syn|
				@name2code_dictionary[syn] = hpo_code
			end
		end
	end

	def translate_names2codes(hpos) 
	  hpo_codes = []
	  rejected_hpos = []
	  hpos.each_with_index do |hpo_name, i|
	    hpo_code = @name2code_dictionary[hpo_name]
	    if hpo_code.nil?
	      rejected_hpos << hpo_name
	    else
	      hpo_codes << hpo_code
	    end
	  end
	  return hpo_codes, rejected_hpos
	end

	def translate_codes2names(terms)
		term_names = []
		rejected_codes = []
	    terms.each do |term|
	      term_data = @ont_data[term]
	      if term_data.nil?
	        STDERR.puts "WARNING: ontology code '#{term}' not exists."
	        rejected_codes << term
	      else
	        term_names << term_data[1]
	      end
	    end
	    return term_names, rejected_codes
	end

	def check_codes(hpos)
	  checked_codes = []
	  rejected_hpos = []
	  hpos.each do |hpo_code|
	    hpo_data = @ont_data[hpo_code]
	    if hpo_data.nil?
	      rejected_hpos << hpo_code
	    else
	      main_hpo_code, name = hpo_data
	      checked_codes << main_hpo_code # change from alternate hpo codes to the main ones
	    end
	  end
	  return checked_codes, rejected_hpos
	end

	def get_parents(term)
		parents = []
		term_data = @ont_data[term]
		term_data[2].each do |par_hpo_code, par_hpo_name|
			parents << par_hpo_code
			parents.concat(get_parents(par_hpo_code))
		end
		return parents
	end

	def get_term_names_from_profiles(profs = [])
	  profs = @profiles if profs.empty?
	  profiles_with_names = []
	  profs.each do |profile|
	  	names = []
	    profile.each do |hpo|
	      hpo_data = @ont_data[hpo]
	      if hpo_data.nil?
	        STDERR.puts "WARNING: hpo code '#{hpo}' not exists."
	      else
	        names << hpo_data[1]
	      end
	    end
	    profiles_with_names << names
	  end
	  return profiles_with_names
	end

	def get_term_frequency_from_profiles(names=true)
	  stats = Hash.new(0)
	  @profiles.each do |profile|
	    profile.each do |term|
	      stats[term] += 1
	    end
	  end
	  n_profiles = @profiles.length
	  freqs = []
	  stats.each do |term, count|
	  	term = @ont_data[term][1] if names
	    freqs << [term, count.fdiv(n_profiles)*100]
	  end
	  freqs.sort!{|h1, h2| h2[1] <=> h1[1]}
	  return freqs
	end

	def get_more_specific_childs_table(hpo_codes)
	  more_specific_hpo = []
	  hpo_codes.each do |hpo_code|
	    hpo_data = @ont_data[hpo_code]
	    main_hpo_code, name = hpo_data
        childs = @term_parent_child_relations[main_hpo_code]
        childs.nil? ? specific_childs = [] : specific_childs = childs
        more_specific_hpo << [[main_hpo_code, name], specific_childs]
      end
      return more_specific_hpo
	end

	def load_profiles(profiles)
		@profiles.concat(profiles)
	end

	def get_profile_mean_length
		prof_lengths = @profiles.map{|p| p.length}
  		return prof_lengths.inject(0){|sum, n| sum + n}.fdiv(@profiles.length).round(4)
  	end

  	def get_profile_length_at_percentile(perc=50)
  	  percentile = (100 - perc).fdiv(100)
  	  prof_lengths = @profiles.map{|p| p.length}.sort
  	  n_profiles = @profiles.length 
  	  percentile_length = nil
	  rate = 0
	  count = 0
	  while rate <= percentile
	    percentile_length = prof_lengths[count+1]
	    rate = count.fdiv(n_profiles)
	    count += 1
	  end
	  return percentile_length
  	end

	def get_ic_by_onto_and_freq(hpo_file)
	  obof = OBO_Handler.new(hpo_file, true) # Load ontology
	  obof.expand_base
	  @profiles.each do |profile|
	    obof.add_observed_terms(profile)
	  end
	  hpos = @profiles.flatten.uniq
	  onto_ic_values = hpos.map{|code| obof.get_IC(code)}
	  freq_ic_values = hpos.map{|code| obof.get_IC(code, false)}
	  hpos.each_with_index do |code, i|
	    @onto_ic[code] = onto_ic_values[i]
	  end 
	  hpos.each_with_index do |code, i|
	    @freq_ic[code] = freq_ic_values[i]
	  end 
	  return @onto_ic, @freq_ic
	end

	def get_ic_profile_by_onto_and_freq
	  @profiles.each do |profile|
	    pf_len = profile.length
	    @onto_ic_profile << (profile.map{|code| @onto_ic[code]}.inject(0){|sum, val| sum +val}).fdiv(pf_len)
	    @freq_ic_profile << (profile.map{|code| @freq_ic[code]}.inject(0){|sum, val| sum +val}).fdiv(pf_len)
	  end
	  return @onto_ic_profile, @freq_ic_profile
	end

	private

	def search_for_parentals(parental_id, counter, hpo_id)
	  hpo_data = @ont_data[parental_id] 
	  if !hpo_data.nil?
	    parental_id = hpo_data[2].first
	    counter += 1
	    search_for_parentals(parental_id, counter,  hpo_id)
	  else
	    query = @term_level[counter]
	    if query.nil?
	      @term_level[counter] = [hpo_id]
	    else
	      query << hpo_id
	    end
	  end
	end
end
