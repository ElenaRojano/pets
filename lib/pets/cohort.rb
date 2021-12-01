require 'semtools'

class Cohort
	@@ont = {}
	class << self # https://www.ruby-forum.com/t/attr-accessor-for-class-variable/136693
		attr_accessor :act_ont # Which ontology use for ont related operations
	end

	attr_accessor :profiles

	def self.get_ontology(ont_id)
		return @@ont[ont_id]
	end

	def self.load_ontology(ont_name, ont_file, excluded_terms_file)
		ont = nil
		if !ont_file.include?('.json')
			if !excluded_terms_file.nil?
				ont = Ontology.new(file: ont_file, load_file: true, removable_terms: read_excluded_ont_file(excluded_terms_file))
			else
				ont = Ontology.new(file: ont_file, load_file: true)
			end
		else
			ont = Ontology.new
			ont.read(ont_file)
			if !excluded_terms_file.nil?
				ont.add_removable_terms(read_excluded_ont_file(excluded_terms_file))
				ont.remove_removable()
				ont.build_index()
			end
		end
		@@ont[ont_name] = ont
	end

	def self.read_excluded_ont_file(file)
		excluded_hpo = []
		File.open(file).each do |line|
			excluded_hpo << line.chomp
		end
		return excluded_hpo
	end

	def initialize()
		@profiles = {}
		@vars = {}
		@var_idx = Genomic_Feature.new([])
	end

	def add_record(rec) #[id, [profile], [[chr1, start1, stop1],[chr1, start1, stop1]]]
		id, profile, vars = rec
		@profiles[id] = profile.map{|t| t.to_sym} if !profile.nil? 
		add_gen_feat(id, vars) if !vars.nil?
	end

	def add_gen_feat(id, feat_array) # [[chr1, start1, stop1],[chr1, start1, stop1]]
		@vars[id] = Genomic_Feature.new(feat_array)
	end

	def get_profile(id)
		return @profiles[id]
	end

	def get_vars(id)
		return @vars[id]
	end

	def each_profile()
		@profiles.each do |id, profile|
			yield(id, profile)
		end
	end

	def each_var()
		@vars.each do |id, var_info|
			yield(id, var_info)
		end
	end

	def check(hard=false) # OLD format_patient_data
		ont = @@ont[Cohort.act_ont]
		rejected_terms = []
		rejected_recs = []
		@profiles.each do |id, terms|
			if hard
				terms = ont.clean_profile_hard(terms)
				rejec_terms = []
			else
				terms, rejec_terms = ont.check_ids(terms)
			end
			if !rejec_terms.empty?
				STDERR.puts "WARNING: record #{id} has the unknown CODES '#{rejec_terms.join(',')}'. Codes removed."
				rejected_terms.concat(rejec_terms)
			end
			if terms.empty?
				rejected_recs << id
			else
				@profiles[id] = terms
			end
		end
		@profiles.select!{|id, record| !rejected_recs.include?(id)}
		@vars.select!{|id, record| !rejected_recs.include?(id)}
		return rejected_terms.uniq, rejected_recs
	end

	def link2ont(ont_id)
		@@ont[ont_id].load_profiles(@profiles)
	end

	def get_profile_redundancy
		ont = @@ont[Cohort.act_ont]
		profile_sizes, parental_terms_per_profile = ont.get_profile_redundancy
		return profile_sizes, parental_terms_per_profile
	end

	def get_profiles_terms_frequency(options={})
		ont = @@ont[Cohort.act_ont]
		term_stats = ont.get_profiles_terms_frequency(**options) #https://www.ruby-lang.org/en/news/2019/12/12/separation-of-positional-and-keyword-arguments-in-ruby-3-0/
		return term_stats
	end

	def compute_term_list_and_childs()
		ont = @@ont[Cohort.act_ont]
		suggested_childs, term_with_childs_ratio = ont.compute_term_list_and_childs()
	end

	def get_profile_ontology_distribution_tables()
		ont = @@ont[Cohort.act_ont]
		ontology_levels, distribution_percentage = ont.get_profile_ontology_distribution_tables
		ontology_levels.unshift(["level", "ontology", "cohort"])
		distribution_percentage.unshift(["level", "ontology", "weighted cohort", "uniq terms cohort"])
		return ontology_levels, distribution_percentage
	end

	def get_ic_analysis()
		ont = @@ont[Cohort.act_ont]
		onto_ic, freq_ic = ont.get_observed_ics_by_onto_and_freq # IC for TERMS
		onto_ic_profile, freq_ic_profile = ont.get_profiles_resnik_dual_ICs # IC for PROFILES
		return onto_ic, freq_ic, onto_ic_profile, freq_ic_profile
	end

	def get_profiles_mean_size
		ont = @@ont[Cohort.act_ont]
		profile_mean_size = ont.get_profiles_mean_size
		return profile_mean_size
	end

	def get_profile_length_at_percentile(perc=50, increasing_sort: false)
		ont = @@ont[Cohort.act_ont]
		length_percent = ont.get_profile_length_at_percentile(perc=perc, increasing_sort: increasing_sort)
		return length_percent
	end

	def get_dataset_specifity_index(type)
		ont = @@ont[Cohort.act_ont]
		dsi = ont.get_dataset_specifity_index(type)
		return dsi
	end

	def compare_profiles(options={})
		ont = @@ont[Cohort.act_ont]
		similarities = ont.compare_profiles(**options)
		return similarities
	end

	def index_vars # equivalent to process_patient_data
		each_var do |id, var|
			@var_idx.merge(var, id)
		end
	end

	def get_vars_sizes(summary=false)
		if summary
			return @var_idx.get_summary_sizes
		else
			return @var_idx.get_sizes
		end
	end

	def generate_cluster_regions(meth, tag, lim)
		@var_idx.generate_cluster_regions(meth, tag, lim)
	end
end
