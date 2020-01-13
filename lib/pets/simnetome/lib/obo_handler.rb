# @author Fernando Moreno Jabato <jabato(at)uma(dot)es>
# @description Class to handle OBO files (based on OBO 1.4 format)

#########################################################
# AUTHOR NOTES
#########################################################

# 1 - Handle "consider" values
# 2 - Handle "replaced_by" values

class OBO_Handler
	#############################################
	# FIELDS
	#############################################
		# @structureType :: used to know if ontology structure is hierarchical or circular
		# @header :: hash with OBO header info
		# @terms :: hash with terms found
		# @typedefs :: hash with typedefs found
		# @instances :: hash with instances found
		# @parents :: hash with parentals for each term found
		# @obsoletes :: array with obsolote term ids
		# @expansions :: hash of hashes with expansion applied by tag for each term found. Updated using expand_terms_by_tag()
		# @frequencies :: hash of arrays with frequencies for each term. Structural = childs + 1; Asigned = given by the user
		# @max_freqs :: array with maximum frequencies observed
		# @alt_ids :: hash ids which correspond to alternative ids
		# @ics :: already calculated ICs
		# @obsoletes :: array with obsoletes IDs

	#############################################
	# CONSTRUCTOR
	#############################################
	def initialize(file = nil, load = false)
		# Initialize object variables
		@structureType = nil
		@header = nil
		@terms = nil
		@typedefs = nil
		@instances = nil
		@parents = nil
		@obsoletes = nil
		@expansions = {}
		@info = nil
		@file = file
		@frequencies = nil
		@alt_ids = {}
		@max_freqs = nil
		@ics = {}
		@obsoletes = []

		# Load if proceed
		info = self.load() if (!file.nil?) & (load)
	end

	

	#############################################
	# CLASS METHODS
	#############################################

	# Class method to transform string with <tag : info> into hash structure
	# Param:
	# +info+:: string with info to be transformed into hash format
	# Return info stored into hash structure
	def self.info2hash(info)
		# Check special cases
		return nil if info.nil?
		return nil if !info.is_a? Array
		return nil if info.length <= 0
		# Load info
		info_hash = {}
		info.each do |tag_tuple|
			# Check special cases
			raise TypeError, 'Info element is NIL' if tag_tuple.nil?
			raise TypeError, 'Info element is not a string' if !tag_tuple.is_a? String
			raise 'Info element is empty string' if tag_tuple.length <= 0
			# Split info
			tag, value = tag_tuple.split(':',2)
			# Check
			raise EncodingError, 'Info element incorrect format' if (tag.nil?) | (value.nil?)
			# Prepare
			tag.lstrip!
			value.lstrip!
			# Store
			if info_hash.keys.include? tag
				if !info_hash[tag].kind_of?(Array)
					info_hash[tag] = [info_hash[tag]]
				end
				info_hash[tag] << value				
			else
				info_hash[tag] = value
			end
		end
		return info_hash
	end


	# Class method to load an OBO format file (based on OBO 1.4 format). Specially focused on load
	# the Header, the Terms, the Typedefs and the Instances.
	# Param:
	# +file+:: OBO file to be loaded
	# Returns the hash with all OBO file stored
	def self.load_obo(file)
		# Check special cases
		return nil if file.nil?
		return nil if !file.is_a? String
		return nil if file.length <= 0

		# Data variables
		header = ""
		terms = {}
		typedefs = {}
		instances = {}
		# Auxiliar variables
		infoType = "Header"
		currInfo = []
		stanzas = ["[Term]","[Typedef]","[Instance]"]
		# Read file
		File.open(file).each do |line|
			line.chomp!
			# Check if new instance is found
			if stanzas.include? line
				currInfo = self.info2hash(currInfo)
				id = currInfo.first[1]
				# Store current info
				case infoType
				when "Header" 
					header = currInfo
				when "Term"
					terms[id] = currInfo
				when "Typedef"
					typedefs[id] = currInfo
				when "Instance"
					instances[id] = currInfo
				end
				# Update info variables
				currInfo = []
				infoType = line.gsub!(/["\[\]"]/,"")
				next
			end
			# Concat info
			currInfo << line unless line.length <= 0
		end
		# Store last loaded info
		if currInfo.length > 0
			currInfo = info2hash(currInfo)
			id = currInfo.first[1]
			# Store current info
			case infoType
			when "Header" 
				header = currInfo
			when "Term"
				terms[id] = currInfo
			when "Typedef"
				typedefs[id] = currInfo
			when "Instance"
				instances[id] = currInfo
			end
		end

		return {"Header" => header, "Term" => terms, "Typedef" => typedefs, "Instance" => instances}
	end
	

	# Expand a (starting) term using a specific tag and return all extended terms into an array and
	# the relationship structuture observed (hierarchical or circular). If circular structure is
	# foumd, extended array will be an unique vector without starting term (no loops) 
	# Param:
	# +start+:: term where start to expand
	# +terms+:: set to be used to expand
	# +target_tag+:: tag used to expand
	# +eexpansion+:: already expanded info
	# +split_info_char+:: special regex used to split info (if it is necessary)
	# +split_info_indx+:: special index to take splitted info (if it is necessary)
	# +alt_ids+:: set of alternative IDs
	# Returns a vector with the observed structure (string) and the array with extended terms
	# Note: we extremly recomend use expand_by_tag function instead of it (directly)
	def self.expand_tag(start,terms,target_tag,expansion = {}, split_info_char = " ! ", split_info_indx = 0, alt_ids = [])
		# Check
		return nil if start.nil?
		return nil if terms.nil?
		return nil if !terms.is_a? Hash
		return nil if terms.length <= 0
		return nil if target_tag.nil?
		return nil if !target_tag.is_a? String
		return nil if target_tag.length <= 0
		return nil if expansion.nil?
		raise ArgumentError, 'Info_index can not be a negative number' if split_info_indx < 0
		# Take start term available info and already accumulated info
		current_expanded = expansion[start]
		current_expanded = [] if current_expanded.nil?
		start_expansion = terms[start][target_tag]
		return ["Source_Term",[]] if start_expansion.nil?
		start_expansion = start_expansion.clone
		start_expansion = Array.new(1,start_expansion) unless start_expansion.kind_of?(Array)

		# Prepare auxiliar variables
		visited = []
		struc = "Hierarchical"

		# Study direct extensions
		while start_expansion.length > 0
			# Take current element
			id = start_expansion.shift
			id = id.split(split_info_char)[split_info_indx]	
			id = alt_ids[id] if alt_ids.include? id # NOTE: if you want to persist current ID instead source ID, re-implement this
			# Handle
			if current_expanded.include? id # Check if already have been included into this expansion
				struct = "Circular" 
			elsif expansion.include? id # Check if current already has been expanded
				# Concat
				current_expanded << id 
				current_expanded = current_expanded | expansion[id]
				# Check circular case
				if current_expanded.include? start
					struct = "Circular"
					[id,start].each do |repeated| current_expanded.delete(repeated) end
				end	
			else # Expand
				# Add current
				current_expanded << id
				expansion[start] = current_expanded
				# Expand current
				structExp, expansionNew = OBO_Handler.send :expand_tag,id,terms,target_tag,expansion,split_info_char,split_info_indx,alt_ids
				# Concat
				current_expanded = current_expanded | expansionNew unless expansionNew.nil?
				# Check struct
				struct = "Circular" if structExp == "Circular"
				# Check circular case
				if (current_expanded.include? start)
					struct = "Circular"
					current_expanded.delete(start)
				end
			end
		end

		# Update
		expansion[start] = current_expanded

		# Return
		return struct, current_expanded
	end


	# Expand terms using a specific tag and return all extended terms into an array and
	# the relationship structuture observed (hierarchical or circular). If circular structure is
	# foumd, extended array will be an unique vector without starting term (no loops) 
	# Param:
	# +terms+:: set to be used to expand
	# +target_tag+:: tag used to expand
	# +split_info_char+:: special regex used to split info (if it is necessary)
	# +split_info_indx+:: special index to take splitted info (if it is necessary)
	# +alt_ids+:: set of alternative IDs
	# Returns a vector with the observed structure (string) and the hash with extended terms
	def self.expand_by_tag(terms,target_tag, split_info_char = " ! ", split_info_indx = 0, alt_ids = [])
		# Check special cases
		return nil if terms.nil?
		return nil if !terms.is_a? Hash
		return nil if terms.length <= 0
		return nil if target_tag.nil?
		return nil if !target_tag.is_a? String
		return nil if target_tag.length <= 0
		raise ArgumentError, 'Info_index can not be a negative number' if split_info_indx < 0

		# Define structure type
		structType = "Hierarchical"
		expansion = {}
		terms.each do |id,tags|
			# Check
			next if tags.nil?
			raise TypeError, 'Tags of term (#{id}) is not a hash' if !tags.is_a? Hash 
			# Check if target tag is defined
			if tags.keys.include? target_tag
				id = id.split(split_info_char)[split_info_indx]
				# Obtain related terms
				set_structure, related_ids = OBO_Handler.send :expand_tag, id, terms, target_tag, expansion, split_info_char, split_info_indx, alt_ids
				# Check structure
				if(set_structure == "Circular")
					structType = "Circular"
				end
				# Update Expansion info
				expansion[id] = related_ids
			end
		end

		# Check special case
		structType = "Atomic" if expansion.length <= 0
		structType = "Sparse" if (expansion.length > 0) & ((terms.length - expansion.length) >= 2)

		# Return type and hash with expansion
		return structType, expansion
	end


	#############################################
	# GENERAL METHODS
	#############################################
	
	# Increase the arbitrary frequency of a given term 
	# Params:
	# +term+:: to be updated
	# +increase+:: amount to be increased
	# Returns true if process ends without errors and false in other cases
	def add_observed_term(term,increase = 1.0)
		# Check
		raise ArgumentError, "Term given is NIL" if term.nil?
		return false if self.exist_term(term) == 0
		return false if @frequencies.nil?
		# Add frequency
		@frequencies[term] = [-1.0,0.0] if !@frequencies.include? term
		@frequencies[term][1] = 0 if @frequencies[term][1] == -1
		@frequencies[term][1] += increase
		# Check maximum frequency
		@max_freqs = [-1,@frequencies[term][1]] if @max_freqs.nil?
		@max_freqs[1] = @frequencies[term][1] if @frequencies[term][1] > @max_freqs[1]
		return true
	end


	# Increase the arbitrary frequency of a given term set 
	# Params:
	# +terms+:: set of terms to be updated
	# +increase+:: amount to be increased
	# Returns true if process ends without errors and false in other cases
	def add_observed_terms(terms, increase = 1.0)
		# Check
		raise ArgumentError, 'Terms array given is NIL' if terms.nil?
		raise ArgumentError, 'Terms given is not an array' if !terms.is_a? Array
		# Add observations
		checks = terms.map{|id| self.add_observed_term(id,increase)}
		return checks
	end


	# Check if a given term is handled into this OBO handler object
	# Param:
	# +term+:: to be checked
	# Returns non-zero element if exists and zero if is not contained. Negative value means obsolete term
	def exist_term(term)
		raise ArgumentError, 'Term given is nil' if term.nil?
		return 0 if !@terms.include? term
		return -1 if @obsoletes.include? term
		return 1
	end


	# Expand alternative IDs arround all already stored terms
	# Params:
	# +alt_tag+:: tag used to expand alternative IDs
	# Returns true if process ends without errors and false in other cases
	def expand_alternatives(alt_tag = "alt_id")
		# Check input
		return false if @terms.nil?
		# Take all alternative IDs
		terms_copy = @terms.keys
		terms_copy.each do |id|
			tags = @terms[id]
			next if tags.nil?
			if tags.keys.include? alt_tag
				# Check alternative ids
				alt_ids = tags[alt_tag]
				alt_ids = Array.new(1,alt_ids) if !alt_ids.kind_of? Array
				# Update info
				alt_ids.each do |alt_term|
					@alt_ids[alt_term] = id
					@terms[alt_term] = @terms[id]
					if !@parents.nil?
						@parents[alt_term] = parents[id] if parents.include? id
					end
				end
			end
		end
		self.expand_frequencies
		# Everything ok
		return true
	end


	# Executes basic expansions of tags (alternatives, obsoletes and parentals) with default values
	# Returns :: VOID method
	def expand_base()
		self.expand_alternatives
		self.expand_obsoletes
		self.expand_parentals
	end


	# Calculates regular frequencies based on ontology structure (using parentals)
	# Returns :: VOID method
	def expand_frequencies()
		# Check
		@frequencies = {} if @frequencies.nil?
		# Reset
		@frequencies.each do |id, freqs|
			next if freqs.nil?
			freqs[0] = 0
		end
		# Per each term, add frequencies
		@terms.each do |id, tags|
			next if @alt_ids.include? id
			# Add if it's necessary
			@frequencies[id] = [0.0,-1.0] unless @frequencies.include? id
			@frequencies[id][0] = 0.0 if @frequencies[id][0] < 0
			# Increase current frequencies
			@frequencies[id][0] += 1.0
			# Increase parental frequencies
			next if @parents.nil?
			if @parents.include? id
				@parents[id].each do |parent_id|
					@frequencies[parent_id] = [0.0,-1.0] unless @frequencies.include? parent_id
					@frequencies[parent_id][0] = 0.0 if @frequencies[parent_id][0] < 0
					@frequencies[parent_id][0] += 1.0						
				end
			end
		end
		# Found maximum frequency
		max_freq = -1.0
		@frequencies.each do |id,freqs|
			next if @alt_ids.include? id
			max_freq = freqs[0] if freqs[0] > max_freq
		end
		# Store
		@max_freqs = [-1.0,-1.0] if @max_freqs.nil? 
		@max_freqs[0] = max_freq
		# Add alternative IDs
		if @alt_ids.length > 0
			@alt_ids.each do |alt,source|
				@frequencies[alt] = @frequencies[source]
			end
		end
	end


	# Expand obsoletes set and link info to their alternative IDs
	# Params:
	# +obs_tag+::
	# ++::
	# ++::
	# Returns true if process ends without errors and false in other cases
	def expand_obsoletes(obs_tag ="is_obsolete",alt="replaced_by",reset_obsoletes=true)
		# Check
		return false if @terms.nil?
		# Reset
		@obsoletes = [] if reset_obsoletes
		# Check obsoletes
		@terms.each do |id|
			tags = @terms[id]
			next if tags.nil?
			if tags.keys.include? obs_tag
				# Check obsolete
				next if !tags[obs_tag]
				alt_id = nil
				# Check if alternative value is available
				alt_id = tags[alt] if tags.keys.include? alt
				# Store
				@alt_ids[id] = alt_id
				@obsoletes << id				
			end
		end
		# END
		return true
	end


	#
	#
	#
	def expand_parentals(tag = "is_a",split_info_char = " ! ", split_info_indx = 0)
		# Check
		return false if @terms.nil?
		# Expand
		structType, parentals = self.class.expand_by_tag(@terms,tag,split_info_char,split_info_indx,@alt_ids)
		# Check
		return false if (structType.nil?) | parentals.nil?
		# Store
		@parents = parentals
		@structureType = structType
		# Expand frequencies
		self.expand_frequencies
		# Finish		
		return true
	end


	#
	#
	#
	def get_ancestors(term)
		# Check
		raise ArgumentError, 'Term specified is NIL' if term.nil?
		return false if @parents.nil?

		# Find into parentals
		return @parents[term]		
	end


	#
	#
	#
	def get_IC(term,regular=true,force=false)
		# Check
		return nil if @frequencies.nil?
		return nil if @max_freqs.nil? # It shouldn't happen if frequencies are calculated
		raise ArgumentError, 'Term specified is NIL' if term.nil?
		return -1 if (@max_freqs[0] <= 0) & regular
		return -1 if (@max_freqs[1] <= 0) & !regular
		# Calculate
		return @ics[term][0] if (@ics.include? term) & regular & !force
		return @ics[term][1] if (@ics.include? term) & !regular & !force
		return -1 if !@frequencies.include? term
		term_ic_reg = -Math.log10(@frequencies[term][0].fdiv(@max_freqs[0])) # Natural logarithm
		term_ic_arb = -Math.log10(@frequencies[term][1].fdiv(@max_freqs[1])) # Natural logarithm
		# Store
		to_store = term
		to_store = @alt_ids[term] if @alt_ids.include? term
		@ics[to_store] = [term_ic_reg,term_ic_arb]
		@ics[term] = @ics[to_store] if @alt_ids.include? term
		# Return
		return term_ic_reg if regular
		return term_ic_arb
	end


	#
	#
	#
	def get_ICMICA(termA,termB,regular=true)
		mica = self.get_MICA(termA,termB,regular)
		return false if !mica
		return mica[1]
	end


	#
	#
	#
	def get_frequency(term,regular=true)
		# Check
		raise ArgumentError, 'Term specified is NIL' if term.nil?
		return false if @frequencies.nil?
		# Return
		return @frequencies[term][0] if regular
		return @frequencies[term][1]
	end


	#
	#
	#
	def get_MICA(termA,termB,regular=true)
		# Obtain ancestors (include itselfs too)
		anc_A = self.get_ancestors(termA) | [termA]
		anc_B = self.get_ancestors(termB) | [termB]
		# Check
		return false if (!anc_A) | (!anc_B)
		# Find shared ancestors
		shared_ancestors = anc_A & anc_B
		return [nil,-1] if shared_ancestors.length <= 0
		# Find MICA
		mica = [nil,-1]
		shared_ancestors.each do |anc|
			# Obtain IC
			ic = self.get_IC(anc,regular)
			# Check
			mica = [anc,ic] if ic > mica[1]
		end
		# Return value calculated
		return mica
	end

	#
	#
	#
	def get_similarity(termsA, termsB, regular=true)
		# Check
		return nil if termsA.nil? | termsB.nil?
		return nil if (!termsA.is_a? Array) | (!termsB.is_a? Array)
		return nil if (termsA.length <= 0) | (termsB.length <= 0) 
		return false if @frequencies.nil?
		# Prepare necessary values
		micasA = []
		micasB = []
		# Launch comparissons A -> B
		termsA.each do |tA|
			micas = termsB.map{|tB| self.get_ICMICA(tA,tB,regular)}
			# Remove special cases
			[false,nil].each do |err_value| micas.delete(err_value) end
			# Obtain maximum value
			micasA << micas.max if micas.length > 0
			micasA << 0 if micas.length <= 0
		end
		# Launch comparissons B -> A
		termsB.each do |tB|
			micas = termsA.map{|tA| self.get_ICMICA(tA,tB,regular)}
			# Remove special cases
			[false,nil].each do |err_value| micas.delete(err_value) end
			# Obtain maximum value
			micasB << micas.max if micas.length > 0
			micasB << 0 if micas.length <= 0
		end
		# Calculate means
		means_sim = [micasA.inject{ |sum, el| sum + el }.to_f / micasA.size , micasB.inject{ |sum, el| sum + el }.to_f / micasB.size]
		# Return similarity
		return means_sim.inject{ |sum, el| sum + el }.to_f / means_sim.size
	end


	# Method used to load information stored into an OBO file and store it into this object.
	# If a file is specified by input parameter, current @file value is updated
	# Param
	# +file+:: optional file to update object stored file
	# Return true if process ends without errors, false in other cases
	def load(file = nil)
		# Update
		@file = file unless file.nil?
		# Check special cases
		return false if @file.nil?
		# Load
		info = self.class.load_obo(@file)
		# Check special cases
		return false if info.nil?
		# Store
		@info = info
		@header = info["Header"]
		@terms = info["Term"]
		@typedefs = info["Typedef"]
		@instances = info["Instance"]
		# Return
		return true
	end


	#############################################
	# ACCESS CONTROL
	#############################################
	## METHODS
	# public 
	# private
	private_class_method :expand_tag

	## ATTRIBUTES
	attr_reader :info, :file, :structureType, :header, :terms, :typedefs, :instances, :parents, :obsoletes, :expansions, :frequencies, :alt_ids, :max_freqs, :ics, :obsoletes
	# attr_writer 

end