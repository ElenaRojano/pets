#Common methods for predictors
#Training file example = 9  131371492   131375954   HP:0010974  2.41161970596	9.3.A.5
#1. Indexing by chr (region)

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

def load_hpo_dictionary_name2code(hpo_file)
	storage = {}
	#STDERR.puts hpo_file.inspect
	File.open(hpo_file).each do |line|
		line.chomp!
		fields = line.split("\t")
		hpo_code = fields.shift
		phenotype = fields.shift
		synonyms = fields.shift
		storage[phenotype] = hpo_code
		if !synonyms.nil?
			synonyms.split('|').each do |syn|
				storage[syn] = hpo_code
			end
		end
	end
	return storage
end

def load_hpo_metadata(hpo_file)
	#input_file: hpo2name.txt
	storage = {}
	File.open(hpo_file).each do |line|
		line.chomp!
		fields = line.split("\t")
		hpo_code = fields.shift
		phenotype = fields.shift
		synonyms = fields.shift 
		relations = []
		fields.each_slice(2) do |pair|
			#pair = HPO code, phenotype
			relations << pair
		end
		storage[hpo_code] = [phenotype, relations]
	end
	return storage
end

def inverse_hpo_metadata(hpo_metadata)
	# for getting hpo childs
	storage_child = {}
	hpo_metadata.each do |hpo_code, hpo_data|
		hpo_name, parents = hpo_data
		parents.each do |par_hpo_code, par_hpo_name|
			query = storage_child[par_hpo_code]
			hpo_child = [hpo_code, hpo_name]
			if query.nil?
				storage_child[par_hpo_code] = [par_hpo_name, [hpo_child]]
			else
				query.last << hpo_child
			end
		end
	end
	return storage_child
end

def load_hpo_ci_values(information_coefficient_file)
	hpos_ci_values = {}
	File.open(information_coefficient_file).each do |line|
		line.chomp!
		hpo_code, ci = line.split("\t")
		hpos_ci_values[hpo_code] = ci.to_f
	end
	return hpos_ci_values
end