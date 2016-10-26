#Common methods for predictors
#Training file example = 9  131371492   131375954   HP:0010974  2.41161970596
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

def load_training_file4HPO(training_file)
	training_set = {}	
	posInfo = loadFile(training_file)
	posInfo.each do |info|
		hpoCode = info.delete_at(3)
		query = training_set[hpoCode]
		if query.nil?
			training_set[hpoCode] = [info]
		else
			query << info
		end
	end
	return training_set
end

#3. Load info file:

def loadFile(file)
	information = []
	File.open(file).each do |line|
		line.chomp!
		allInfo = line.split("\t")
		chr = allInfo[0]
		startPos = allInfo[1].to_i
		stopPos = allInfo[2].to_i
		hpoCode = allInfo[3]
		associationValue = allInfo[4].to_f
		information << [chr, startPos, stopPos, hpoCode, associationValue]
	end
	return information
end




