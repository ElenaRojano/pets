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