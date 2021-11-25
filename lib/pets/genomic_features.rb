class Genomic_Feature
	def initialize(feat_array) # [[chr1, start1, stop1],[chr1, start1, stop1]]
		@regions = []
		load_features(feat_array)
	end

	def load_features(feat_array)
		feat_array.each do |chr, start, stop|
			@regions << {chr: chr.to_sym, start: start, stop: stop}
		end
	end

	def each()
		@regions.each do |region|
			yield(region)
		end
	end
end
