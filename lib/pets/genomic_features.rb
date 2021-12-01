class Genomic_Feature
	#If any method use gen_fet as name is a Genomic_Feature object
	def initialize(feat_array) # [[chr1, start1, stop1],[chr1, start1, stop1]]
		@regions = {}
		@reg_id = -1
		load_features(feat_array)
	end

	def load_features(feat_array)
		feat_array.each do |chr, start, stop|
			chr = chr.to_sym
			region = {start: start, stop: stop, to: @reg_id +=1 }
			add_record(@regions, chr, region)
		end
	end

	def each()
		@regions.each do |chr, regs|
			regs.each do |region|
				yield(chr, region)
			end
		end
	end

	def get_chr
		return @regions.keys
	end

	def get_sizes
		sizes = []
		each do |chr, region|
			size = region[:stop] - region[:start] + 1
			sizes << size
		end
		return sizes
	end

	def get_summary_sizes
		sizes = Hash.new(0)
		each do |chr, region|
			size = region[:stop] - region[:start] + 1
			sizes[size] += 1
		end
		return sizes.to_a.sort!{|s| s[1] <=> s[1] }
	end

	def merge(gen_fet, to = nil) # 'to' the regions must be connected "to" given id
		gen_fet.each do |chr, region|
			to.nil? ? region[:to] = @reg_id +=1 : region[:to] = to # handle id custom or default
			add_record(@regions, chr, region)
		end
	end

	def get_reference_overlaps(genomic_ranges, reference)
		overlaps = []
		reference.each do |start, stop|
			reg_ids = []
			genomic_ranges.each do |reg|
				reg_ids << reg[:to] if coor_overlap?(start, stop, reg)
			end
			overlaps << reg_ids.uniq
		end
		return overlaps
	end

	def generate_cluster_regions(meth, tag, ids_per_reg = 1)
		compute_windows(meth) # Get putative genome windows
		patients_out_of_cluster = 0
		ids_by_cluster = {}
		annotated_full_ref = [] # All reference windows wit uniq id and chr tagged
		@regions.each do |chr, regs|
			reference = @windows[chr]
			overlaps = get_reference_overlaps(regs, reference) # See what patient has match with a overlap region
			clust_numb = 0
			reference.each_with_index do |ref, i|
				current_ids = overlaps[i]
				if current_ids.length > ids_per_reg
					clust_id = "#{chr}.#{clust_numb +=1}.#{tag}.#{current_ids.length}"
					current_ids.each do |curr_id|
						add_record(ids_by_cluster, curr_id, clust_id, true)
					end
					annotated_full_ref << ref.dup.concat([chr, clust_id])
				end
			end
		end
		return ids_by_cluster, annotated_full_ref
	end

	def compute_windows(meth)
		@windows = {}
		@regions.each do |chr, regs|
			chr_windows = nil
			if meth == :reg_overlap
				chr_windows = compute_region_overlap_windows(regs)
			end
			@windows[chr] = chr_windows
		end
	end

	private
	
	def add_record(hash, key, record, uniq=false)
		query = hash[key]
		if query.nil?
			hash[key] = [record]
		elsif !uniq # We not take care by repeated entries
			query << record
		elsif !query.include?(record) # We want uniq entries
			query << record
		end
	end

	def compute_region_overlap_windows(genomic_ranges)
		reference = []
		reference.concat(genomic_ranges.map{|gr| gr[:start]})# get start
		reference.concat(genomic_ranges.map{|gr| gr[:stop]})# get stop
		reference.uniq!
		reference.sort!
		#Define overlap ranges
		final_reference = []
		reference.each_with_index do |coord,i|
			next_coord = reference[i + 1]
			final_reference << [coord, next_coord] if !next_coord.nil? 
		end
		return final_reference
	end

	def coor_overlap?(start, stop, reg)
		overlap = false
		reg_start = reg[:start]
		reg_stop = reg[:stop]
		if (start <= reg_start && stop >= reg_stop) ||
			(start > reg_start && stop < reg_stop) ||
			(stop > reg_start && stop <= reg_stop) ||
			(start >= reg_start && start < reg_stop)
			overlap = true
		end
		return overlap
	end
end
