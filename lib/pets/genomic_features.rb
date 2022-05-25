class Genomic_Feature
	def self.array2genomic_feature(arr)
		new(arr.map{|r| yield(r)})
	end

	def self.hash2genomic_feature(h)
		vars = []
		h.each do |h, v|
			vars << yield(h, v)
		end
		new(vars)
	end

	#If any method use gen_fet as name is a Genomic_Feature object
	def initialize(feat_array) # [[chr1, start1, stop1],[chr1, start1, stop1]]
		@regions = {}
		@reg_by_to = {}
		@reg_id = -1
		load_features(feat_array)
	end

	def load_features(feat_array)
		feat_array.each do |chr, start, stop, to|
			chr = chr.to_sym
			@reg_id +=1
			id = to.nil? ? @reg_id : to
			region = {chr: chr, start: start, stop: stop, to: id }
			@reg_by_to[id] = region
			add_record(@regions, chr, region)
		end
	end

	def length
		return @regions.length
	end

	def each_chr()
		@regions.each do |chr, regs|
			yield(chr, regs)
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

	def get_chr_regs(chr)
		return @regions[chr]
	end

	def region_by_to(to)
		return @reg_by_to[to]
	end

	def get_sizes
		sizes = []
		each do |chr, region|
			size = region[:stop] - region[:start] + 1
			sizes << size
		end
		return sizes
	end

	def match(other_gen_feat)
		all_matches = {}
		each_chr do |chr, regs|
			other_regs = other_gen_feat.get_chr_regs(chr)
			next if other_regs.nil?
			local_matches = []
			regs.each do |reg|
				start = reg[:start] 
				stop = reg[:stop] 
				other_regs.each do |other_reg|
					local_matches << other_reg[:to] if coor_overlap?(start, stop, other_reg)
				end
				all_matches[reg[:to]] = local_matches if !local_matches.empty?
			end
		end
		return all_matches
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
				overlap = coor_overlap?(start, stop, reg)
				reg_ids << reg[:to] if overlap
			end
			overlaps << reg_ids.uniq
		end
		return overlaps
	end

	def generate_cluster_regions(meth, tag, ids_per_reg = 1, obj = false)
		compute_windows(meth) # Get putative genome windows
		ids_by_cluster = {}
		annotated_full_ref = [] # All reference windows wit uniq id and chr tagged
		@regions.each do |chr, regs|
			reference = @windows[chr]
			overlaps = get_reference_overlaps(regs, reference)
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
		annotated_full_ref = Genomic_Feature.array2genomic_feature(annotated_full_ref){|r| [r[2], r[0], r[1], r[3]]} if obj
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
		single_nt = []
		genomic_ranges.each do |gr|
			start = gr[:start]
			stop = gr[:stop]
			if stop - start > 0
				reference << start # get start
				reference << stop # get stop
			else # Build a window of at least one nt for snv
				single_nt << start
			end
		end
		reference.uniq!
		single_nt.each do |snt| # add start stop for snv
			reference << snt
			reference << snt
		end
		reference.sort!
		#Define overlap ranges
		final_reference = []
		last_len = 1
		reference.each_with_index do |coord,i|
			next_coord = reference[i + 1]
			if !next_coord.nil?
				current_len = next_coord - coord
				coord = coord + 1 if last_len == 0 # Separate SNV window from others
				if current_len == 0 && last_len > 0 && !final_reference.empty?
					final_reference.last[1] -= 1 # Separate SNV window from others
				end
				final_reference << [coord, next_coord]  
				last_len = current_len
			end
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
