def process_patient_data(patient_data)
	parsed_patient_data = {}
	patient_data.each do |patientID, metadata|
		phenotypes, chr, start, stop = metadata
		info = [patientID, start, stop]
		query = parsed_patient_data[chr]
		if query.nil?
			parsed_patient_data[chr] = [info]
		else
			query << info
		end
	end
	return parsed_patient_data
end

def get_final_coverage(raw_coverage, bin_size)
	coverage_to_plot = []
	raw_coverage.each do |chr, coverages|
		coverages.each do |start, stop, coverage|
			bin_start = start - start % bin_size
			bin_stop = stop - stop%bin_size
			while bin_start < bin_stop
				coverage_to_plot << [chr, bin_start, coverage]
				bin_start += bin_size
			end
		end
	end
	return coverage_to_plot
end

def get_sor_length_distribution(raw_coverage)
	all_cnvs_length = []
	cnvs_count = []
	raw_coverage.each do |chr, coords_info|
		coords_info.each do |start, stop, pat_records|
			region_length = stop - start + 1
			all_cnvs_length << [region_length, pat_records]
		end
	end
	all_cnvs_length.sort!{|c1, c2| c1[1] <=> c2[1]}
	return all_cnvs_length
end

def get_cnvs_length(patient_data)
	length_stats = Hash.new(0)
	patient_data.each do |pat_id, patient_record|
    	string_hpos, chr, start, stop = patient_record
    	length_stats[stop - start] += 1
    end
    return length_stats.to_a.sort!{|stat| stat[1] <=> stat[1] }
end


def calculate_coverage(regions_data, delete_thresold = 0)
	raw_coverage = {}
	n_regions = 0
	patients = 0
	nt = 0
	regions_data.each do |start, stop, chr, node|
		number_of_patients = node.split('.').last.to_i
		if number_of_patients <= delete_thresold
			number_of_patients = 0
		else
			n_regions += 1
			patients += number_of_patients
			nt += stop - start			
		end
		coords = [start, stop, number_of_patients]
		query = raw_coverage[chr]
		if query.nil?
			raw_coverage[chr] = [coords]
		else
			query << coords
		end
	end
	return raw_coverage, n_regions, nt, patients.fdiv(n_regions)
end