#! /usr/bin/env ruby

require 'optparse'
require 'semtools'

######################################################################################
## METHODS
######################################################################################
def load_tabular_file(file, skip = 0)
  records = []
  File.open(file).each do |line|
    line.chomp!
    fields = line.split("\t")
    records << fields
  end
  records.shift(skip) unless skip == 0
  return records
end

def pairs2profiles(file)
  id2prof = {}
  file.each do |fields|
    load_value(id2prof, fields.first, fields.last)
  end
  return id2prof
end

def load_value(hash_to_load, key, value, unique = true)
   	query = hash_to_load[key]
    if query.nil?
        value = [value] if value.class != Array
        hash_to_load[key] = value
    else
        if value.class == Array
            query.concat(value)
        else
            query << value
        end
        query.uniq! unless unique == nil
    end
end

def clean_profile(profile, hpo, options)
	cleaned_profile = hpo.clean_profile_hard(profile)
	unless options[:term_filter].nil?
		cleaned_profile.select! {|hp| hpo.get_ancestors(hp).include?(options[:term_filter])}
	end	
	return cleaned_profile
end

def clean_profiles(profiles, hpo, options)
	removed_profiles = []
	profiles.each do |id, phen_profile|
		cleaned_profile = clean_profile(phen_profile, hpo, options)
		profiles[id] = cleaned_profile
		removed_profiles << id if cleaned_profile.empty?
	end
	removed_profiles.each{|rp| profiles.delete(rp)}
	return removed_profiles
end

def clean_PACO(paco, hpo, options)
	removed_profiles = []
	idx = []
	count = 0
	paco.map! do |record|
		cleaned_profile = clean_profile(record.last, hpo, options)
		if cleaned_profile.empty?
			removed_profiles << record.first
			idx << count
		else
			record[-1] = cleaned_profile
		end	
		count += 1
		record
	end
	idx.reverse.each{|i| paco.delete_at(i)}
	return removed_profiles		
end	

####################################################################################
## OPTPARSE
####################################################################################
options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{File.basename(__FILE__)} [options]"

  options[:input_file] = nil
  opts.on("-i", "--input_file PATH", "Filepath of profile data") do |item|
    options[:input_file] = item
  end

  options[:ontology_file] = nil
  opts.on("-o", "--ontology_file PATH", "Path to ontology file") do |item|
  	options[:ontology_file] = item
  end

  options[:term_filter] = nil
  opts.on("-t", "--term_filter STRING", "If specified, only terms that are descendants of the specified term will be kept on a profile") do |item|
  	options[:term_filter] = item.to_sym
  end
  
  options[:file_output] = nil
  opts.on("-f", "--file_output PATH", "Output filepath") do |item|
  	options[:file_output] = item
  end

  options[:pairs_clean] = false
  opts.on("-p", "--pairs_clean", "Remove ancestor-child hp pairs from netanalyzer output") do
  	options[:pairs_clean] = true
  end

  options[:removed_path] = nil
  opts.on("-r", "--removed_path PATH", "Desired path to write removed profiles file") do |item|
  	options[:removed_path] = item
  end 

  options[:PACO_clean] = false
  opts.on("-P", "--PACO_clean", "Clean PACO files") do
  	options[:PACO_clean] = true
  end			

end.parse!

if options[:removed_path].nil?
	options[:removed_path] = File.dirname(options[:file_output])
end	
####################################################################################
## MAIN
####################################################################################
hpo = Ontology.new(file: options[:ontology_file], load_file: true)

if options[:pairs_clean]
	output_type = 'profile'
	data = load_tabular_file(options[:input_file])
	data.map!{|pair| [pair[0].to_sym, pair[1].to_sym, pair[2]]}
	data.select! do |pair|
		term_a, term_b, score = pair
		!hpo.get_ancestors(term_a).include?(term_b) && !hpo.get_descendants(term_a).include?(term_b)
	end
elsif options[:PACO_clean]
	output_type = 'profile'
	data = load_tabular_file(options[:input_file])
	header = data.shift
	data.map! do |values| 
		profile = values.last.split("|").map {|x| x.to_sym}
		values[-1] = profile
		values
	end
	removed_profiles = clean_PACO(data, hpo, options)
	data.map! do |rec|
		rec[-1] = rec.last.map{|term| term.to_s}.join("|")
		rec
	end
	data.unshift(header)
else
	output_type = "pairs"	
	pairs = load_tabular_file(options[:input_file])
	data = pairs2profiles(pairs)
	data.transform_values!{|values| values.map{|x| x.to_sym}}
	removed_profiles = clean_profiles(data, hpo, options)
end

File.open(options[:file_output], 'w') do |file|
	if output_type == "pairs"
		data.each do |record|
			record.last.each do |hp|
				file.puts([record.first, hp].join("\t"))
			end
		end
	elsif output_type == "profile"
		data.each do |record|
			file.puts record.join("\t")
		end
	end		
end


if !removed_profiles.nil? && !removed_profiles.empty?
	rejected_file = File.basename(options[:input_file], ".*")+'_excluded_patients'
	file = File.join(options[:removed_path], rejected_file)
	File.open(file, 'w') do |f|
		removed_profiles.each do |profile|
			f.puts profile
		end
	end
end	