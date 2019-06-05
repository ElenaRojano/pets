#! /usr/bin/env ruby

N_CUTS = 200

pp_limit = ARGV[1].to_f if !ARGV[1].nil?

count = 0
values = []
File.open(ARGV[0]).each do |line|
	line.chomp!
	count += 1
	next if count == 1
	fields = line.split("\t")
	values << [fields.first, fields.last.to_f]
end
values.sort!{|v1, v2| v1.last <=> v2.last}
max_score = values.last.last
min_score = values.first.last

interval = (max_score - min_score).fdiv(N_CUTS)
cuts = []
current = min_score
while cuts.length < N_CUTS
	cuts << current
	current += interval
end

header =  %w[tp tn fp fn cut pre rec]
header << 'group' if !ARGV[2].nil?
puts header.join("\t")
last_pre = 1
last_rec = 0
change_data = false
pre_range = 0
rec_range = 0
all_weigths = []
total_weigth = 0
cuts.reverse.each_with_index do |cut, i|
	tp = 0
	tn = 0
	fp = 0
	fn = 0
	values.each do |label, score|
		if score >= cut
			if label == 'in'
				tp += 1
			else
				fp += 1
			end
		else
			if label == 'out'
				tn += 1
			else
				fn += 1
			end
		end
	end
	pre = tp.fdiv(tp+fp)
	rec = tp.fdiv(tp+fn)
	if !ARGV[1].nil? && pp_limit > 0
		if !change_data
			pp = pre/last_pre 
			if pp >= pp_limit
				change_data = true
				#pre_range = last_pre/(N_CUTS - i)
				pre_range = last_pre#/(N_CUTS - i)
				rec_range = (1 - last_rec)/(N_CUTS - i)
				(N_CUTS - i).times do |n|
					all_weigths << (n+1)**8
					total_weigth += (n+1)**8
				end
			else
				last_pre = pre
				last_rec = rec
			end
		end
		if change_data
			tp = tn = fp = fn = 0
			last_pre -= pre_range * all_weigths.pop.fdiv(total_weigth)
			last_rec += rec_range
			pre = last_pre
			pre = 0 if pre < 0
			rec = last_rec
		end
	end
	row = [tp, tn, fp, fn, cut, pre, rec]
	row << ARGV[2] if !ARGV[2].nil?
	puts row.join("\t")
end

