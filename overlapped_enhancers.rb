require_relative 'statistics/fisher_table'

def fill_gaps_cmd(pattern_size: 200, genome_sizes_fn: 'genome_sizes/mm10.genome')
  [
    "bedtools slop -g #{genome_sizes_fn} -b #{pattern_size}", # dilute
    "bedtools merge",
    "bedtools slop -g #{genome_sizes_fn} -b #{-pattern_size}", # erode
  ].join('|')
end

def peaks_cmd(filename)
  "cat #{filename}" + '|' + fill_gaps_cmd
end

#ToDo: size_factor do not work
# Peak of size `size_factor` times larger than original, shifted to `shift_factor` times size of peak
def peaks_shifted_left_cmd(filename, size_factor: 1, shift_factor: 20, genome_sizes_fn: 'genome_sizes/mm10.genome')
  [
    peaks_cmd(filename),
    "bedtools slop  -g #{genome_sizes_fn} -pct -l #{shift_factor - 1}   -r 0",
    "bedtools flank -g #{genome_sizes_fn} -pct -l #{1.0 / shift_factor} -r 0",
  ].join('|')
end

# Peak of size `size_factor` times larger than original, shifted to `shift_factor` times size of peak
def peaks_shifted_right_cmd(filename, size_factor: 1, shift_factor: 20, genome_sizes_fn: 'genome_sizes/mm10.genome')
  [
    peaks_cmd(filename),
    "bedtools slop  -g #{genome_sizes_fn} -pct -l 0 -r #{shift_factor - 1}",
    "bedtools flank -g #{genome_sizes_fn} -pct -l 0 -r #{1.0 / shift_factor}",
    ].join('|')
end

# Peak of size `size_factor` times larger than original, shifted to `shift_factor` times size of peak
def peaks_shifted_leftright_cmd(filename, size_factor: 1, shift_factor: 20, genome_sizes_fn: 'genome_sizes/mm10.genome')
  doubled_size = 2.0 * shift_factor - 1.0
  [
    peaks_cmd(filename),
    "bedtools slop  -g #{genome_sizes_fn} -pct -l #{shift_factor - 1} -r #{shift_factor - 1}",
    "bedtools flank -g #{genome_sizes_fn} -pct -l #{1.0 / doubled_size} -r #{1.0 / doubled_size}",
    ].join('|')
end

def num_peaks_overlapped_cmd(overlapped_by_fn)
  [
    "bedtools intersect -a - -b #{overlapped_by_fn} -u",
    "wc -l",
  ].join('|')
end

def num_peaks_nonoverlapped_cmd(overlapped_by_fn)
  [
    "bedtools intersect -a - -b #{overlapped_by_fn} -v",
    "wc -l",
  ].join('|')
end

###################
raise 'Specify cell-line filename'  unless cell_line_enhancers_fn = ARGV[0] # 'merged_mm10_TS_strong_enhancers/Liver.bed'

raise 'Specify shift factor' unless shift_factor = ARGV[1]
shift_factor = Integer(shift_factor)

raise 'Specify folder with bed files of DNA-regions bound by different TFs'  unless bound_regions_folder = ARGV[2]
# bound_regions_folder = 'gtrd/adaptive_quality/'
# bound_regions_folder = 'gtrd/confirmed_by_motif/'
##################

cmd = peaks_shifted_leftright_cmd(cell_line_enhancers_fn, shift_factor: shift_factor) + ' | ' + \
                                  num_peaks_overlapped_cmd(" <( #{peaks_cmd(cell_line_enhancers_fn)} )")
puts "echo -n #{File.basename(cell_line_enhancers_fn)} \" \"; #{peaks_cmd(cell_line_enhancers_fn)} | wc -l | xargs echo -n ; echo -n \" \"; #{cmd}"
#puts [`#{peaks_cmd(cell_line_enhancers_fn)} | wc -l`, `#{cmd}`].join("\t")

=begin

headers =  [
  'TF', 'Odds ratio', 'Significance',
  'CRM covered by ChIPSeq peak', 'CRM uncovered by ChIPSeq peak',
  'Control covered by ChIPSeq peak', 'Control uncovered by ChIPSeq peak'
]

puts headers.join("\t")


Dir.glob(File.join(bound_regions_folder, '*_MOUSE.bed')).reject{|tf_chipseq_fn|
  File.size(tf_chipseq_fn) == 0
}.sort.map{|tf_chipseq_fn|
  fisher_table = FisherTable.by_two_classes(
    class_a_positive: `#{peaks_cmd(cell_line_enhancers_fn) + '|' + num_peaks_overlapped_cmd(tf_chipseq_fn)}`.to_i,
    class_a_negative: `#{peaks_cmd(cell_line_enhancers_fn) + '|' + num_peaks_nonoverlapped_cmd(tf_chipseq_fn)}`.to_i,
    class_b_positive: `#{peaks_shifted_leftright_cmd(cell_line_enhancers_fn, shift_factor: shift_factor) + '|' + num_peaks_overlapped_cmd(tf_chipseq_fn)}`.to_i,
    class_b_negative: `#{peaks_shifted_leftright_cmd(cell_line_enhancers_fn, shift_factor: shift_factor) + '|' + num_peaks_nonoverlapped_cmd(tf_chipseq_fn)}`.to_i
  )

  tf = File.basename(tf_chipseq_fn, '.bed')
  [tf, fisher_table]
}.sort_by{|tf, fisher_table|
  fisher_table.a_to_b_positive_rate_ratio || -1
}.reverse.each{|tf, fisher_table|
  fisher_values = [
    fisher_table.class_a_positive,
    fisher_table.class_a_negative,
    fisher_table.class_b_positive,
    fisher_table.class_b_negative
  ]
  puts [tf, fisher_table.a_to_b_positive_rate_ratio || '--', fisher_table.significance.to_f, *fisher_values].join("\t")
}

=end
