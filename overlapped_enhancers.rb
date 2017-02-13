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

raise 'Specify cell-line filename'  unless cell_line_enhancers_fn = ARGV[0] # 'merged_mm10_TS_strong_enhancers/Liver.bed'
raise 'Specify shift factor' unless shift_factor = ARGV[1]
shift_factor = Integer(shift_factor)

puts ['TF', 'Rate ratio', 'Significance', 'Peaks covered', 'Peaks uncovered', 'Control covered', 'Control uncovered'].join("\t")
# tf_chipseq_fn = 'gtrd/adaptive_quality/ANDR_MOUSE.bed'
Dir.glob('gtrd/adaptive_quality/*_MOUSE.bed').sort.map{|tf_chipseq_fn|
  fisher_table = FisherTable.by_two_classes(
    class_a_positive: `#{peaks_cmd(cell_line_enhancers_fn) + '|' + num_peaks_overlapped_cmd(tf_chipseq_fn)}`.to_i,
    class_a_negative: `#{peaks_cmd(cell_line_enhancers_fn) + '|' + num_peaks_nonoverlapped_cmd(tf_chipseq_fn)}`.to_i,
    class_b_positive: `#{peaks_shifted_left_cmd(cell_line_enhancers_fn, shift_factor: shift_factor) + '|' + num_peaks_overlapped_cmd(tf_chipseq_fn)}`.to_i,
    class_b_negative: `#{peaks_shifted_left_cmd(cell_line_enhancers_fn, shift_factor: shift_factor) + '|' + num_peaks_nonoverlapped_cmd(tf_chipseq_fn)}`.to_i
  )

  tf = File.basename(tf_chipseq_fn, '_MOUSE.bed')
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
  puts [tf, fisher_table.a_to_b_positive_rate_ratio, fisher_table.significance.to_f, *fisher_values].join("\t")
}
