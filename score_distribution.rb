# require 'shellwords'
# require 'tempfile'

def peaks_cmd(filename)
  "cat #{filename}" + ' | ' + fill_gaps_cmd
end

def fill_gaps_cmd(pattern_size: 200, genome_sizes_fn: 'genome_sizes/mm10.genome')
  [
    "bedtools slop -g #{genome_sizes_fn} -b #{pattern_size}", # dilute
    "bedtools sort",
    "bedtools merge",
    "bedtools slop -g #{genome_sizes_fn} -b #{-pattern_size}", # erode
  ].join(' | ')
end


raise 'Specify cell-line filename'  unless cell_line_enhancers_fn = ARGV[0] # 'merged_mm10_TS_strong_enhancers/Liver.bed'
raise 'Specify chipseq file'  unless chipseq_fn = ARGV[1]
raise 'Specify motif' unless motif = ARGV[2]
raise 'Specify pvalue cutoff'  unless pvalue_cutoff = ARGV[3]
raise 'Specify enhancer type'  unless enh_type = ARGV[4]

raise 'Specify output bed file for enhancers-on-ChIPseq'  unless chipseq_enhancers_fn = ARGV[5]
# chipseq_enhancers_file = Tempfile.new
# chipseq_enhancers_file.close
# chipseq_enhancers_fn = chipseq_enhancers_file.path

raise 'Specify output file with raw logpvalues'  unless raw_pvals_fn = ARGV[6]
raise 'Specify output file with pvalues stats'  unless output_fn = ARGV[7]

pvalue_cutoff = Float(pvalue_cutoff)
motif_fn = "motif_collection/#{motif}.pwm"
thresholds_fn = "motif_thresholds/#{motif}.thr"

cmd_1 = [
  peaks_cmd(cell_line_enhancers_fn),
  "bedtools intersect -u -a - -b #{chipseq_fn}",
].join(' | ') + " > #{chipseq_enhancers_fn}"

cmd_2 = [
  "bedtools getfasta -fi genomes/mouse/mm10.fa -bed #{chipseq_enhancers_fn}",
  "java -cp sarus.jar ru.autosome.SARUS - #{motif_fn} #{pvalue_cutoff} --pvalues-file #{thresholds_fn} --output-scoring-mode logpvalue --suppress --threshold-mode pvalue --output-bed",
  "bedtools intersect -loj -a #{chipseq_enhancers_fn} -b -",
  "bedtools groupby -c 8 -o collapse",
].join(' | ') + " > #{raw_pvals_fn}"

cmd_3 = "cat #{raw_pvals_fn} | ruby interval_scores.rb #{pvalue_cutoff} | awk -e '{print \$0 \"\\t#{enh_type}\"}' > #{output_fn}"

puts "#{cmd_1}; #{cmd_2}; #{cmd_3}"
# system(cmd_1)
# system(cmd_2)
# system(cmd_3)

# chipseq_enhancers_file.unlink

#  'R --slave -e \'x <- read.csv(file="stdin", sep="\t"); summary(x)\'',
