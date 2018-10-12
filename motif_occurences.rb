require 'fileutils'

def motifs_in_region_cmd(region_bed_fn:, motif_fn:, thresholds_fn:, pvalue_cutoff:)
  [
    "bedtools getfasta -bed #{region_bed_fn} -fi all_chromosomes.fa",
    "java -cp sarus-2.0.1.jar ru.autosome.SARUS - #{motif_fn} #{pvalue_cutoff} --pvalues-file #{thresholds_fn} --output-scoring-mode logpvalue --threshold-mode pvalue --output-bed --precision 2",
    "cut -f 1,2,3"
  ].join(' | ')
end

FileUtils.mkdir_p 'motif_occurences_in_cistrome'

tfs = Dir.glob('gtrd/adaptive_quality/*.bed').map{|fn| File.basename(fn, '.bed') }
tfs.each do |tf|
  motif_names = Dir.glob("motif_collection/#{tf}.*.pwm").map{|fn| File.basename(fn, '.pwm') }
  next  if motif_names.empty?
  cmds = motif_names.map{|motif|
    motifs_in_region_cmd(region_bed_fn: "gtrd/adaptive_quality/#{tf}.bed",
                         motif_fn: "motif_collection/#{motif}.pwm",
                         thresholds_fn: "motif_thresholds/#{motif}.thr",
                         pvalue_cutoff: 0.0005)
  }
  cmd = '( ' + cmds.join('; ') + ' ) | cut -f 1,2,3 | bedtools sort | bedtools merge ' + "> motif_occurences_in_cistrome/#{tf}.bed"
  puts cmd
end
