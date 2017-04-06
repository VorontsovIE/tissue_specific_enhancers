require 'set'
require_relative 'statistics/fisher_table'

raise 'Specify combined chipseq-validation results for a tissue' unless tissue_chipseq_fn = ARGV[0]
raise 'Specify tissue-specific uniprots' unless tissue_specific_fn = ARGV[1]

tissue_specific_uniprots = File.readlines(tissue_specific_fn).map(&:chomp).to_set

tf_infos = File.readlines(tissue_chipseq_fn).drop(1).map{|l|
  l.split("\t")
}

all_tfs = tf_infos.map{|tf,*rest| "#{tf}_MOUSE" }.to_set

uniprots_of_interest = tf_infos.map{|row|
  tf = row[0]
  rate_1 = row[1].to_f
  significance_1 = row[2].to_f
  rate_2 = row[7].to_f
  significance_2 = row[8].to_f
  [tf, [rate_1, significance_1], [rate_2, significance_2]]
}.select{|tf,(r1,s1),(r2,s2)|
  alpha = 0.05 / tf_infos.size # Bonferroni correction
  r1 > 1 && s1 < alpha || r2 > 1 && s2 < alpha
}.map{|tf,*rest|
  "#{tf}_MOUSE"
}.to_set

# p tissue_specific_uniprots
# p uniprots_of_interest

fisher_table = FisherTable.by_two_classes(
  class_a_positive: (tissue_specific_uniprots & uniprots_of_interest).size,
  class_a_negative: (tissue_specific_uniprots - uniprots_of_interest).size,
  class_b_positive: ((all_tfs - tissue_specific_uniprots) & uniprots_of_interest).size,
  class_b_negative: ((all_tfs - tissue_specific_uniprots) - uniprots_of_interest).size
)

# tissue = File.basename(tissue_chipseq_fn, File.extname(tissue_chipseq_fn))
puts "#{tissue_chipseq_fn}\tRate: #{fisher_table.a_to_b_positive_rate_ratio.round(3)}\tSignificance: #{fisher_table.significance}\tFisher table: #{fisher_table}"
