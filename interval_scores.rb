require_relative 'homotypic_score'

pvalue_cutoff = Float(ARGV[0])

puts ['chr', 'from', 'to', 'MaxScore', 'NumSites', 'Len', 'Lambda', 'HomotypicScore',].join("\t")

$stdin.readlines.map{|l|
  l.chomp.split("\t")
}.map{|chr, from,to, pvalues|
  len = (to.to_i - from.to_i)
  lambda_coeff = len * pvalue_cutoff
  pvalues = pvalues.split(",").map(&:to_f)
  max = pvalues.max

  next [chr, from, to, 0, 0, len, lambda_coeff, 0]  if max == -1
  [chr, from, to, max, pvalues.size, len, lambda_coeff, homotypic_score(lambda_coeff, pvalues.size)]
}.each{|infos|
  puts infos.join("\t")
}
