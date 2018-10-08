require_relative 'homotypic_score'

pvalue_cutoff = Float(ARGV[0])

header = [
  'chr', 'from', 'to',
  'MaxScore', 'NumSites',
  'Cluster length', 'Cluster lambda', 'Cluster homotypic score',
  'Interval length', 'Interval lambda', 'Interval homotypic score',
  # 'Corrected cluster homotypic score',
]
puts header.join("\t")

$stdin.readlines.map{|l|
  l.chomp.split("\t")
}.map{|chr, from,to, starts,ends,pvalues|
  starts = starts.split(",").map(&:to_i)
  ends = ends.split(",").map(&:to_i)
  pvalues = pvalues.split(",").map(&:to_f)
  max_pvalue = pvalues.max

#  cluster_length = ends.max - starts.min + 1
  cluster_length = starts.max - starts.min + 1
  len = (to.to_i - from.to_i)
  lambda_coeff = cluster_length * pvalue_cutoff
  cluster_lambda_coeff = len * pvalue_cutoff

# multiple testing gives larger probability to meet homotypic cluster of a given size
#  correction = Math.log10(len - cluster_length + 1)

  if max_pvalue == -1
    infos = [
      chr, from, to, 
      0, 0, 
      cluster_length, cluster_lambda_coeff, 0,
      len, lambda_coeff, 0,
      # -correction,
    ]
  else
    infos = [
      chr, from, to, max_pvalue, pvalues.size, 
      cluster_length, cluster_lambda_coeff, homotypic_score(cluster_lambda_coeff, pvalues.size),
      len, lambda_coeff, homotypic_score(lambda_coeff, pvalues.size),
      # homotypic_score(cluster_lambda_coeff, pvalues.size) - correction,
    ]
  end
  infos
}.each{|infos|
  puts infos.join("\t")
}
