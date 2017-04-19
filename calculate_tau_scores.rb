File.readlines('Chipseq_enrichment/MATRIX_rpkm_pseudoCount.txt').drop(1).each{|line|
  gene, *expressions = line.chomp.split("\t")
  expressions.map!(&:to_f)
  max_expression = expressions.max
  normed_expressions = expressions.map{|expr| expr / max_expression }
  n = expressions.size
  tau = normed_expressions.map{|normed_expr| 1 - normed_expr }.inject(0.0, &:+) / (n - 1)
  puts "#{gene}\t#{tau}"
}
