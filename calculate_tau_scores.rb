def tau_score(expressions)
  max_expression = expressions.max.to_f
  return 0  if max_expression == 0
  normed_expressions = expressions.map{|expr| expr / max_expression }
  n = expressions.size
  normed_expressions.map{|normed_expr| 1 - normed_expr }.inject(0.0, &:+) / (n - 1)
end

def hg_score(expressions)
  sum_expr = expressions.inject(0.0, &:+)
  return 0  if sum_expr == 0
  n = expressions.length
  normed_expr = expressions.map{|expr| expr / sum_expr }
  entropy = normed_expr.map{|xi| xi == 0 ? 0 : Math.log2(xi) * xi  }.inject(&:+)
  (Math.log2(n) + entropy) / Math.log2(n)
end

File.readlines('Chipseq_enrichment/MATRIX_rpkm.txt').drop(1).each{|line|
  gene, *expressions = line.chomp.split("\t")
  expressions = expressions.map(&:to_f)
  next  if expressions.all?(&:zero?)
  expressions = expressions.map{|x| x + 1.0 } ## pseudocount addition
  lognormed_expressions = expressions.map{|x| x <= 1 ? 0.0 : Math.log2(x) }
  puts [gene, tau_score(expressions), tau_score(lognormed_expressions), hg_score(expressions), hg_score(lognormed_expressions)].join("\t")
}
