def detect_columns(table, patterns)
  patterns.map{|pattern| detect_column(table, pattern) }
end

def detect_column(table, pattern)
  case pattern
  when Integer
    pattern
  else
    table[0].index{|cell, idx| pattern === cell }
  end
end

def take_column(table, column)
  table.map{|row| row[column] }
end

def take_columns(table, columns)
  table.map{|row| row.values_at(*columns) }
end

def without_element(list, idx_to_reject)
  list.each_with_index.reject{|el, idx| idx == idx_to_reject }.map{|el, idx| el }
end

def join_tables(table_1, table_2, column_1: 0, column_2: 0, has_header_1: true, has_header_2: true, header_1: nil, header_2: nil)
  if has_header_1
    header_1 = table_1.first
    table_1 = table_1.drop(1)
  end
  if has_header_2
    header_2 = without_element(table_2.first, column_2)
    table_2 = table_2.drop(1)
  end
  content = table_1.map{|row|
    row_to_join = table_2.detect{|row_2| row[column_1] == row_2[column_2] }
    if row_to_join
      row + without_element(row_to_join, column_2)
    else
      $stderr.puts "Skip #{row}: `#{row[column_1]}` not found."
      nil
    end
  }.compact
  
  header = []
  
  header += header_1 || ['--'] * table_1.first.length
  header += header_2 || ['--'] * (table_2.first.length - 1)

  (has_header_1 || has_header_2) ? [header] + content : content
end

def table_from_stream(stream)
  stream.readlines.map{|l| l.chomp.split("\t") }
end

def table_from_file(filename)
  File.readlines(filename).map{|l| l.chomp.split("\t") }
end

def print_table(table)
  table.each{|row| puts row.join("\t") }
end

def tau_score(expressions)
  max_expression = expressions.max.to_f
  return 0  if max_expression == 0
  normed_expressions = expressions.map{|expr| expr / max_expression }
  n = expressions.size
  normed_expressions.map{|normed_expr| 1 - normed_expr }.inject(0.0, &:+) / (n - 1)
end

raise 'Specify tissue'  unless tissue = ARGV[0]

rpkms = table_from_file('Chipseq_enrichment/MATRIX_rpkm_pseudoCount.txt')
uniprot_genename = table_from_file('uniprot_main_gene_name.tsv')

uniprot_genename.first(10)

overrepr = table_from_stream($stdin)

tissue_column = detect_column(rpkms, tissue) - 1 # first columns is "Gene"

result = join_tables(overrepr, uniprot_genename, has_header_2: false, header_2: ['Gene'])
genes = take_column(result, detect_column(result, 'Gene')).drop(1)

tau_scores = File.readlines('Chipseq_enrichment/quant.norm.csv').drop(1).map{|line|
  ensg, *logexpressions = line.chomp.split(",") # already log-normed
  logexpressions = logexpressions.map(&:to_f)
  expressions = logexpressions.map{|x| 2 ** x }
  expr = expressions[tissue_column].to_f 
  rel_expr = expr / expressions.inject(0.0, &:+)
  relmax_expr = expr / expressions.max
  tau = tau_score(logexpressions)
  [ensg, expr, tau, rel_expr * tau, relmax_expr * tau]
}

tau_scores_by_ensg = tau_scores.map{|ensg, expr, tau, rel_tau, relmax_tau| [ensg, [expr, tau, rel_tau, relmax_tau]] }.to_h

ensgs_by_gene = File.readlines('Chipseq_enrichment/ensembl_67.txt').drop(1).map{|l|
  l.chomp.split("\t").first(2)
}.group_by(&:last).map{|gene, ensg_gene_pairs|
  [gene, ensg_gene_pairs.map(&:first)]
}.to_h

tau_scores_by_gene = genes.map{|gene|
  ensgs = ensgs_by_gene[gene]
  tau_infos = ensgs.map{|ensg|
    tau_scores_by_ensg[ensg]
  }.compact
  if !tau_infos.empty?
    exprs, taus, rel_taus, relmax_taus = tau_infos.transpose
    [gene, ensgs.join(';'), exprs.max, taus.max, rel_taus.max, relmax_taus.max]
  else
    [gene, '', 0, 0, 0, 0]
  end
}

tau_table = [['Gene', 'Ensemble gene ID', 'Expression', 'Tau score', 'Tau fraction', 'Tau fraction from max']] + tau_scores_by_gene


#tissue_expressions = take_columns(rpkms, detect_columns(rpkms, ['Gene', tissue]))
#result = join_tables(result, tissue_expressions, column_1: detect_column(result, 'Gene'))

result = join_tables(result, tau_table, column_1: detect_column(result, 'Gene'))
print_table(result)
