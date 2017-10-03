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

def join_tables(table_1, table_2, column_1: 0, column_2: 0, has_header_1: true, has_header_2: true, header_1: nil, header_2: nil, replacement_for_skip: false)
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
      if replacement_for_skip
        $stderr.puts "Replace #{row} with blank cells: `#{row[column_1]}` not found."
        row + replacement_for_skip
      else
        $stderr.puts "Skip #{row}: `#{row[column_1]}` not found."
        nil
      end
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

# only the first part (before dot) of specified filename will be treated as a tissue name
raise 'Specify filename starting with tissue name'  unless tissue = ARGV[0]
tissue = tissue.split(".")[0]

uniprot_genes_fn = 'uniprot_main_gene_name.tsv'
annotation_fn = 'TSRE/Tau_expression/annotation.txt'
rpkms_fn = 'TSRE/Tau_expression/MATRIX_rpkm.txt'
tau_fn = 'TSRE/Tau_expression/tau_fraction_max_matrix.txt'

uniprot_genename = table_from_file(uniprot_genes_fn)
overrepr = table_from_stream($stdin)

result = join_tables(overrepr, uniprot_genename, has_header_2: false, header_2: ['Gene'])

annotation_table = [['GeneId', 'Gene', 'Length']] + table_from_file(annotation_fn)
annotation_table = take_columns(annotation_table, detect_columns(annotation_table, ['Gene', 'GeneId']))

rpkm_table = table_from_file(rpkms_fn)
rpkm_table = take_columns(rpkm_table, detect_columns(rpkm_table, ['Gene', tissue]))
rpkm_table.shift
rpkm_table.unshift(['Gene', 'RPKM'])

tau_table = table_from_file(tau_fn)
tau_table = [tau_table.first + ['Tau']] + tau_table.drop(1).map{|row| tau = row.drop(2).map{|x| Float(x) }.max; [*row, tau] }
tau_table = take_columns(tau_table, detect_columns(tau_table, ['GeneName', 'GeneId', tissue, 'Tau']))
tau_table.shift
tau_table.unshift(['GeneName', 'GeneId', 'TauFraction', 'Tau'])

result = join_tables(result, annotation_table, column_1: detect_column(result, 'Gene'), column_2: detect_column(annotation_table, 'Gene'))
result = join_tables(result, rpkm_table, column_1: detect_column(result, 'GeneId'), column_2: detect_column(rpkm_table, 'Gene'))
result = join_tables(result, tau_table, column_1: detect_column(result, 'GeneId'), column_2: detect_column(tau_table, 'GeneId'), replacement_for_skip: ['--', '--'])
print_table(result)
