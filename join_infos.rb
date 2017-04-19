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

raise 'Specify tissue'  unless tissue = ARGV[0]

rpkms = table_from_file('Chipseq_enrichment/MATRIX_rpkm_pseudoCount.txt')
uniprot_genename = table_from_file('uniprot_main_gene_name.tsv')
# tau_scores = table_from_file('tau_scores.txt')

overrepr = table_from_stream($stdin)

# p overrepr.first(2)
# p uniprot_genename.first(2)
# p rpkms.first(2)

result = join_tables(overrepr, uniprot_genename, has_header_2: false, header_2: ['Gene'])
tissue_expressions = take_columns(rpkms, detect_columns(rpkms, ['Gene', tissue]))
result = join_tables(result, tissue_expressions, column_1: detect_column(result, 'Gene'))
print_table(result)
