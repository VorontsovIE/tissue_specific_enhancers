$stdin.lazy.reject{|l|
  l.start_with?('#')
}.map{|l|
  chr, source, type, from, to, score, strand, phase, attribs = l.chomp.split("\t")
  attribs = attribs.split(';').map{|attr| attr.split('=') }.to_h
  {chr: chr, from: from, to: to, strand: strand, type: type, gene_name: attribs['gene_name'], transcript_name: attribs['transcript_name']}
}.select{|data|
  data[:type] == 'transcript'
}.each{|data|
  puts data.values_at(:chr, :from, :to, :gene_name, :transcript_name, :strand).join("\t")
}
