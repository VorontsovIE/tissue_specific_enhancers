mouse_tfs = File.readlines('mouse_tfs.tsv').map{|l|
  l.chomp.split("\t",4)
}.map{|name, term_id, uniprot_acs, uniprot_id|
  [name, term_id, uniprot_acs.split(','), uniprot_id]
}

# mgi_by_uniprotac = File.readlines('MRK_Sequence.tsv').map{|l|
#   l.chomp.split("\t")
# }.map{|l|
#   l.values_at(0,1,14,13,17)
# }.drop(2).reject{|l| l[2].empty? }.flat_map{|l|
#   l[2].split('|').map{|ac| [ac, l] }
# }.to_h; mgi_by_uniprotac.size

# mouse_tfs.map{|name, termid, acs, id|
#   [name, termid, acs, id, acs.map{|ac| mgi_by_uniprotac[ac] }.compact.uniq]
# }.select{|name, termid, acs, id, mgis| ! mgis.empty?  }.first(3).each{|x| p x }; nil

uniprots = File.readlines('uniprot_mgi.tsv').drop(1).map{|l|
  l.chomp.split("\t")
}.map{|ac, id, name, sym, mgi, ensembl|
  [id, [ac, id, mgi, ensembl]]
}.to_h

entrez_by_mgi = File.readlines('MGI_EntrezGene.tsv').drop(1).map{|l|
  l.chomp.split("\t")
}.map{|l|
  [l[0], l[2], l[8]]
}.reject{|mgi, status, entrez|
  status == 'W' # withdrawn
}.map{|mgi, status, entrez|
  [mgi, entrez]
}.to_h

symbol_by_mgi = File.readlines('MGI_EntrezGene.tsv').drop(1).map{|l|
  l.chomp.split("\t")
}.map{|l|
  [l[0], l[2], l[1]]
}.reject{|mgi, status, symbol|
  status == 'W' # withdrawn
}.map{|mgi, status, symbol|
  [mgi, symbol]
}.to_h


# File.readlines('MRK_Sequence.tsv').map{|l|
#   l.chomp.split("\t")
# }.map{|l|
#   l #.values_at(0,1,14,13,17)
# }.select{|l|
#   l[1] && l[1].start_with?('Mef2b')
# }[0].each_with_index{|x, ind| p "#{ind}) #{x}" }; nil


# File.readlines('MRK_Sequence.tsv').first.chomp.split("\t").each_with_index{|x, ind| puts "#{ind}) #{x}" }

enst_to_ensg = File.readlines('enst_ensg.tsv').drop(1).map{|l|
  ensg, enst = l.chomp.split("\t")
  [enst, ensg]
}.to_h

final_infos = mouse_tfs.map{|name, termid, acs, id|
  [name, termid, acs, id, uniprots[id]]
}.reject{|name, termid, acs, id, infos|
  infos.nil?
}.map{|name, termid, acs, id, infos|
  _uniprot_ac, _uniprot_id, mgis, ensts = *infos
  mgis = (mgis || "").strip.split(';').reject(&:empty?).map{|mgi| "MGI:#{mgi}" }
  ensts = (ensts || "").strip.split(';').reject(&:empty?)
  ensgs = ensts.map{|enst| enst_to_ensg[enst] }.uniq
  entrezes = mgis.map{|mgi| entrez_by_mgi[mgi] }.compact.uniq
  symbols = mgis.map{|mgi| symbol_by_mgi[mgi] }.compact.uniq
  [name, symbols, termid, acs, id, mgis, ensgs, entrezes, ensts]
}

File.open('mouse_tf_infos.tsv', 'w'){|fw|
  fw.puts ['Name', 'Symbol', 'TFClass id', 'Uniprot ACs', 'Uniprot ID', 'MGIs', 'ENSGs', 'EntrezGenes', 'ENSTs'].join("\t")
  final_infos.each{|name, symbols, termid, acs, id, mgis, ensgs, entrezes, ensts|
    fw.puts [name, symbols.join(';'), termid, acs.join(';'), id, mgis.join(';'), ensgs.join(';'), entrezes.join(';'), ensts.join(';')].join("\t")
  }
}
