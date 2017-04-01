module Enumerable
  def maximums_by(&block)
    el_vals = map{|el| [el, block.call(el)]}
    max = el_vals.map{|el, val| val }.max
    el_vals.select{|el, val| val == max }.map{|el, val| el }
  end
end

UniprotInfo = Struct.new(:ac, :id, :status, :protein_names, :gene_names, :annotation_score)
class UniprotInfo
  def initialize(ac, id, status, protein_names, gene_names, annotation_score)
    self.ac = ac
    self.id = id
    self.status = status
    self.protein_names = protein_names
    self.gene_names = gene_names
    self.annotation_score = annotation_score
  end

  def self.from_string(str)
    ac, id, status, protein_names, gene_names, annotation_score = str.chomp.split("\t")
    self.new(ac, id, status.to_sym, protein_names, gene_names.split, annotation_score.split.first.to_i)
  end
end

infos = File.readlines("uniprot_gene_names.tsv").drop(1).map{|line|
  UniprotInfo.from_string(line)
}.select{|ui| ui.status == :reviewed }

ui_by_genename = infos.flat_map{|ui|
  ui.gene_names.each_with_index.map{|gene_name, ind|
    [gene_name, ui]
  }
}.group_by{|gene_name, ui|
  gene_name
}.map{|gene_name, name_uis|
  uis = name_uis.map{|name, ui| ui }
  [gene_name, uis]
}.map{|gene_name, uis|
  if uis.size == 1
    [gene_name, uis[0]]
  else
    # uis.count{|ui|  }
    main_name_uis = uis.select{|ui|
      ui.gene_names.index(gene_name) == 0
    }
    if main_name_uis.size == 1
      [gene_name, main_name_uis[0]]
    else
      best_uis = uis.maximums_by(&:annotation_score)
      if best_uis.size == 1
        [gene_name, best_uis[0]]
      else
        nil # ignore
      end
    end
  end
}.compact.to_h

raise 'Specify filename with tissue specificities'  unless filename = ARGV[0]

lines = File.readlines(filename).map(&:chomp)

puts lines.first

lines.drop(1).each{|line|
  gene_name = line.split("\t").last
  uniprot_id = ui_by_genename[gene_name].id rescue ""
  puts line + "\t" + uniprot_id
}
