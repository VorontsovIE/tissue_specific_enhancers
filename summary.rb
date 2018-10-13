def int(x)
  x == '--' ? 0 : Integer(x)
end

def float(x)
  x == '--' ? 0.0 : Float(x)
end

TFRow = Struct.new(:tf,
    :odds_ratio, :significance, :crm_chip, :crm_nochip, :control_chip, :control_nochip,
    :confirmed_odds_ratio, :confirmed_significance, :confirmed_crm_chip, :confirmed_crm_nochip, :confirmed_control_chip, :confirmed_control_nochip,
    :gene, :geneid, :rpkm, :genename, :taufrac, :tau) do
  
  def self.from_string(str)
    tf, \
    odds_ratio, significance, crm_chip, crm_nochip, control_chip, control_nochip, \
    confirmed_odds_ratio, confirmed_significance, confirmed_crm_chip, confirmed_crm_nochip, confirmed_control_chip, confirmed_control_nochip, \
    gene, geneid, rpkm, genename, taufrac, tau = str.chomp.split("\t")
    self.new(tf, \
      *[odds_ratio, significance].map{|x| float(x) }, \
      *[crm_chip, crm_nochip, control_chip, control_nochip].map{|x| int(x) }, \
      *[confirmed_odds_ratio, confirmed_significance].map{|x| float(x) }, \
      *[confirmed_crm_chip, confirmed_crm_nochip,confirmed_control_chip, confirmed_control_nochip].map{|x| int(x) }, \
      gene, geneid, float(rpkm), genename, float(taufrac), float(tau))
  end
  
  def self.each_in_file(fn, &block)
    return enum_for(:each_in_file, fn)  unless block_given?
    File.open(fn){|f|
      f.readline
      f.each_line{|l| yield self.from_string(l) }
    }
  end

  def to_s
    [tf, \
    odds_ratio, significance, crm_chip, crm_nochip, control_chip, control_nochip, \
    confirmed_odds_ratio, confirmed_significance, confirmed_crm_chip, confirmed_crm_nochip, confirmed_control_chip, confirmed_control_nochip, \
    gene, geneid, rpkm, genename, taufrac, tau].join("\t")
  end

def self.header_compact
    ['TF', 
    'Odds ratio', 'Significance', 'Corrected Significance',
    'Odds ratio.1', 'Significance.1', 'Corrected Significance.1',
    'Gene', 'GeneId', 'RPKM', 'GeneName', 'TauFraction', 'Tau'].join("\t")
end
def to_s_compact(bonferroni)
    [tf, \
    odds_ratio, significance, significance * bonferroni,
    confirmed_odds_ratio, confirmed_significance, confirmed_significance * bonferroni,
    gene, geneid, rpkm, genename, taufrac, tau].join("\t")
  end

  def total
    crm_chip + crm_nochip + control_chip + control_nochip
  end
end

folder = 'prefinal_08sep2018/TSRE/results_combined/SuperEnhancers_27ac/Constituent_Enhs/'
bonferroni = Dir.glob("#{folder}/*.tsv").map{|fn|
  TFRow.each_in_file(fn).reject{|infos|
    ['CDX4_MOUSE', 'EVX1_MOUSE', 'EVX2_MOUSE'].include? infos.tf
  }.count
}.inject(0, &:+)
# $stderr.puts bonferroni

puts ['Cell line', 'SE or TE', TFRow.header_compact].join("\t")
Dir.glob("#{folder}/*.tsv").sort.each{|fn|
  cl_enhancer = File.basename(fn, '.tsv')
  cl, enhancer_type, *_rest = cl_enhancer.split('.')
  rows = TFRow.each_in_file(fn).reject{|infos|
    infos.total == 0
  }.select{|infos|
    # infos.significance * bonferroni <= 1e-2
    true
  }.reject{|infos|
    ['CDX4_MOUSE', 'EVX1_MOUSE', 'EVX2_MOUSE'].include? infos.tf
  }.sort_by{|infos|
    # infos.taufrac
    # infos.significance
    infos.tau == 0 ? 0 : (infos.taufrac / infos.tau)
  }.reverse

  top_rows = rows.first(5) + rows.drop(5).take_while{|infos| (infos.tau == 0 ? 0 : (infos.taufrac / infos.tau)) == 1 }
  top_rows.each{|infos|
    puts cl + "\t" + enhancer_type.upcase + "\t" + infos.to_s_compact(bonferroni)
  }
  # .sort_by{|infos| infos[:odds_ratio] }.reverse
  # .sort_by{|infos| infos[:significance] }
  # .first(10).map(&:tf)
  # puts [cl, tfs.join(',')]
  # puts
}
