def float(x)
  x == 'inf' ? Float::INFINITY : Float(x)
end
TfTissueInfo = Struct.new(:tf, :cl, :logpvalue, :se_te_means_ratio, :se_te_medians_ratio, :num_te, :num_se, :te_mean, :se_mean, :te_median, :se_median, keyword_init: true) do
  def self.from_string(str)
    tf, cl, logpvalue, se_te_means_ratio, se_te_medians_ratio, num_te, num_se, te_mean, se_mean, te_median, se_median = str.chomp.split("\t")
    self.new(tf: tf, cl: cl, logpvalue: float(logpvalue),
             se_te_means_ratio: float(se_te_means_ratio), se_te_medians_ratio: float(se_te_medians_ratio),
             num_te: Integer(num_te), num_se: Integer(num_se),
             te_mean: Float(te_mean), se_mean: Float(se_mean), te_median: Float(te_median), se_median: Float(se_median))
  end
  def self.each_in_file(fn, &block)
    return enum_for(:each_in_file, fn)  unless block_given?
    File.open(fn){|f|
      f.readline # drop header
      f.each_line{|l| yield self.from_string(l) }
    }
  end
end

header = ['-log10(P-value) bin', 'TF-tissue pairs with dense TE-s', 'TF-tissue pairs with dense SE-s', 'SE/TE dense tissue pairs ratio']
puts header.join("\t")
infos = TfTissueInfo.each_in_file('u-mann-whitney.tsv')
((0..3).step(0.5).to_a + [Float::INFINITY]).each_cons(2).map{|logpval_from, logpval_to|
  infos_in_bin = infos.select{|info| logpval_from <= info.logpvalue && info.logpvalue < logpval_to }
  num_te_win = infos_in_bin.count{|info| info.se_te_medians_ratio < 1 }
  num_se_win = infos_in_bin.count{|info| info.se_te_medians_ratio > 1 }
  row = ["#{logpval_from} - #{logpval_to}", num_te_win, num_se_win, num_se_win.to_f / num_te_win]
  puts row.join("\t")
}
