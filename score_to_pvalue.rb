require 'optparse'

def read_bsearch_table(filename)
  File.readlines(filename).map{|l| l.chomp.split("\t").map(&:to_f) }
end

def pvalue_by_score(requested_score, bsearch_table)
  (bsearch_table.bsearch{|score, pvalue| score >= requested_score } || bsearch_table.last).last
end

precision = nil
OptionParser.new{|opts|
  opts.on('--precision VALUE'){|val|
    precision = Integer(val)
  }
}.parse!(ARGV)

threshold_fn = ARGV.shift
threshold_pvalue_table = read_bsearch_table(threshold_fn)

$stdin.each_line{|line|
  chr, from, to, score = line.chomp.split("\t")
  score = score.to_f
  pvalue = pvalue_by_score(score, threshold_pvalue_table)
  logP = - Math.log(pvalue)
  logP_formatted = precision ? logP.round(precision) : logP
  puts [chr, from, to, logP_formatted].join("\t")
}
