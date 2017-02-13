puts ARGF.map{|l| s,f = l.split("\t")[1,2].map(&:to_i); f-s }.sum
