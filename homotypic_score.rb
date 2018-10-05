def fact(k)
  k == 0 ? 1 : fact(k-1) * k
end

# P(k >= n) -- probability to have no less than n occurences
def prob_more_occurs(lambda, n)
  #res = 1.0 - (0...n).map{|k| lambda**k / fact(k) }.inject(0.0, &:+) * Math.exp(-lambda)
  res = 1.0 - (0...n).map{|k|
    Math.exp(k * Math.log(lambda) - Math.lgamma(k + 1).first)
  }.inject(0.0, &:+) * Math.exp(-lambda)
  [res, 0].max
end

def homotypic_score(lambda, num_occurences)
  #$stderr.puts [lambda, num_occurences, prob_more_occurs(lambda, num_occurences), prob_more_occurs(lambda, 1), -Math.log10(prob_more_occurs(lambda, num_occurences) / prob_more_occurs(lambda, 1))].inspect  if prob_more_occurs(lambda, num_occurences) == 0
  return 0  if num_occurences == 1
  -Math.log10(prob_more_occurs(lambda, num_occurences) / prob_more_occurs(lambda, 1))
end
