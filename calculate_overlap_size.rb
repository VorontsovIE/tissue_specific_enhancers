$stdin.each_line{|l|
  *interval, starts, finishes = l.chomp.split("\t")
  starts = starts.split(',').reject{|x| x == '-1' }.map(&:to_i)
  finishes = finishes.split(',').reject{|x| x == '-1' }.map(&:to_i)
  overlap_size = starts.zip(finishes).map{|start, finish| finish - start }.inject(0, &:+)
  puts [*interval, overlap_size].join("\t")
}
