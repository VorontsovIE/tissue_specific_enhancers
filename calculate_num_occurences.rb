$stdin.each_line{|l|
  *interval, occurences = l.chomp.split("\t")
  num_occurences = occurences.split(',').reject{|x| x == '-1' }.count
  puts [*interval, num_occurences].join("\t")
}

