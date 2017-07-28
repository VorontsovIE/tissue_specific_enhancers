fn_1, fn_2 = ARGV[0,2]
line_by_tf_1 = File.readlines(fn_1).drop(1).map(&:chomp).map{|line| [line.split("\t")[0], line] }.to_h
line_by_tf_2 = File.readlines(fn_2).drop(1).map(&:chomp).map{|line| [line.split("\t")[0], line] }.to_h

#puts(File.readlines(fn_1).first.chomp + "\t" + File.readlines(fn_2).first.chomp)
header_1 = File.readlines(fn_1).first.chomp
header_2 = File.readlines(fn_2).first.chomp.split("\t")
header_2 = (header_2.first(1) + header_2[1...-6].map{|el| "#{el}.1" } + header_2[-6,6]).join("\t")
puts(header_1 + "\t" + header_2)

line_by_tf_1.each_key{|tf|
  if line_by_tf_2.has_key?(tf)
    puts(line_by_tf_1[tf] + "\t" + line_by_tf_2[tf])
  else
    infos_1 = line_by_tf_1[tf].split("\t")
    puts line_by_tf_1[tf] + "\t" + (['--'] * (infos_1.size - 4) + infos_1[-4,4]).join("\t")
  end
}
