def float(x)
  if x == 'Infinity'
    Float::INFINITY
  elsif x == '-Infinity'
    -Float::INFINITY
  else
    Float(x)
  end
end

def median(xs)
  xs.sort[xs.length/2]
end

def mean(xs)
  xs.inject(0.0, &:+) / xs.size
end

def vals(fn, col_idx: 1)
  File.readlines(fn).drop(1).map{|l| float( l.chomp.split("\t")[col_idx] )}
end

def norm_numsites(fn)
  File.readlines(fn).drop(1).map{|l| num, lambda_coeff = l.chomp.split("\t").values_at(4,6).map{|x| float(x) }; num / lambda_coeff }
end

fold = 2.5
['SuperEnhancers/SuperEnhs/', 'SuperEnhancers/Constituent_Enhs/', 'SuperEnhancers_27ac/SuperEnhs/', 'SuperEnhancers_27ac/Constituent_Enhs/',].each do |folder|
  puts
  puts "-- #{folder} --"
  Dir.glob("score_stats/#{folder}/*.te.*.txt").map{|fn|
    te_basename = File.basename(fn)
    se_basename = te_basename.sub(/\.te\./, '.se.')
    # te_vals = vals("score_stats/#{folder}/#{te_basename}", col_idx: 4)
    # se_vals = vals("score_stats/#{folder}/#{se_basename}", col_idx: 4)
    te_vals = norm_numsites("score_stats/#{folder}/#{te_basename}")
    se_vals = norm_numsites("score_stats/#{folder}/#{se_basename}")
    te_num = te_vals.size
    se_num = se_vals.size

    median_homoscore_te = median(te_vals)
    median_homoscore_se = median(se_vals)
    if [te_num, se_num].min >= 20 && [median_homoscore_te, median_homoscore_se].min > 0 && ((median_homoscore_se / median_homoscore_te > fold) || (median_homoscore_te / median_homoscore_se > fold))
      puts [te_basename, "#TE: #{te_num}", "#SE: #{se_num}", median_homoscore_te, median_homoscore_se].join("\t")
    end
  }
end
