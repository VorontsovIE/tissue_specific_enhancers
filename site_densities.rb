require 'fileutils'
require 'statsample'

def densities(motif_coverages, overlap_sizes)
  motif_coverages.zip(overlap_sizes).each{|(interval_1, motif_coverage), (interval_2, overlap)|
    raise "Intervals `#{interval_1.join(':')}` and `#{interval_2.join(':')}` don't match"  unless interval_1 == interval_2
  }

  motif_coverages.zip(overlap_sizes).select{|(interval_1, motif_coverage), (interval_2, overlap)|
    overlap > 0
  }.map{|(interval_1, motif_coverage), (interval_2, overlap)|
    {motif_coverage: motif_coverage, cistrome_overlap: overlap, motif_density: motif_coverage.to_f / overlap}
  }
end

enhancers = Dir.glob("TSRE/glued/SuperEnhancers_27ac/Constituent_Enhs/*.bed").map{|glued_enhancers_fn|
  File.basename(glued_enhancers_fn, '.const.txt.glued.bed').sub(/\.[st]e$/, '')
}.uniq

header = ['TF', 'Cell line', '-log10(P-value)', 'SE/TE mean densities ratio', 'SE/TE median density ratio', '#TE', '#SE', 'TE mean', 'SE mean', 'TE median', 'SE median']
puts header.join("\t")

Dir.glob("motif_occurences_in_cistrome/*.bed").each do |tf_fn|
  tf = File.basename(tf_fn, '.bed')
  FileUtils.mkdir_p "densities/#{tf}"
  enhancers.each do |enhancer|
    se_te_infos = ['se', 'te'].map{|enhancer_type|
      motif_coverages = File.readlines("num_occurences_per_enhancer/#{tf}/#{tf}.#{enhancer}.#{enhancer_type}.const.txt.bed").map{|l| chr,from,to, coverage = l.chomp.split("\t"); [[chr,from,to], coverage.to_i] }
      overlap_sizes = File.readlines("cistrome_overlap_per_enhancer/#{tf}/#{tf}.#{enhancer}.#{enhancer_type}.const.txt.bed").map{|l| chr,from,to, overlap = l.chomp.split("\t"); [[chr,from,to], overlap.to_i] }
      infos = densities(motif_coverages, overlap_sizes)
      output_fn = "densities/#{tf}/#{tf}.#{enhancer}.#{enhancer_type}.const.txt.tsv"
      [infos, enhancer_type, output_fn]
    }
    se_te_infos = se_te_infos.each{|infos, enhancer_type, output_fn|
=begin
      File.open(output_fn, 'w') do |fw|
        fw.puts ['coverage', 'overlap', 'density'].join("\t")
        infos.each{|enhancer_info|
          fw.puts enhancer_info.values_at(:motif_coverage, :cistrome_overlap, :motif_density).join("\t")
        }
      end
=end
    }
    se_densities, te_densities = se_te_infos.map{|infos, enh_type, output_fn|
                                   infos
                                 }.map{|info|
                                   info.map{|enh_info|
                                     enh_info[:motif_density]
                                   }
                                 }
    begin
    se_densities.reject!(&:zero?)
    te_densities.reject!(&:zero?)
    if (se_densities.size > 2) && (te_densities.size > 2)
      pvalue = Statsample::Test.u_mannwhitney(
                 Daru::Vector.new(se_densities),
                 Daru::Vector.new(te_densities)
               ).probability_exact
#      pvalue = Statsample::Test.t_two_samples_independent(
#                 Daru::Vector.new(se_densities),
#                 Daru::Vector.new(te_densities)
#               ).probability_equal_variance
      te_mean = te_densities.inject(0.0, &:+) / te_densities.size
      se_mean = se_densities.inject(0.0, &:+) / se_densities.size
      te_median = te_densities.sort[te_densities.size / 2]
      se_median = se_densities.sort[se_densities.size / 2]
      infos = [
        tf, enhancer, -Math.log10(pvalue),
        (te_mean != 0) ? se_mean / te_mean : (se_mean == 0 ? 'N/A' : 'inf'),
        (te_median != 0) ? se_median / te_median : (se_median == 0 ? 'N/A' : 'inf'),
        te_densities.size, se_densities.size,
        te_mean, se_mean,
        te_median, se_median,
      ]
      puts infos.join("\t") #if pvalue <= 0.001
    end
    rescue
      $stderr.puts [tf, enhancer, se_densities.inspect, te_densities.inspect].join("\t")
    end
  end
end
