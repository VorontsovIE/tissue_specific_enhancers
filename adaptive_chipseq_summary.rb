require 'fileutils'
require 'shellwords'

qual_names = {'highest'=>'A','high'=>'B','medium'=>'C','low'=>'D', 'lowest' => 'E'}

def num_lines(fn)
  # cmd = ['wc', '-l', fn].shelljoin
  # `#{cmd}`.to_i
  File.foreach(fn).count
end

FileUtils.mkdir_p 'gtrd/adaptive_quality'
# take `quality_default_level` best datasets. If resulting file is empty, add dataset with the next quality level
quality_default_level = 2
quality_levels = ['highest', 'high', 'medium', 'low'] #, 'lowest'
tf_list = Dir.glob('gtrd/highest/*.bed').map{|fn| File.basename(fn, '.bed')[0..-3] }
tf_list.each do |tf|
  tf_chipseq_filelist = (quality_default_level .. quality_levels.size).map{|quality_level|
    quality_levels.first(quality_level).flat_map{|level|
      Dir.glob("gtrd/#{level}/#{tf}.?.bed")
    }
  }.detect{|filelist|
    filelist.map{|fn| num_lines(fn) }.inject(0, &:+) >= 100
  }
  tf_chipseq_filelist ||= quality_levels.flat_map{|level| Dir.glob("gtrd/#{level}/#{tf}.?.bed") }

  quals = tf_chipseq_filelist.map{|fn| File.basename(File.dirname(fn)) }
  peaks_fn = "gtrd/adaptive_quality/#{tf}.bed"
  peaks_confirmed_fn = "gtrd/confirmed_by_motif/#{tf}.bed"
  infos = [
    tf, quals.map{|qual| qual_names[qual] }.join(','), 
    *[peaks_fn, peaks_confirmed_fn].map{|fn|
      File.exist?(fn) ? File.readlines(fn).map(&:strip).reject(&:empty?).size : '-'
    }
  ]
  puts infos.join("\t")
end
