require 'fileutils'
require 'shellwords'

def num_lines(fn)
  # cmd = ['wc', '-l', fn].shelljoin
  # `#{cmd}`.to_i
  File.foreach(fn).count
end

FileUtils.mkdir_p 'gtrd/adaptive_quality'
# take `quality_default_level` best datasets. If resulting file is empty, add dataset with the next quality level
quality_default_level = 2
quality_levels = ['highest', 'high', 'medium', 'low', 'lowest']
tf_list = Dir.glob('gtrd/highest/*.bed').map{|fn| File.basename(fn, '.bed')[0..-3] }
tf_list.each do |tf|
  tf_chipseq_filelist = (quality_default_level .. quality_levels.size).map{|quality_level|
    quality_levels.first(quality_level).flat_map{|level|
      Dir.glob("gtrd/#{level}/#{tf}.?.bed")
    }
  }.detect{|filelist|
    filelist.any?{|fn| num_lines(fn) > 0 }
  }

  pipeline = [
    ['cat', *tf_chipseq_filelist].shelljoin,
     'sort -k1,1 -k2,2n',
     'bedtools merge',
  ]
  cmd = pipeline.join(' | ') + " > gtrd/adaptive_quality/#{tf}.bed"
  # system(cmd)
  puts cmd
end
