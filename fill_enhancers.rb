require 'fileutils'

def file_with_filled_gaps_cmd(filename, pattern_size: 200, genome_sizes_fn: 'genome_sizes/mm10.genome')
  [
    "cat #{filename}",
    "bedtools sort",
    "bedtools slop -g #{genome_sizes_fn} -b #{pattern_size}", # dilute
    "bedtools merge",
    "bedtools slop -g #{genome_sizes_fn} -b #{-pattern_size}", # erode
  ].join('|')
end

enhancers_folder = 'TSRE/SuperEnhancers_27ac/Constituent_Enhs/'
glued_enhancers_folder = 'TSRE/glued/SuperEnhancers_27ac/Constituent_Enhs/'
FileUtils.mkdir_p  glued_enhancers_folder

Dir.glob(File.join(enhancers_folder, '*.bed')).each{|enhancers_fn|
  glued_fn = File.join(glued_enhancers_folder, File.basename(enhancers_fn, '.bed') + '.glued.bed')
  cmd = file_with_filled_gaps_cmd(enhancers_fn) + " > #{glued_fn}"
  puts cmd
}

