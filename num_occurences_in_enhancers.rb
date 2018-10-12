require 'fileutils'
FileUtils.mkdir_p 'num_occurences_per_enhancer'

Dir.glob("motif_occurences_in_cistrome/*.bed").each{|tf_fn|
  tf = File.basename(tf_fn, '.bed')
  FileUtils.mkdir_p "num_occurences_per_enhancer/#{tf}"

  Dir.glob("TSRE/glued/SuperEnhancers_27ac/Constituent_Enhs/*.bed").each{|glued_enhancers_fn|
    enhancer = File.basename(glued_enhancers_fn, '.glued.bed')
    output_fn = "num_occurences_per_enhancer/#{tf}/#{tf}.#{enhancer}.bed"
    cmd = [
      "bedtools intersect -sorted -loj -a #{glued_enhancers_fn} -b <( echo $'chr1\\t1\\t1'; cat #{tf_fn} )",
      "bedtools groupby -c 5,6 -o collapse,collapse",
      "ruby calculate_overlap_size.rb",
    ].join(' | ') + " > #{output_fn}"
    puts cmd  unless File.exist?(output_fn) && File.size(output_fn) > 0
  }
}
