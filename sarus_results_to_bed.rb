require 'bioinform'

def is_number?(str)
  !!Float(str) rescue false
end

def load_dipwm_matrix(filename)
  rows = File.readlines(filename).map(&:strip).reject(&:empty?).map{|line| line.split }
  if rows[0].size == 16 && rows[0].all?{|el| is_number?(el) } # No header line
    raise unless rows.all?{|row| row.all?{|el| is_number?(el) } }
    rows.map{|row|
      row.map{|el| Float(el) }
    }
  else # drop header
    raise unless rows.drop(1).all?{|row| row.all?{|el| is_number?(el) } }
    rows.drop(1).map{|row|
      row.map{|el| Float(el) }
    }
  end
end


raise 'Specify motif file'  unless motif_fn = ARGV.shift
raise 'Specify `mono` or `di` mode'  unless mode = ARGV.shift

case mode.downcase
when 'mono'
  site_length = Bioinform::MotifModel::PWM.from_file(motif_fn).length
when 'di'
  site_length = load_dipwm_matrix(motif_fn).length + 1
else
  raise "Unknown mode `#{mode}` specified. Specify `mono` or `di`."
end

chr = nil
$stdin.each_line{|line|
  if line.start_with?('>')
    chr = line[1..-1].strip
  else
    score, pos, strand = line.chomp.split
    pos = pos.to_i
    puts [chr, pos, pos + site_length, score].join("\t")
  end
}
