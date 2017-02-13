def convert_fasta_to_plain(genome_folder)
  Dir.glob(File.join(genome_folder, '*.fa')).each do |fasta_filename|
    num_entries = 0
    File.open(fasta_filename.ext('.plain'), 'w') do |fw|
      File.open(fasta_filename) do |f|
        f.each_line do |line|
          if line.start_with?('>')
            num_entries += 1
            next
          end
          fw.print(line.strip)
        end
      end
    end
    raise "Error in #{fasta_filename}! More than one entry found!!!"  if num_entries > 1
  end
end

# Formally, we should download genome from SYNAPSE, but it's less trivial
#   while genome assembly seems to be the same
desc 'Download hg19 genome assembly (from ENSEMBL).'
task :download_genome_assembly do
  folder = '/home/ilya/iogen/genome/mm10'
  mkdir_p folder
  sh 'wget', "--directory-prefix=#{folder}", "http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz"
  sh 'tar', '-C', folder, '-zxf', File.join(folder, 'chromFa.tar.gz')
  convert_fasta_to_plain(folder)
  rm_f Dir.glob(File.join(folder, '*.fa'))
end

desc 'Combine all chromosomes into a single fasta file'
task :generate_whole_genome_fasta do
  File.open('all_chromosomes.fa', 'w') do |fw|
    # Mus musculus has only 40 chromosomes!
    ((1..19).map(&:to_s) + ['X', 'Y']).each do |chr|
      fw.puts ">chr#{chr}"
      fw.puts File.read("genome/chr#{chr}.plain")
    end
  end
end
