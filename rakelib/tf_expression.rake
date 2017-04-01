require 'cgi'
require 'set'
require 'rake/clean'

def sample_names_from_peaks_file(tpm_filename)
  sample_names = []
  sample_name_pattern = /##ColumnVariables\[tpm\.(?<sample_name>.+)\]=/

  File.open(tpm_filename) do |f|
    f.each_line do |line|
      if match = line.match( sample_name_pattern )
        sample_name = CGI.unescape(match[:sample_name])
        sample_names << sample_name
      end
    end
  end
  sample_names
end


def gene_expressions_from_peaks_file(tpm_filename, target_hgnc_ids)
  gene_expressions = {} # Hash.new{|hsh, hgnc_id| hsh[hgnc_id] = [] }
  hgnc_pattern = /^HGNC:(?<hgnc>\d+)$/

  File.open(tpm_filename) do |f|
    f.each_line do |line|
      next  unless line.start_with? 'chr'
      $stderr.puts "Line #{f.lineno} "  if f.lineno % 10000 == 0

      row = line.strip.split("\t")
      anno, short_dsc, dsc, assoc, entrez, hgnc, uniprot = row.first(7)

      hgnc_ids = hgnc.split(',').map{|hgnc_term|
        hgnc_term.match(hgnc_pattern)
      }.map{|match|
        match && match[:hgnc]
      }.compact.map(&:to_i)

      next unless Set.new(hgnc_ids).intersect?(target_hgnc_ids)

      peak_expressions = row.drop(7).map(&:to_f)

      hgnc_ids.each do |hgnc_id|
        next unless target_hgnc_ids.include?(hgnc_id)

        gene_expressions[hgnc_id] ||= Array.new(peak_expressions.size, 0)
        peak_expressions.each_with_index do |peak_expression, idx|
          gene_expressions[hgnc_id][idx] += peak_expression
        end
      end
    end
  end
  gene_expressions
end

def extract_gene_expression_into_file(from_filename, to_filename, hgnc_symbols, target_hgnc_ids)
  sample_names = sample_names_from_peaks_file(from_filename)
  gene_expressions = gene_expressions_from_peaks_file(from_filename, target_hgnc_ids)

  File.open(to_filename, 'w') do |fw|
    fw.puts ['HGNC ID', 'HGNC symbol', *sample_names].join("\t")
    gene_expressions.each do |gene, gene_expressions_by_sample|
      fw.puts [gene, hgnc_symbols[gene], *gene_expressions_by_sample].join("\t")
    end
  end
end

namespace :data do
  desc 'Load gene expressions by samples'
  task :load_expressions => ['load_expressions:with_timecourses', 'load_expressions:without_timecourses']

  namespace :load_expressions do
    desc 'Load gene expressions by samples (with timecourses)'
    task :with_timecourses => ['public/public_data/gene_expressions_by_sample_with_timecourses.txt']

    desc 'Load gene expressions by samples (without timecourses)'
    task :without_timecourses => ['public/public_data/gene_expressions_by_sample.txt']

    task :download => ['download:tpm_without_timecourses', 'download:tpm_with_timecourses']
    namespace :download do
      directory 'source_data'
      desc 'Download and unpack TPMs without timecourses'
      task :tpm_without_timecourses => ['source_data/hg19.cage_peak_tpm_ann.osc.txt']

      desc 'Download and unpack TPMs with timecourses'
      task :tpm_with_timecourses => ['source_data/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz']

      file 'source_data/hg19.cage_peak_tpm_ann.osc.txt' => ['source_data/hg19.cage_peak_tpm_ann.osc.txt.gz'] do |t|
        sh 'gzip', '--keep', '-d', t.prerequisites.first
      end

      file 'source_data/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt' => ['source_data/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz'] do |t|
        sh 'gzip', '--keep', '-d', t.prerequisites.first
      end

      file 'source_data/hg19.cage_peak_tpm_ann.osc.txt.gz' => ['source_data'] do |t|
        sh 'wget', '-O', t.name, 'http://fantom.gsc.riken.jp/5/datafiles/phase1.3/extra/CAGE_peaks/hg19.cage_peak_tpm_ann.osc.txt.gz'
      end
      file 'source_data/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz' => ['source_data'] do |t|
        sh 'wget', '-O', t.name, 'http://fantom.gsc.riken.jp/5/datafiles/phase2.0/extra/CAGE_peaks/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz'
      end

      CLEAN.include('source_data/hg19.cage_peak_tpm_ann.osc.txt')
      CLEAN.include('source_data/hg19.cage_peak_tpm_ann.osc.txt.gz')

      CLEAN.include('source_data/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt')
      CLEAN.include('source_data/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz')
    end

    { 'public/public_data/gene_expressions_by_sample_with_timecourses.txt' => 'source_data/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt',
      'public/public_data/gene_expressions_by_sample.txt' => 'source_data/hg19.cage_peak_tpm_ann.osc.txt'
    }.each do |sample_expressions_file, tpm_filename|
      file sample_expressions_file => [tpm_filename, :load_epigenes, :load_histones] do |t|
        hgnc_symbols = {}
        $epigenes.select{|info| info[:hgnc_id] && info[:hgnc_id] != '-' }.each{|info| hgnc_symbols[info[:hgnc_id]] = info[:hgnc_symbol] }
        $histones.select{|info| info[:hgnc_id] && info[:hgnc_id] != '-' }.each{|info| hgnc_symbols[info[:hgnc_id]] = info[:hgnc_symbol] }
        target_hgnc_ids = Set.new(hgnc_symbols.keys)

        extract_gene_expression_into_file(t.prerequisites.first, t.name, hgnc_symbols, target_hgnc_ids)
      end
      CLOBBER.include(sample_expressions_file)
    end
  end
end
