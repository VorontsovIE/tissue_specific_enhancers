rake download_genome_assembly
ln -s /home/ilya/iogen/genome/mm10/ genome
rake generate_whole_genome_fasta

# mkdir -p gtrd/AB_quality
# find gtrd/ -xtype f -iname '*.bed' | grep -Pe 'gtrd/(high|higest)/' | xargs -n1 basename -s .bed | ruby -e 'puts readlines.map{|l| l.chomp[0..-3] }' | sort | uniq | xargs -I{} -n1 echo 'cat gtrd/highest/{}.A.bed gtrd/high/{}.B.bed | sort -k1,1 -k2,2n | bedtools merge > gtrd/AB_quality/{}.bed' | bash

ruby adaptive_chipseq_combiner.rb | bash

# Анализ покрытия генома чипсеками
find gtrd/adaptive_quality/ -xtype f -iname '*.bed' | xargs -n1 basename -s '.bed' | xargs -n1 -I{} echo "echo -n {} ' '; cat gtrd/adaptive_quality/{}.bed | ruby -e 'puts readlines.map{|l| s,f = l.split[1,2].map(&:to_f); f-s}.inject(0.0, &:+) / 1e6'" | bash | ruby -e 'puts readlines.sort_by{|l| l.split[-1].to_f }'

##############################

FOLDER=TS_active_promoters
# FOLDER=TS_strong_enhancers

# find $FOLDER -xtype f -iname '*.bed' | xargs -n1 basename -s .bed | xargs -n1 -I{} awk -e '{print $1 "\t" $2 "\t" $3 "\t" "{}_"$4}' ${FOLDER}/{}.bed | sort -k1,1 -k2,2n > ${FOLDER}.bed

mkdir -p mm10_${FOLDER}

find $FOLDER -xtype f -iname '*.bed' | xargs -n1 basename -s .bed | xargs -n1 -I{} echo "liftOver ${FOLDER}/{}.bed liftOverChains/mm9ToMm10.over.chain.gz mm10_${FOLDER}/{}.bed mm10_${FOLDER}/{}.bed.unlifted" | bash

mkdir -p "merged_mm10_${FOLDER}"
find $FOLDER -xtype f -iname '*.bed' | xargs -n1 basename -s .bed | xargs -n1 -I{} echo "cat mm10_${FOLDER}/{}.bed | sort -k1,1 -k2,2n | bedtools merge > merged_mm10_${FOLDER}/{}.bed" | bash

# Анализ покрытия генома (в мегабазах)
find merged_mm10_${FOLDER}/ -xtype f -iname '*.bed' | xargs -n1 basename -s '.bed' | xargs -n1 -I{} echo "echo -n {} ' '; cat merged_mm10_${FOLDER}/{}.bed | ruby -e 'puts readlines.map{|l| s,f = l.split[1,2].map(&:to_f); f-s}.inject(0.0, &:+) / 1e6'" | bash | ruby -e 'puts readlines.sort_by{|l| l.split[-1].to_f }'
# Интересный факт:
#   Промотеры в Placenta - всего 0.0174 мегабаз, тогда как промотеры в остальных тканях занимают порядка 0.5 - 5.0 мегабаз. 
#   С энхансерами всё ок.

# for CL in `find merged_mm10_${FOLDER}/ -xtype f -iname '*.bed' | xargs -n1 basename -s '.bed'` ; do
#   for MOTIF in `find gtrd/adaptive_quality/ -xtype f -iname '*_MOUSE.bed' | xargs -n1 basename -s '.bed'` ; do
#     # COVERED_ELEMENTS=`bedtools intersect -c -a gtrd/adaptive_quality/${MOTIF}.bed -b merged_mm10_${FOLDER}/${CL}.bed | cut -f4 | sort | uniq -c | ruby -e 'counts=readlines.map{|l| l.split.map(&:to_i) }; puts counts.map{|a,b| a*b}.sum'`
#     bedtools intersect -wao -a gtrd/adaptive_quality/${MOTIF}.bed -b merged_mm10_${FOLDER}/${CL}.bed | bedtools groupby -c 7 -o sum
#   done
# done

mkdir -p results/${FOLDER}
find merged_mm10_${FOLDER} -iname '*.bed' | xargs -n1 basename -s .bed | sort | xargs -n1 -I{} echo "ruby overlapped_enhancers.rb merged_mm10_${FOLDER}/{}.bed 100 > results/${FOLDER}/{}.txt" | bash



find gtrd/adaptive_quality/ -iname '*_MOUSE.bed' | xargs -n1 basename -s .bed | sort | xargs -n1 -I{} echo 'find motif_collection/ -iname "{}*"' | bash | xargs -n1 basename | xargs -n1 -I MOTIF echo './run_for_single_motif.sh MOTIF motif_collection motif_thresholds ./get_scores.sh sites mono 0.0005' | parallel -j8

mkdir -p gtrd/confirmed_by_motif
find gtrd/adaptive_quality/ -iname '*_MOUSE.bed' | xargs -n1 basename -s .bed | sort | xargs -n1 -I{} echo 'cat `find sites/ -iname "{}*"` | sort -k1,1 -k2,2n | bedtools merge | bedtools intersect -a gtrd/adaptive_quality/{}.bed -b - -u > gtrd/confirmed_by_motif/{}.bed' | parallel -j8
