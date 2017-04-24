rake download_genome_assembly
ln -s /home/ilya/iogen/genome/mm10/ genome
rake generate_whole_genome_fasta

##########################
# Make uniprot --> gene name mapping unique and consistent with expressions matrix
cat uniprot_gene_names.tsv | tail -n+2 | cut -f 2,5 | ruby -e 'puts readlines.map{|l| id,nm = l.chomp.split("\t",2); [id, nm.split(" 
").first].join("\t") }' | sort -k1,1 > uniprot_main_gene_name.tsv
sed -ie 's/^P53_MOUSE\tTp53$/P53_MOUSE\tTrp53/' uniprot_main_gene_name.tsv
sed -ie 's/^SPI1_MOUSE\tSpi1$/SPI1_MOUSE\tSfpi1/' uniprot_main_gene_name.tsv
sed -ie 's/^ZN143_MOUSE\tZnf143$/ZN143_MOUSE\tZfp143/' uniprot_main_gene_name.tsv
sed -ie 's/^ZN335_MOUSE\tZnf335$/ZN335_MOUSE\tZfp335/' uniprot_main_gene_name.tsv
sed -ie 's/^ZN322_MOUSE\tZnf322$/ZN322_MOUSE\tZfp322a/' uniprot_main_gene_name.tsv
sed -ie 's/^KMT2A_MOUSE\tKmt2a$/KMT2A_MOUSE\tMll1/' uniprot_main_gene_name.tsv
sed -ie 's/^KMT2B_MOUSE\tKmt2b$/KMT2B_MOUSE\tMll2/' uniprot_main_gene_name.tsv
sed -ie 's/^KAT8_MOUSE\tKat8$/KAT8_MOUSE\tMyst1/' uniprot_main_gene_name.tsv
##########################

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

find gtrd/adaptive_quality/ -iname '*_MOUSE.bed' | xargs -n1 basename -s .bed | sort | xargs -n1 -I{} echo 'find motif_collection/ -iname "{}*"' | bash | xargs -n1 basename | xargs -n1 -I MOTIF echo './run_for_single_motif.sh MOTIF motif_collection motif_thresholds ./get_scores.sh sites mono 0.0005' | parallel -j8

mkdir -p gtrd/confirmed_by_motif
find gtrd/adaptive_quality/ -iname '*_MOUSE.bed' | xargs -n1 basename -s .bed | sort | xargs -n1 -I{} echo 'cat `find sites/ -iname "{}*"` | sort -k1,1 -k2,2n | bedtools merge | bedtools intersect -a gtrd/adaptive_quality/{}.bed -b - -u > gtrd/confirmed_by_motif/{}.bed' | parallel -j8

# ruby calculate_tau_scores.rb | sort -k1,1 > tau_scores.tsv

##########################################################
# Calculate significances of "enhancer-specific" binding #
##########################################################
mkdir -p results_all_bound/${FOLDER}
find merged_mm10_${FOLDER} -iname '*.bed' | xargs -n1 basename -s .bed | sort | xargs -n1 -I{} echo \
  "ruby overlapped_enhancers.rb merged_mm10_${FOLDER}/{}.bed 100 gtrd/adaptive_quality > results_all_bound/${FOLDER}/{}.tsv" \
  | parallel -j24

mkdir -p results_bound_motif_confirmed/${FOLDER}
find merged_mm10_${FOLDER} -iname '*.bed' | xargs -n1 basename -s .bed | sort | xargs -n1 -I{} echo \
  "ruby overlapped_enhancers.rb merged_mm10_${FOLDER}/{}.bed 100 gtrd/confirmed_by_motif > results_bound_motif_confirmed/${FOLDER}/{}.tsv" \
  | parallel -j24
##########################################################


###################################################################
# Join gene name and tissue-specific expression to significances  #
###################################################################
find results_all_bound/${FOLDER}/ -iname '*.tsv' | xargs -n1 basename -s .tsv | xargs -n1 -I{} echo \
  "cat results_all_bound/${FOLDER}/{}.tsv | ruby join_infos.rb {} | sponge results_all_bound/${FOLDER}/{}.tsv" \
  | parallel -j24

find results_bound_motif_confirmed/${FOLDER}/ -iname '*.tsv' | xargs -n1 basename -s .tsv | xargs -n1 -I{} echo \
  "cat results_bound_motif_confirmed/${FOLDER}/{}.tsv | ruby join_infos.rb {} | sponge results_bound_motif_confirmed/${FOLDER}/{}.tsv" \
  | parallel -j24
###################################################################


#####################################################################################
# переименуем папки результатов, сделанных по чипсекам, поправленным/непоправленным #
# на наличие мотива в `results_shifted100` и `results_shifted100_withMotif`         #
#####################################################################################

mkdir -p results_combined/${FOLDER}
find results_all_bound/${FOLDER}/ -xtype f | xargs -n1 basename | xargs -n1 -I{} echo \
  "ruby glue_results.rb results_all_bound/${FOLDER}/{} results_bound_motif_confirmed/${FOLDER}/{} | cut -f1-7,13- > results_combined/${FOLDER}/{}" \
  | bash

#####################################################################################

find tissue_specificity/ -xtype f -iname '*.tsv' | xargs -n1 -I{} echo 'ruby correct_tissue_specificity_ids.rb "{}" > "{}"_' | bash
find tissue_specificity/ -xtype f -iname '*.tsv_' | xargs -n1 -I{} basename -s .tsv_ "{}" | xargs -n1 -I{} echo 'mv "tissue_specificity/{}.tsv_" "tissue_specificity/{}.tsv"' | bash

# leave only in-gtrd samples
for CL_FN in `find tissue_specificity -xtype f -iname '*.tsv'` ; do
  find gtrd/adaptive_quality/ -xtype f -size +0 | xargs -n1 basename -s '.bed' | sort | awk -e '{print "\\b" $0 "\\b"}' | grep -f - "${CL_FN}" | cut -f5 | ruby -e 'puts ARGF.readlines.map(&:chomp).reject(&:empty?).inject([]){|res, mot| (res.include?(mot) || res.size > 20) ? res : res + [mot] }' > "${CL_FN}_selected"
done

for RES_FOLDER in results_combined/TS_strong_enhancers results_combined/TS_active_promoters ; do
  ruby measure_tissue_specificity.rb  ${RES_FOLDER}/Cerebellum.tsv     tissue_specificity/Cerebellum.FF15-8B2.tsv_selected
  ruby measure_tissue_specificity.rb  ${RES_FOLDER}/Cortex.tsv         tissue_specificity/Cortex.FF12-14D5.tsv_selected
  ruby measure_tissue_specificity.rb  ${RES_FOLDER}/Heart.tsv          tissue_specificity/Heart.FF1390-42D2.tsv_selected
  ruby measure_tissue_specificity.rb  ${RES_FOLDER}/Kidney.tsv         tissue_specificity/Kidney.FF1385-42H1.tsv_selected
  ruby measure_tissue_specificity.rb  ${RES_FOLDER}/Liver.tsv          tissue_specificity/Liver.FF1382-42D1.tsv_selected
  ruby measure_tissue_specificity.rb  ${RES_FOLDER}/Lung.tsv           tissue_specificity/Lung.FF28-22B1.tsv_selected
  ruby measure_tissue_specificity.rb  ${RES_FOLDER}/Placenta.tsv       tissue_specificity/Placenta.FF577-18G3.tsv_selected
  ruby measure_tissue_specificity.rb  ${RES_FOLDER}/SmallIntestine.tsv tissue_specificity/Small_intestine.FF790-21I1.tsv_selected
  ruby measure_tissue_specificity.rb  ${RES_FOLDER}/Spleen.tsv         tissue_specificity/Spleen.FF25-2G2.tsv_selected
  ruby measure_tissue_specificity.rb  ${RES_FOLDER}/Testis.tsv         tissue_specificity/Testis.FF57-7G5.tsv_selected
  ruby measure_tissue_specificity.rb  ${RES_FOLDER}/Thymus.tsv         tissue_specificity/Thymus.FF38-12B5.tsv_selected
done
