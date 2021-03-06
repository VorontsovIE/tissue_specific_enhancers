wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M15/gencode.vM15.annotation.gff3.gz
zcat gencode.vM15.annotation.gff3.gz | ruby extract_promoters.rb | bedtools flank -g ~/genome_sizes/mm10.genome -s -l 500 -r 0 | bedtools slop -g ~/genome_sizes/mm10.genome -s -l 0 -r 100 | bedtools sort > mm10_gencode_promoters.bed
cat pc_genes_coordinates_mm9.txt | awk -e '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $1 "\t" $5}' > pc_genes_coordinates_mm9.bed
./liftOver pc_genes_coordinates_mm9.bed liftOverChains/mm9ToMm10.over.chain.gz pc_genes_coordinates_mm10.bed pc_genes_coordinates_mm10.bed.unlifted
cat pc_genes_coordinates_mm10.bed | bedtools flank -g ~/genome_sizes/mm10.genome -s -l 500 -r 0 | bedtools slop -g ~/genome_sizes/mm10.genome -s -l 0 -r 100 | bedtools sort > mm10_pc_genes_promoters.bed

wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod +x liftOver

sudo yum install moreutils

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
sed -ie 's/^FOXF1_MOUSE\tFoxf1$/FOXF1_MOUSE\tFoxf1a/' uniprot_main_gene_name.tsv
sed -ie 's/^ZN281_MOUSE\tZnf281$/ZN281_MOUSE\tZfp281/' uniprot_main_gene_name.tsv
sed -ie 's/^ZN431_MOUSE\tZnf431$/ZN431_MOUSE\tZfp932/' uniprot_main_gene_name.tsv
##########################

# mkdir -p gtrd/AB_quality
# find gtrd/ -xtype f -iname '*.bed' | grep -Pe 'gtrd/(high|higest)/' | xargs -n1 basename -s .bed | ruby -e 'puts readlines.map{|l| l.chomp[0..-3] }' | sort | uniq | xargs -I{} -n1 echo 'cat gtrd/highest/{}.A.bed gtrd/high/{}.B.bed | sort -k1,1 -k2,2n | bedtools merge > gtrd/AB_quality/{}.bed' | bash

ruby adaptive_chipseq_combiner.rb | bash

# Анализ покрытия генома чипсеками
find gtrd/adaptive_quality/ -xtype f -iname '*.bed' | xargs -n1 basename -s '.bed' | xargs -n1 -I{} echo "echo -n {} ' '; cat gtrd/adaptive_quality/{}.bed | ruby -e 'puts readlines.map{|l| s,f = l.split[1,2].map(&:to_f); f-s}.inject(0.0, &:+) / 1e6'" | bash | ruby -e 'puts readlines.sort_by{|l| l.split[-1].to_f }'

##############################
# Find site occurences
find gtrd/adaptive_quality/ -iname '*_MOUSE.bed' | xargs -n1 basename -s .bed | sort | xargs -n1 -I{} echo 'find motif_collection/ -iname "{}*"' | bash | xargs -n1 basename | xargs -n1 -I MOTIF echo './run_for_single_motif.sh MOTIF motif_collection motif_thresholds ./get_scores.sh sites mono 0.0005' | parallel -j8

# Find GTRD peaks confirmed by motif occurences
mkdir -p gtrd/confirmed_by_motif
find gtrd/adaptive_quality/ -iname '*_MOUSE.bed' | xargs -n1 basename -s .bed | sort | xargs -n1 -I{} echo 'cat `find sites/ -iname "{}*"` | sort -k1,1 -k2,2n | bedtools merge | bedtools intersect -a gtrd/adaptive_quality/{}.bed -b - -u > gtrd/confirmed_by_motif/{}.bed' | parallel -j8

##############################
for FOLDER in  Active_promoters  Strong_enhancers  SuperEnhancers/Constituent_Enhs  SuperEnhancers/SuperEnhs  SuperEnhancers_27ac/Constituent_Enhs  SuperEnhancers_27ac/SuperEnhs ; do
  echo "${FOLDER}"
  # find $FOLDER -xtype f -iname '*.bed' | xargs -n1 -I{} basename -s .bed "{}" | xargs -n1 -I{} awk -e '{print $1 "\t" $2 "\t" $3 "\t" "{}_"$4}' ${FOLDER}/{}.bed | sort -k1,1 -k2,2n > ${FOLDER}.bed

  mkdir -p "TSRE/mm10/${FOLDER}"

  find "TSRE/${FOLDER}" -xtype f -iname '*.bed' | xargs -n1 -I{} basename -s .bed "{}" | xargs -n1 -I{} echo "./liftOver 'TSRE/${FOLDER}/{}.bed' liftOverChains/mm9ToMm10.over.chain.gz 'TSRE/mm10/${FOLDER}/{}.bed' 'TSRE/mm10/${FOLDER}/{}.bed.unlifted'" | bash

  mkdir -p "TSRE/mm10_merged/${FOLDER}"
  find "TSRE/mm10/${FOLDER}" -xtype f -iname '*.bed' | xargs -n1 -I{} basename -s .bed "{}" | xargs -n1 -I{} echo "cat 'TSRE/mm10/${FOLDER}/{}.bed' | sort -k1,1 -k2,2n | bedtools merge > 'TSRE/mm10_merged/${FOLDER}/{}.bed'" | bash

  # Analysis of genome coverage (in megabases)
  find "TSRE/mm10_merged/${FOLDER}/" -xtype f -iname '*.bed' | xargs -n1 -I{} basename -s .bed "{}" | xargs -n1 -I{} echo "echo -n {} ' '; cat 'TSRE/mm10_merged/${FOLDER}/{}.bed' | ruby -e 'puts readlines.map{|l| s,f = l.split[1,2].map(&:to_f); f-s}.inject(0.0, &:+) / 1e6'" | bash | ruby -e 'puts readlines.sort_by{|l| l.split[-1].to_f }'
  # Interesting fact:
  #   Promoters in Placenta cover only 0.0174 Mbp while promoters in other tissues cover approximately 0.5 - 5.0 Mbp. So possibly, it's just a bad sample
  #   Enhancers in Placenta are ok

  # for CL in `find "mm10_merged/${FOLDER}/" -xtype f -iname '*.bed' | xargs -n1 -I{} basename -s .bed "{}"` ; do
  #   for MOTIF in `find gtrd/adaptive_quality/ -xtype f -iname '*_MOUSE.bed' | xargs -n1 -I{} basename -s .bed "{}"` ; do
  #     # COVERED_ELEMENTS=`bedtools intersect -c -a "gtrd/adaptive_quality/${MOTIF}.bed" -b "mm10_merged/${FOLDER}/${CL}.bed" | cut -f4 | sort | uniq -c | ruby -e 'counts=readlines.map{|l| l.split.map(&:to_i) }; puts counts.map{|a,b| a*b}.sum'`
  #     bedtools intersect -wao -a "gtrd/adaptive_quality/${MOTIF}.bed" -b "mm10_merged/${FOLDER}/${CL}.bed" | bedtools groupby -c 7 -o sum
  #   done
  # done
done

# ruby calculate_tau_scores.rb | sort -k1,1 > tau_scores.tsv

for FOLDER in  Active_promoters  Strong_enhancers  SuperEnhancers/Constituent_Enhs  SuperEnhancers/SuperEnhs  SuperEnhancers_27ac/Constituent_Enhs  SuperEnhancers_27ac/SuperEnhs ; do
  ##########################################################
  # Calculate significances of "enhancer-specific" binding #
  ##########################################################
  mkdir -p "TSRE/results_all_bound/intermediate/${FOLDER}"
  find "TSRE/mm10_merged/${FOLDER}" -iname '*.bed' | xargs -n1 -I{} basename -s .bed "{}" | sort | xargs -n1 -I{} echo \
    "ruby overlapped_enhancers.rb 'TSRE/mm10_merged/${FOLDER}/{}.bed' 100 gtrd/adaptive_quality > 'TSRE/results_all_bound/intermediate/${FOLDER}/{}.tsv'" \
    | parallel -j8

  mkdir -p "TSRE/results_bound_motif_confirmed/intermediate/${FOLDER}"
  find "TSRE/mm10_merged/${FOLDER}" -iname '*.bed' | xargs -n1 -I{} basename -s .bed "{}" | sort | xargs -n1 -I{} echo \
    "ruby overlapped_enhancers.rb 'TSRE/mm10_merged/${FOLDER}/{}.bed' 100 gtrd/confirmed_by_motif > 'TSRE/results_bound_motif_confirmed/intermediate/${FOLDER}/{}.tsv'" \
    | parallel -j8
  ##########################################################
done

for FOLDER in  Active_promoters  Strong_enhancers  SuperEnhancers/Constituent_Enhs  SuperEnhancers/SuperEnhs  SuperEnhancers_27ac/Constituent_Enhs  SuperEnhancers_27ac/SuperEnhs ; do
  ###################################################################
  # Join gene name and tissue-specific expression to significances  #
  ###################################################################
  mkdir -p "TSRE/results_all_bound/${FOLDER}"
  find "TSRE/results_all_bound/intermediate/${FOLDER}/" -iname '*.tsv' | xargs -n1 -I{} basename "{}" | xargs -n1 -I{} echo \
    "cat 'TSRE/results_all_bound/intermediate/${FOLDER}/{}' | ruby join_infos.rb '{}' | sponge 'TSRE/results_all_bound/${FOLDER}/{}'" \
    | parallel -j8

  mkdir -p "TSRE/results_bound_motif_confirmed/${FOLDER}"
  find "TSRE/results_bound_motif_confirmed/intermediate/${FOLDER}/" -iname '*.tsv' | xargs -n1 -I{} basename "{}" | xargs -n1 -I{} echo \
    "cat 'TSRE/results_bound_motif_confirmed/intermediate/${FOLDER}/{}' | ruby join_infos.rb '{}' | sponge 'TSRE/results_bound_motif_confirmed/${FOLDER}/{}'" \
    | parallel -j8
  ###################################################################
done

for FOLDER in  Active_promoters  Strong_enhancers  SuperEnhancers/Constituent_Enhs  SuperEnhancers/SuperEnhs  SuperEnhancers_27ac/Constituent_Enhs  SuperEnhancers_27ac/SuperEnhs ; do
  ##########################################################
  # Join motif confirmed results to motif independent ones #
  ##########################################################

  mkdir -p "TSRE/results_combined/${FOLDER}"
  find "TSRE/results_all_bound/${FOLDER}/" -xtype f | xargs -n1 basename | xargs -n1 -I{} echo \
    "ruby glue_results.rb 'TSRE/results_all_bound/${FOLDER}/{}' 'TSRE/results_bound_motif_confirmed/${FOLDER}/{}' | cut -f1-7,15- > 'TSRE/results_combined/${FOLDER}/{}'" \
    | bash
  #####################################################################################
done


# find tissue_specificity/ -xtype f -iname '*.tsv' | xargs -n1 -I{} echo 'ruby correct_tissue_specificity_ids.rb "{}" > "{}"_' | bash
# find tissue_specificity/ -xtype f -iname '*.tsv_' | xargs -n1 -I{} basename -s .tsv_ "{}" | xargs -n1 -I{} echo 'mv "tissue_specificity/{}.tsv_" "tissue_specificity/{}.tsv"' | bash

# # leave only in-gtrd samples
# for CL_FN in `find tissue_specificity -xtype f -iname '*.tsv'` ; do
#   find gtrd/adaptive_quality/ -xtype f -size +0 | xargs -n1 basename -s '.bed' | sort | awk -e '{print "\\b" $0 "\\b"}' | grep -f - "${CL_FN}" | cut -f5 | ruby -e 'puts ARGF.readlines.map(&:chomp).reject(&:empty?).inject([]){|res, mot| (res.include?(mot) || res.size > 20) ? res : res + [mot] }' > "${CL_FN}_selected"
# done

# for RES_FOLDER in results_combined/TS_strong_enhancers results_combined/TS_active_promoters ; do
#   ruby measure_tissue_specificity.rb  ${RES_FOLDER}/Cerebellum.tsv     tissue_specificity/Cerebellum.FF15-8B2.tsv_selected
#   ruby measure_tissue_specificity.rb  ${RES_FOLDER}/Cortex.tsv         tissue_specificity/Cortex.FF12-14D5.tsv_selected
#   ruby measure_tissue_specificity.rb  ${RES_FOLDER}/Heart.tsv          tissue_specificity/Heart.FF1390-42D2.tsv_selected
#   ruby measure_tissue_specificity.rb  ${RES_FOLDER}/Kidney.tsv         tissue_specificity/Kidney.FF1385-42H1.tsv_selected
#   ruby measure_tissue_specificity.rb  ${RES_FOLDER}/Liver.tsv          tissue_specificity/Liver.FF1382-42D1.tsv_selected
#   ruby measure_tissue_specificity.rb  ${RES_FOLDER}/Lung.tsv           tissue_specificity/Lung.FF28-22B1.tsv_selected
#   ruby measure_tissue_specificity.rb  ${RES_FOLDER}/Placenta.tsv       tissue_specificity/Placenta.FF577-18G3.tsv_selected
#   ruby measure_tissue_specificity.rb  ${RES_FOLDER}/SmallIntestine.tsv tissue_specificity/Small_intestine.FF790-21I1.tsv_selected
#   ruby measure_tissue_specificity.rb  ${RES_FOLDER}/Spleen.tsv         tissue_specificity/Spleen.FF25-2G2.tsv_selected
#   ruby measure_tissue_specificity.rb  ${RES_FOLDER}/Testis.tsv         tissue_specificity/Testis.FF57-7G5.tsv_selected
#   ruby measure_tissue_specificity.rb  ${RES_FOLDER}/Thymus.tsv         tissue_specificity/Thymus.FF38-12B5.tsv_selected
# done
