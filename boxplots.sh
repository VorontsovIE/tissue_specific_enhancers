TFS=`find gtrd/adaptive_quality/ -xtype f -iname '*_MOUSE.bed' -size +0 | xargs -rn1 basename -s.bed`

for CL in BAT CH12 Kidney MEF Spleen Cortex MEL BmarrowDm Limb Testis Esb4 OlfactoryBulb Bmarrow Liver Thymus Es-E14 Placenta Cerebellum Lung Wbrain Heart SmallIntestine; do

# for folder in 'SuperEnhancers/Constituent_Enhs' 'SuperEnhancers_27ac/Constituent_Enhs' ; do
for folder in 'SuperEnhancers_27ac/Constituent_Enhs' ; do
  mkdir -p score_stats/$folder
  mkdir -p raw_scores/$folder
  mkdir -p sites_in_enhancers/$folder
  mkdir -p joined_score_stats/$folder
  for TF in $TFS ; do
    MOTIFS=`find motif_collection/ -xtype f -iname "${TF}*" | xargs -n1 --no-run-if-empty basename -s .pwm`
    for MOTIF in $MOTIFS; do
      ruby score_distribution.rb TSRE/mm10/$folder/${CL}.te.const.txt.bed gtrd/adaptive_quality/${TF}.bed $MOTIF  0.001  TE \
           sites_in_enhancers/$folder/${CL}.te.const.${MOTIF}.txt \
           raw_scores/$folder/${CL}.te.const.${MOTIF}.txt  score_stats/$folder/${CL}.te.const.${MOTIF}.txt
      ruby score_distribution.rb TSRE/mm10/$folder/${CL}.se.const.txt.bed gtrd/adaptive_quality/${TF}.bed $MOTIF  0.001  SE \
           sites_in_enhancers/$folder/${CL}.se.const.${MOTIF}.txt  \
           raw_scores/$folder/${CL}.se.const.${MOTIF}.txt  score_stats/$folder/${CL}.se.const.${MOTIF}.txt
#      ( \
#        echo $'chr\tfrom\tto\tMaxScore\tNumSites\tLen\tLambda\tHomotypicScore\tEnhancerType'; \
#        tail -n+2  score_stats/$folder/${CL}.te.const.${MOTIF}.txt; \
#        tail -n+2  score_stats/$folder/${CL}.se.const.${MOTIF}.txt; \
#      ) > joined_score_stats/$folder/${CL}.${MOTIF}.txt
    done
  done
done


# for folder in 'SuperEnhancers/SuperEnhs' 'SuperEnhancers_27ac/SuperEnhs' ; do
#   mkdir -p score_stats/$folder
#   mkdir -p raw_scores/$folder
#   mkdir -p chipseq_enhancer_bed/$folder
#   mkdir -p joined_score_stats/$folder
#   for TF in $TFS ; do
#     MOTIFS=`find motif_collection/ -xtype f -iname "${TF}*" | xargs -n1 --no-run-if-empty basename -s .pwm`
#     for MOTIF in $MOTIFS; do
#       #ruby score_distribution.rb TSRE/mm10/$folder/${CL}.te.bed gtrd/adaptive_quality/${TF}.bed $MOTIF  0.001 chipseq_enhancer_bed/$folder/${CL}.te.${MOTIF}.txt  raw_scores/$folder/${CL}.te.${MOTIF}.txt  score_stats/$folder/${CL}.te.${MOTIF}.txt TE
#       #ruby score_distribution.rb TSRE/mm10/$folder/${CL}.se.bed gtrd/adaptive_quality/${TF}.bed $MOTIF  0.001 chipseq_enhancer_bed/$folder/${CL}.se.${MOTIF}.txt  raw_scores/$folder/${CL}.se.${MOTIF}.txt  score_stats/$folder/${CL}.se.${MOTIF}.txt SE
#       ( \
#         echo $'chr\tfrom\tto\tMaxScore\tNumSites\tLen\tLambda\tHomotypicScore\tEnhancerType'; \
#         tail -n+2  score_stats/$folder/${CL}.te.${MOTIF}.txt; \
#         tail -n+2  score_stats/$folder/${CL}.se.${MOTIF}.txt; \
#       ) > joined_score_stats/$folder/${CL}.${MOTIF}.txt
#     done
#   done
# done

done


# R --slave -e 'x <- read.csv(file="stdin", sep="\t"); c(median(x$MaxScore), median(x$HomotypicScore))'
# R --slave -e 'x <- read.csv(file="stdin", sep="\t"); summary(x)'
