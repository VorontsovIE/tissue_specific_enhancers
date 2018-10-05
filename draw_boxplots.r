# join <( cat top_3enh_list.txt | sort -k1,1 ) <( find motif_collection/ -xtype f | xargs -rn1 basename -s .pwm | ruby -e 'readlines.each{|l| m=l.chomp; puts [m.split(".").first, m].join("\t")}' | sort -k1,1 ) -t $'\t' > top_3enh_list_w_motif.txt
#
# find boxplot_lists -xtype f | xargs -n1 -I{} echo "join <( cat {} | cut -f2,3 | tail -n+2 | sed -re 's/\"//g' | sort -k2,2 ) <( find motif_collection/ -xtype f | xargs -rn1 basename -s .pwm | ruby -e 'readlines.each{|l| m=l.chomp; puts [m.split(\".\").first, m].join(\"\\t\")}' | sort -k1,1 ) -t \$'\\t'  -1 2  -2 1 | sponge {}"
#
#

#cl_motif_triples_fn <- 'top_3enh_list_w_motif.txt'
#enh_type <- 'Constituent_Enhs/'
#output_folder <- 'boxplot_Sid/'

##cl_motif_triples_fn <- 'constenhancers_interesting.txt'
##enh_type <- 'Constituent_Enhs/'
#cl_motif_triples_fn <- 'superenhancers_interesting.txt'
#enh_type <- 'SuperEnhs/'
#output_folder <- 'boxplot_interesting/'

boxplot_configs <- list(
  list('SpE27acCnE_TF_se_top3.tsv', 'Constituent_Enhs/', 'SuperEnhancers_27ac/', 'boxplots_top_factors/SuperEnhancers_27ac/Constituent_Enhs/se_top3/'),
  list('SpE27acCnE_TF_te_top3.tsv', 'Constituent_Enhs/', 'SuperEnhancers_27ac/', 'boxplots_top_factors/SuperEnhancers_27ac/Constituent_Enhs/te_top3/'),

  list('SpE27acCnE_TF_se_top5.tsv', 'Constituent_Enhs/', 'SuperEnhancers_27ac/', 'boxplots_top_factors/SuperEnhancers_27ac/Constituent_Enhs/se_top5/'),
  list('SpE27acCnE_TF_te_top5.tsv', 'Constituent_Enhs/', 'SuperEnhancers_27ac/', 'boxplots_top_factors/SuperEnhancers_27ac/Constituent_Enhs/te_top5/'),
  list('SuECnE_TF_se_top5.tsv', 'Constituent_Enhs/', 'SuperEnhancers/', 'boxplots_top_factors/SuperEnhancers/Constituent_Enhs/se_top5/'),
  list('SuECnE_TF_te_top5.tsv', 'Constituent_Enhs/', 'SuperEnhancers/', 'boxplots_top_factors/SuperEnhancers/Constituent_Enhs/te_top5/'),

  list('SpE27acSpE_TF_se_top5.tsv', 'SuperEnhs/', 'SuperEnhancers_27ac/', 'boxplots_top_factors/SuperEnhancers_27ac/SuperEnhs/se_top5/'),
  list('SpE27acSpE_TF_te_top5.tsv', 'SuperEnhs/', 'SuperEnhancers_27ac/', 'boxplots_top_factors/SuperEnhancers_27ac/SuperEnhs/te_top5/'),
  list('SuESpE_TF_se_top5.tsv', 'SuperEnhs/', 'SuperEnhancers/', 'boxplots_top_factors/SuperEnhancers/SuperEnhs/se_top5/'),
  list('SuESpE_TF_te_top5.tsv', 'SuperEnhs/', 'SuperEnhancers/', 'boxplots_top_factors/SuperEnhancers/SuperEnhs/te_top5/'),

  list('top_3enh_list_w_motif.txt', 'Constituent_Enhs/', 'SuperEnhancers_27ac/', 'boxplots_top_3enh_list/SuperEnhancers_27ac/Constituent_Enhs/'),
  list('constenhancers_interesting.txt', 'Constituent_Enhs/', 'SuperEnhancers_27ac/', 'boxplots_interesting/SuperEnhancers_27ac/Constituent_Enhs/'),
  list('superenhancers_interesting.txt', 'SuperEnhs/', 'SuperEnhancers_27ac/', 'boxplots_interesting/SuperEnhancers_27ac/SuperEnhs/')
)

for (boxplot_config in boxplot_configs) {
  cl_motif_triples_fn <- boxplot_config[[1]]
  enh_type <- boxplot_config[[2]]
  detection_type <- boxplot_config[[3]]
  output_folder <- boxplot_config[[4]]

  dir.create(paste(output_folder, 'numsites/', sep=''), showWarnings = FALSE, recursive = TRUE)
  dir.create(paste(output_folder, 'homoscore/', sep=''), showWarnings = FALSE, recursive = TRUE)

  factor_list <- read.csv(file=paste('boxplot_lists/', cl_motif_triples_fn, sep=''), sep="\t", header=FALSE, col.names=c('TF', 'CL', 'MOTIF'))

  for (i in 1:length(factor_list$TF)) {
    CL <- factor_list$CL[i]
    # TF <- factor_cl$TF
    MOTIF <- factor_list$MOTIF[i]

    sites <- read.csv(file=paste('joined_score_stats/', detection_type, enh_type, CL, '.', MOTIF, '.txt', sep=''), sep="\t")

    png(file = paste(output_folder, 'homoscore/', CL, '.', MOTIF, '.png', sep=''))
    boxplot(HomotypicScore ~ EnhancerType, sites, varwidth=TRUE, notch=FALSE)
    dev.off()

    png(file = paste(output_folder, 'numsites/', CL, '.', MOTIF, '.png', sep=''))
    boxplot(NumSites ~ EnhancerType, sites, varwidth=TRUE, notch=FALSE)
    dev.off()
  }
}
