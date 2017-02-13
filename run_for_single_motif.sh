MotifName=$1
MotifCollectionFolder=$2
MotifThresholdsFolder=$3
ScoringScript=$4
OutputFolder=$5
MotifMode=$6
PvalueCutoff=${7:-0.0005}

MotifFn=${MotifCollectionFolder}/${MotifName}
ThresholdFn=${MotifThresholdsFolder}/${MotifName%.*pwm}.thr

mkdir -p ${OutputFolder}
${ScoringScript} ${MotifFn} ${PvalueCutoff} | \
  ruby sarus_results_to_bed.rb ${MotifFn} ${MotifMode} | \
  ruby score_to_pvalue.rb ${ThresholdFn} --precision 2 \
  > ${OutputFolder}/${MotifName%.*pwm}.txt
