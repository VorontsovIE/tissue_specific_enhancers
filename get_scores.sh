MotifFn=$1
PvalueCutoff=${2:-0.0005}

Threshold=`java -cp ape.jar ru.autosome.ape.FindThreshold ${MotifFn}  ${PvalueCutoff} --background uniform | grep -v -P '^#' | cut -f 3`

java -Xmx2G -Xms2G -cp sarus.jar ru.autosome.SARUS all_chromosomes.fa ${MotifFn} ${Threshold} skipn
