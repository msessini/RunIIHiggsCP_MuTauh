 qsub -q    -o  -e  -F "$1 $2" Set_.sh
 qsub -q  \${que}  -o \${output} -e \${error} -F "$1 $2" Set_.sh
 qsub -q  ${que}  -o ${output} -e ${error} -F "$1 $2" Set_.sh
