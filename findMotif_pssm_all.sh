#! /bin/bash

if [[ $# -ne 4 ]]
then
    echo 'Usage: ./'$0 ' infile1 infile2 infile3 outfile'  
    exit 1
fi 

infile1=$1
infile2=$2
infile3=$3
outfile=$4

echo "plotting Motif for current task"
Rscript $(dirname "$0")/findMotif_pssm_all.r -f $infile1 -g $infile2 -i $infile3 -o $outfile
wait

file_list=$(ls -1v *_${outfile}*_all.pdf|tr '\n' ' ')
mutool merge -o ${outfile}_all.pdf $file_list #&> /dev/null
mupdf-x11 ${outfile}_all.pdf &> /dev/null

metagene_list=$(ls -1v *_${outfile}*_mRNA_metagene.pdf|tr '\n' ' ')
mutool merge -o ${outfile}_mRNA_metagene.pdf $metagene_list #&> /dev/null
mupdf-x11 ${outfile}_mRNA_metagene.pdf &> /dev/null

echo -e "Finished: findMotif_pssm done!\n" 