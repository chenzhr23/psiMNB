#! /bin/bash

if [[ $# -ne 2 ]]
then
    echo 'Usage: ./'$0 ' bedfile total_polyA_mix_data'  
    exit 1
fi 

bedfile=$1
total_polyA_mix_data=$2

# nohup perl $(dirname "$0")/motifFinding.pl --conf $(dirname "$0")/motifConfigure_hg38.xml --bed $bedfile --output ${bedfile%.bed}_findMotif_out --pre $(basename ${bedfile%.bed}) > ${bedfile%.bed}.log 2>&1 &
perl $(dirname "$0")/motifFinding.pl --conf $(dirname "$0")/motifConfigure_hg38.xml --bed $bedfile --output ${bedfile%.bed}_findMotif_out --pre $(basename ${bedfile%.bed})

echo "visualizing by motifStack"
Rscript $(dirname "$0")/motifStack_vis.r -f $total_polyA_mix_data -g ${bedfile%.bed}_findMotif_out/full_gene_HomerGenomeResults/$(basename ${bedfile%.bed}).full_gene.significant.positionSpecifc.motifs.sequences.txt -i ${bedfile%.bed}_findMotif_out/full_gene_HomerGenomeResults/$(basename ${bedfile%.bed}).full_gene.significant.positionSpecifc.motifs.txt -o ${bedfile%.bed}_findMotif_out/full_gene_HomerGenomeResults/$(basename ${bedfile%.bed})

# file_list=$(ls -1v ${bedfile%.bed}_findMotif_out/full_gene_HomerGenomeResults/$(basename ${bedfile%.bed})*_all.pdf|tr '\n' ' ')
# mutool merge -o ${bedfile%.bed}_findMotif_out/full_gene_HomerGenomeResults/$(basename ${bedfile%.bed})_all_motif.pdf $file_list #&> /dev/null
# mupdf-x11 ${bedfile%.bed}_findMotif_out/full_gene_HomerGenomeResults/$(basename ${bedfile%.bed})_all_motif.pdf &> /dev/null

# metagene_list=$(ls -1v ${bedfile%.bed}_findMotif_out/full_gene_HomerGenomeResults/$(basename ${bedfile%.bed})*_mRNA_metagene.pdf|tr '\n' ' ')
# mutool merge -o ${bedfile%.bed}_findMotif_out/full_gene_HomerGenomeResults/$(basename ${bedfile%.bed})_all_metagene.pdf $metagene_list #&> /dev/null
# mupdf-x11 ${bedfile%.bed}_findMotif_out/full_gene_HomerGenomeResults/$(basename ${bedfile%.bed})_all_metagene.pdf &> /dev/null

echo -e "Finished: findMotif done!\n"