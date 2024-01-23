#!/usr/bin/env Rscript
options(warn = -1)
suppressMessages(library("optparse"))
suppressMessages(library("openxlsx"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggplot2"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("optparse"))
suppressMessages(library("motifStack"))
suppressMessages(library("gridGraphics"))
suppressMessages(library("stringr"))


option_list = list(
  make_option(c("-f", "--infile1"), type = "character", default = NULL,
              help = "input mixed data of totalRNA and polyARNA psiFinder information [file]", metavar = "character"),
  make_option(c("-g", "--infile2"), type = "character", default = NULL,
              help = "findMotif_pssm.file [file]", metavar = "character"),
  make_option(c("-i", "--infile3"), type = "character", default = NULL,
              help = "findMotif_pssm.file2 [file]", metavar = "character"),
  make_option(c("-o", "--outfile"), type = "character", default = NULL,
              help = "output file name [default= %default]", metavar = "character")
);
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if (is.null(opt$infile1) || is.null(opt$infile2) || is.null(opt$infile3) || is.null(opt$outfile)) {
  print_help(opt_parser);
  stop("Please provide -f infile1 (Day0_common_rep1_anno_group_redundance.xlsx) -g infile2 (totalRNA_Day0_common_motifs.full_gene.significant.positionSpecifc.motifs.sequences.txt) -i infile3 (ePSI_seq_total_polyA_Day0_mix_motifs.full_gene.significant.positionSpecifc.motifs.txt) and -o outfile (Day0_common_rep1_findMotif_pssm) option", call. = FALSE);
}
input.mixed.file = opt$infile1
findMotif_pssm.file = opt$infile2
findMotif_pssm.file2 = opt$infile3
outFile = opt$outfile

print(input.mixed.file)
print(findMotif_pssm.file)
print(findMotif_pssm.file2)
print(outFile)

input <- read.delim(input.mixed.file, header = T)
positionSpecifc_motifs_sequences <- read.delim(findMotif_pssm.file, header = T)
colnames(positionSpecifc_motifs_sequences)[1] <- "motifName"
positionSpecifc_motifs_sequences$motif_col <- paste(positionSpecifc_motifs_sequences$motifName, positionSpecifc_motifs_sequences$motifSeq, positionSpecifc_motifs_sequences$distType, positionSpecifc_motifs_sequences$dist, sep = "_")
positionSpecifc_motifs_pval <- read.delim(findMotif_pssm.file2, header = T)
colnames(positionSpecifc_motifs_pval)[1] <- "motifName"
positionSpecifc_motifs_pval$motif_col <- paste(positionSpecifc_motifs_pval$motifName, positionSpecifc_motifs_pval$motifSeq, positionSpecifc_motifs_pval$distType, positionSpecifc_motifs_pval$dist, sep = "_")

positionSpecifc_motifs_sequences <- positionSpecifc_motifs_sequences %>% left_join(positionSpecifc_motifs_pval, by = c("motif_col" = "motif_col"))
positionSpecifc_motifs_sequences <- positionSpecifc_motifs_sequences %>% dplyr::select(-c(11:14))
colnames(positionSpecifc_motifs_sequences)[1:4] <- str_replace(colnames(positionSpecifc_motifs_sequences)[1:4], ".x", "")
colnames(positionSpecifc_motifs_sequences)[11] <- "negLg_pval"

positionSpecifc_motifs_sequences_append_info <- positionSpecifc_motifs_sequences %>% left_join(input, by = c("seqName" = "name"))
sncRNA_index <- which(positionSpecifc_motifs_sequences_append_info$gene_biotype == "sncRNA")
positionSpecifc_motifs_sequences_append_info$gene_biotype[sncRNA_index] <- positionSpecifc_motifs_sequences_append_info$gene_feature[sncRNA_index]
positionSpecifc_motifs_sequences_append_info <- positionSpecifc_motifs_sequences_append_info %>% arrange(negLg_pval)
write.xlsx(positionSpecifc_motifs_sequences_append_info, paste(outFile, "_findMotif_pssm_append_info.xlsx", sep = ""), overwrite = TRUE)
positionSpecifc_motifs_sequences_append_info$motif_col <- factor(positionSpecifc_motifs_sequences_append_info$motif_col, levels = unique(positionSpecifc_motifs_sequences_append_info$motif_col))
motifName_list <- split(positionSpecifc_motifs_sequences_append_info, positionSpecifc_motifs_sequences_append_info$motif_col)

lapply(seq_along(motifName_list), function(x) {
  if(any(unique(motifName_list[[x]]$gene_biotype) %in% "tRNA")){
    print(paste(names(motifName_list[x])," found tRNA gene_biotype",sep=""))
    # write.xlsx(motifName_list[[x]], paste(outFile, "_", names(motifName_list)[x], ".xlsx", sep = ""), overwrite = TRUE)
  }else{
    print(paste(names(motifName_list[x])," found no tRNA gene_biotype",sep=""))
  }
})


invisible(lapply(seq_along(motifName_list), function(x) {
  print(paste(paste("loopx ", x, ":", sep = ""), " processing ", paste(dirname(outFile),"/homerResults/", str_replace_all(names(motifName_list)[x], "_[A-Z]*_[0-9]*pDist_[0-9]*", ""), sep = ""), sep = ""))
  #homerResults

  if (str_detect(motifName_list[[x]]$motifInSeq[1], "T")) {
    if (strsplit(motifName_list[[x]]$sequence[1], split = "")[[1]][11] == "T") {
      file <- paste(dirname(outFile),"/homerResults/", str_replace_all(names(motifName_list)[x], "_[A-Z]*_[0-9]*pDist_[0-9]*", ""), sep = "")
      motiffile <- read.table(file, header = F, skip = 1)

      names(motiffile) <- c("A", "C", "G", "U")
      motiffile <- data.frame(motiffile, Y = 0)
      motifInSeq_pos <- str_locate(motifName_list[[x]]$sequence[1], motifName_list[[x]]$motifInSeq[1])[1, 1]

      if (motifInSeq_pos > 11) {
        print(paste("motifInSeq position outrange 11 for ", paste(dirname(outFile),"/homerResults/", str_replace_all(names(motifName_list)[x], "_[A-Z]*_[0-9]*pDist_[0-9]*", ""), sep = ""), sep = ""))
      } else {
        index <- 11 - motifInSeq_pos + 1
        if (is.na(motiffile$Y[index])) {
          print(paste("motifInSeq position U base get NA for ", paste(dirname(outFile),"/homerResults/", str_replace_all(names(motifName_list)[x], "_[A-Z]*_[0-9]*pDist_[0-9]*", ""), sep = ""), sep = ""))
        } else {
          motiffile$Y[index] <- motiffile$U[index]

          if (motiffile$Y[index] > 0.7) {
            motiffile$U[index] <- 0
            A <- motiffile[, 1]
            C <- motiffile[, 2]
            G <- motiffile[, 3]
            U <- motiffile[, 4]
            Y <- motiffile[, 5]
            data <- rbind(A, C, G, U, Y)
            pcm <- data[, 1:ncol(data)]
            rownames(pcm) <- c("A", "C", "G", "U", "Ψ")
            motif <- new("pcm", mat = as.matrix(pcm), name = paste("Negative log(p-value): ", motifName_list[[x]]$negLg_pval[1], sep = ""))

            #gene_biotype pie plot
            motifName_list_tmp <- motifName_list[[x]] %>% distinct(extendSeq, .keep_all = TRUE)
            pseudoU_anno_genetype_num <- arrange(as.data.frame(table(motifName_list_tmp$gene_biotype)), desc(Freq))

            if (sum(pseudoU_anno_genetype_num$Freq) > 15) {
              percent <- round(100 * pseudoU_anno_genetype_num$Freq / sum(pseudoU_anno_genetype_num$Freq), 2)
              percent <- paste('(', percent, "%", ", ", pseudoU_anno_genetype_num$Freq, ')', sep = "")
              pseudoU_anno_genetype_num$Var1 <- paste(pseudoU_anno_genetype_num$Var1, percent, sep = '')
              pseudoU_anno_genetype_num$Var1 <- factor(pseudoU_anno_genetype_num$Var1, levels = pseudoU_anno_genetype_num$Var1)
              mycol = c(brewer.pal(12, "Set3"), brewer.pal(12, "Paired"), brewer.pal(11, "Spectral"), brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(8, "Pastel1"), brewer.pal(8, "Pastel2"), brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"))
              pie_plot <- ggplot(data = pseudoU_anno_genetype_num, mapping = aes(x = 'Content', y = Freq, fill = Var1)) +
                    geom_bar(stat = 'identity', position = 'stack', width = 1) + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') +
                      guides(fill = guide_legend(title = "gene biotype")) +
                      scale_fill_manual(values = mycol) +
                      ggtitle(paste(paste(x, "_", names(motifName_list)[x], ":", sep = ""), motifName_list[[x]]$motifSeq[1], sep = " ")) +
                      theme(
                      plot.title = element_text(hjust = 0.5),
                      panel.grid = element_blank(),
                      panel.background = element_blank(),
                      axis.text = element_blank(),
                      axis.title = element_blank(),
                      axis.ticks = element_blank(),
                      plot.margin = unit(c(2, 2, 4, 2), "cm"),
                      legend.position = "right"
                      )

              motifName_list_tmp <- motifName_list[[x]] %>% distinct(extendSeq, .keep_all = TRUE)
              motifName_list_tmp_mRNA <- motifName_list_tmp[motifName_list_tmp$gene_biotype == "mRNA",]

              if (!is.null(motifName_list_tmp_mRNA)) {
                if (length(motifName_list_tmp_mRNA$gene_biotype) > 2) {
                  motifName_list_tmp_mRNA$seqName <- paste(motifName_list_tmp_mRNA$seqName, motifName_list_tmp_mRNA$gene_feature, sep = "_")
                  mRNA_bed <- motifName_list_tmp_mRNA %>% select(chrom, chromStart, chromEnd, seqName, foldChange.y, strand)
                  write.table(mRNA_bed, paste(outFile, "_", x, "_", names(motifName_list)[x], "_mRNA.bed", sep = ""), row.names = F, col.names = F, quote = F)
                  # system(paste("nohup bash /public/home/chenzr/PSI_Seq_brainCell/A1-A12-totalRNA-result/psiFinder_SVM_res/findMotif_pssm/metagene.sh -i ", paste(outFile, "_", x, "_", names(motifName_list)[x], "_mRNA.bed", sep = ""), "-g /public/home/chenzr/software/qt_project/PseudoTB_dev/build-psiFinder-Desktop_Qt_5_9_1_static_libpng-Release/snakemake/genome/hg38.fa -a /public/home/chenzr/software/qt_project/PseudoTB_dev/build-psiFinder-Desktop_Qt_5_9_1_static_libpng-Release/script/gencode.v32.chr_patch_hapl_scaff.annotation.gtf -A /public/home/chenzr/software/qt_project/PseudoTB_dev/build-psiFinder-Desktop_Qt_5_9_1_static_libpng-Release/script/hg38.genecode.v30.tRNA.snoRNA.miRNA.rmsk.exonFeatures.bed6", "-m", paste(paste(x, "_", names(motifName_list)[x], ":", sep = ""), motifName_list_tmp$motifSeq[1], sep = ""), ">", paste(outFile, "_", names(motifName_list)[x], "_metagene.log", sep = ""), "2>&1 &", sep = " "))
                } else {
                  print("mRNA number lower than 5")
                }
              } else {
                print("No genetype number summary!")
              }

              plot(motif)
              grab <- grid.grab()

              pdf(paste(outFile, "_", x, "_", names(motifName_list)[x], "_", motifName_list[[x]]$motifSeq[x], "_all", ".pdf", sep = ""), width = 7, height = 7)
              print(pie_plot)
              pushViewport(viewport(x = .45, y = .2, height = .25, width = .4))
              grid.draw(grab)
              invisible(dev.off())
            } else {
              print("total gene_biotype number less than 15!")
            }
          } else {
            print("Ψ percentage less than 0.8")
          }



        }
      }

    } else {
      print(paste("position 11 (21 nt) not U base for ", paste(dirname(outFile),"/homerResults/", str_replace_all(names(motifName_list)[x], "_[A-Z]*_[0-9]pDist_[0-9]", ""), sep = ""), sep = ""))
    }
  } else {
    print(paste("motifInSeq not contain U base for ", paste(dirname(outFile),"/homerResults/", str_replace_all(names(motifName_list)[x], "_[A-Z]*_[0-9]pDist_[0-9]", ""), sep = ""), sep = ""))
  }


}))



