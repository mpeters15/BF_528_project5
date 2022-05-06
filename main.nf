#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process getFastqs {
  input:
  val sample_ids_ch

  output:
  path '*fastq.gz', emit: fastqs
  tuple val(sample_ids_ch), path('*_1.fastq.gz'), path('*_2.fastq.gz'), emit: fastq_pairs

  """
  ln -s /project/bf528/project_3/samples/${sample_ids_ch}* .
  """
}


process runFastqc {
  beforeScript 'module load fastqc'

  executor "sge"
  clusterOptions "-P bf528 -pe omp 16"

  input:
  path fastqs
    
  output:
  path "*_fastqc.*", emit: fastqc_dir

  """
  fastqc -t 16 -o . ./*
  """
}


process runSTAR {
  publishDir 'nf_out/bams', pattern: "*.bam"//, mode: 'copy'

  beforeScript 'source /etc/bashrc; module load star/2.6.0c'

  errorStrategy 'retry'
  maxRetries 3

  executor "sge"
  clusterOptions "-P bf528 -pe omp 18"

  input:
  tuple val(sample_id), path(fp1), path(fp2)
    
  output:
  path "${sample_id}_*", emit: star_dir
  tuple val(sample_id), path("*${sample_id}*.bam"), emit: bams

  """
  #!/bin/bash

  STAR --genomeDir /project/bf528/project_3/reference/rn4_STAR/ \
  --readFilesIn ${fp1} ${fp2} \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outFileNamePrefix ${sample_id}_ \
  --outSAMtype BAM SortedByCoordinate
  """
}


process runMultiQC_Fastqc_STAR {
  publishDir 'nf_out/multiQC/FastqcAndStar', mode: 'copy'

  beforeScript 'module purge; module load python3/3.7.9; module load multiqc/1.10.1'

  input:
  path fastqc_dir
  path star_dir
    
  output:
  path "*" , emit: runMultiQC_FeatCounts_dir

  """
  multiqc -o . ${fastqc_dir} ${star_dir} 
  """
}





process runFeatureCounts {
  publishDir 'nf_out/featureCounts/', mode: 'copy'
  
  beforeScript 'module load subread/1.6.2'

  input:
  tuple val(sample_id), path(sample_bam)

  output:
  path "featureCounts_*" , emit: featureCounts_dir
  
  """

  featureCounts \
  -T 16 \
  -a /project/bf528/project_3/reference/rn4_refGene_20180308.gtf \
  -o featureCounts_${sample_id}.txt ${sample_bam}
  """
}


process runMultiQC_FeatCounts {
  publishDir 'nf_out/multiQC/featureCounts', mode: 'copy'

  beforeScript 'module purge; module load python3/3.7.9; module load multiqc/1.10.1'

  input:
  path featureCounts_dir
    
  output:
  path "*" , emit: runMultiQC_FeatCounts_dir

  """
  multiqc -o . ${featureCounts_dir} 
  """
}


process combineFeatureCounts {
  publishDir 'nf_out/results/', mode: 'copy'

  beforeScript 'module load R'

  input:
  path featureCounts_dir
    
  output:
  path "featureCounts_combined.csv", emit: featureCounts_combined

  """
  #!/usr/bin/env Rscript
  library(tidyverse)
  library(dplyr)

  all_files <- list.files(path = "./", pattern="*.txt\$")
  num_files <- length(all_files)

  for (i in 1:num_files){
    my_file <- all_files[i]
    lastcolname <- gsub( ".txt", "",  gsub("featureCounts_", "", all_files[i])  )
    colnames = c('Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length', lastcolname)

    datatable <- read.table(my_file, header=TRUE, col.names=colnames)
    
    # Get first and last columns
    datatable <- datatable[,c(1, ncol(datatable))]
    
    # Convert to tibble and merge all tibbles
    myTibble <- as_tibble(datatable)
    if (i==1){ combinedTibble <- myTibble }
    else{ combinedTibble <- full_join(combinedTibble, myTibble, by="Geneid") }
  }
  
  combinedTibble <- combinedTibble %>% select(Geneid, everything())

  names(combinedTibble) <- sub("Geneid", "gene_id", names(combinedTibble))

  write.csv(combinedTibble, file="featureCounts_combined.csv", row.names=FALSE)
  """
}


process makeBoxPlots {
  publishDir 'nf_out/results/', mode: 'copy'

  beforeScript 'module load R'

  input:
  path featureCounts_combined
    
  output:
  path "*boxplot*"

  """
  #!/usr/bin/env Rscript  
  library(tidyverse)
  library(reshape2)
  library(ggplot2)

  myTibble <- as_tibble(read.csv("featureCounts_combined.csv", header=TRUE))
  colnames(myTibble)
  t1 <- melt(myTibble, id.vars='gene_id', 
          measure.vars=c('SRR1177963', 'SRR1177964', 'SRR1177965', 
                          'SRR1177997', 'SRR1177999', 'SRR1178002', 
                          'SRR1178014', 'SRR1178021', 'SRR1178047'))

  t1\$moa[t1\$variable == "SRR1177963"] <- 'PPARA'
  t1\$moa[t1\$variable == "SRR1177964"] <- 'PPARA'
  t1\$moa[t1\$variable == "SRR1177965"] <- 'PPARA'
  t1\$moa[t1\$variable == "SRR1177997"] <- 'AhR'
  t1\$moa[t1\$variable == "SRR1177999"] <- 'AhR'
  t1\$moa[t1\$variable == "SRR1178002"] <- 'AhR'
  t1\$moa[t1\$variable == "SRR1178014"] <- 'CAR/PXR'
  t1\$moa[t1\$variable == "SRR1178021"] <- 'CAR/PXR'
  t1\$moa[t1\$variable == "SRR1178047"] <- 'CAR/PXR'


  ggplot(t1) +
  geom_boxplot(aes(x=variable, y=value, fill=moa)) +
          scale_y_continuous(trans='log10') +
          theme_classic() +
          labs(
          title = "Distribution of Gene Counts",
          x = "Sample",
          y = "Gene Counts (log10)",
          fill="Mechanism of Action"
          )

  ggsave( paste("boxplot", "jpeg", sep = "."), 
          plot = last_plot(), device = jpeg, width = 10, height = 10, units = 'in', dpi = 300, )
  """
}

process lastPart {
  publishDir 'nf_out/DESEQ', mode: 'copy'

  beforeScript 'module load R'

  input:
  path featureCounts_combined

  output:
  path "*"

  """
  #!/usr/bin/env Rscript

  library(tidyverse)
  library(DESeq2)
  library(ggplot2)


  controls <- read.table('/project/bf528/project_3/samples/control_counts.csv', header = TRUE, sep = ',')
  group_6_rna_info <- read.csv('/project/bf528/project_3/groups/group_6_rna_info.csv', header = TRUE)
  controls_group_6 <- controls[ , group_6_rna_info\$Run[ group_6_rna_info\$mode_of_action=='Control' ] ]
  controls_group_6\$gene_id <- controls\$Geneid
  

  treatment_counts <- read.csv('featureCounts_combined.csv')
  #combine controls with treatments
  counts_combined <- merge(controls_group_6, treatment_counts, by='gene_id') %>% column_to_rownames(var='gene_id')

  # sample information
  group_6_rna_info1 <- group_6_rna_info[ group_6_rna_info\$vehicle == 'CMC_.5_%', ] %>% 
    subset(mode_of_action == 'AhR' | mode_of_action == 'Control')
  dna_subset <- counts_combined[ , group_6_rna_info1\$Run ]

  group_6_rna_info2 <- group_6_rna_info[ group_6_rna_info\$vehicle == 'CORN_OIL_100_%', ] %>%
    subset(mode_of_action == 'CAR/PXR' | mode_of_action == 'Control')
  car_pxr_subset<- counts_combined[ , group_6_rna_info2\$Run ]

  group_6_rna_info3 <- group_6_rna_info[ group_6_rna_info\$vehicle == 'CMC_.5_%' , ] %>% 
    subset(mode_of_action == 'PPARA' | mode_of_action == 'Control')
  ppara_subset<- counts_combined[ , group_6_rna_info3\$Run ]


  DESeqDataSet1 <- DESeqDataSetFromMatrix(
    colData = group_6_rna_info1,
    countData = dna_subset,
    design= ~ mode_of_action
  )

  DESeqDataSet2 <- DESeqDataSetFromMatrix(
    colData = group_6_rna_info2,
    countData = car_pxr_subset,
    design= ~ mode_of_action
  )

  DESeqDataSet3 <- DESeqDataSetFromMatrix(
    colData = group_6_rna_info3,
    countData = ppara_subset,
    design= ~ mode_of_action
  )


  DESeqDataSet1\$mode_of_action <- relevel( DESeqDataSet1\$mode_of_action, ref='Control' )
  DESeqDataSet2\$mode_of_action <- relevel( DESeqDataSet2\$mode_of_action, ref='Control' )
  DESeqDataSet3\$mode_of_action <- relevel( DESeqDataSet3\$mode_of_action, ref='Control' )


  DESeqDataSet1 <- DESeq(DESeqDataSet1)
  DESeqDataSet2 <- DESeq(DESeqDataSet2)
  DESeqDataSet3 <- DESeq(DESeqDataSet3)



  res_AhR <- results(DESeqDataSet1, contrast=c('mode_of_action','AhR','Control'))
  res_AhR <- lfcShrink(DESeqDataSet1, 2)
  write.csv(res_AhR, 'DESeq_AhR_Rresults.csv')
  write.csv(res_AhR[ order(res_AhR\$padj), ], 'AhR_resultsPADJ.csv')
  print("Number of significant genes at threshold for AhR:")
  print( length(which(res_AhR\$padj < 0.05)) )


  res_CAR_PXR <- results(DESeqDataSet2, contrast=c('mode_of_action','CAR/PXR','Control'))
  res_CAR_PXR <- lfcShrink(DESeqDataSet2, 2)
  write.csv(res_CAR_PXR, 'DESeq_CAR_PXR_results.csv')
  write.csv(res_CAR_PXR[ order(res_CAR_PXR\$padj), ], 'CAR_PXR_resultsPADJ.csv')
  print("Number of significant genes at threshold for CAR/PXR:")
  print( length(which(res_CAR_PXR\$padj < 0.05)) )


  res_PPARA <- results(DESeqDataSet3, contrast=c('mode_of_action','PPARA','Control'))
  res_PPARA <- lfcShrink(DESeqDataSet3, 2)
  write.csv(res_PPARA, 'DESeq_PPARA_results.csv')
  write.csv(res_PPARA[ order(res_PPARA\$padj), ], 'PPARA_resultsPADJ.csv')
  print("Number of significant genes at threshold for PPARA:")
  print( length(which(res_PPARA\$padj < 0.05)) )



  res_AhR = res_AhR[ order(res_AhR\$padj), ]
  res_AhR = res_AhR[ complete.cases(res_AhR), ]
  res_CAR_PXR = res_CAR_PXR[ order(res_CAR_PXR\$padj), ]
  res_CAR_PXR = res_CAR_PXR[ complete.cases(res_CAR_PXR), ]
  res_PPARA = res_PPARA[ order(res_PPARA\$padj), ]
  res_PPARA = res_PPARA[ complete.cases(res_PPARA), ]



  AhR_dim = dim(res_AhR[ res_AhR\$padj<0.05, ])[1]
  AhR_top10 = row.names(res_AhR)[1:10]
  CAR_PXR_dim = dim(res_CAR_PXR[ res_CAR_PXR\$padj<0.05, ])[1]
  CAR_PXR_top10 = row.names(res_CAR_PXR)[1:10]
  PPARA_dim = dim(res_PPARA[ res_PPARA\$padj<0.05, ])[1]
  PPARA_top10 = row.names(res_PPARA)[1:10]

  final_df = data.frame(
    col1=c(AhR_dim,AhR_top10),
    col2=c(CAR_PXR_dim,CAR_PXR_top10),
    col3=c(PPARA_dim,PPARA_top10)
    )
    
  names(final_df) = c("AhR", "CAR/PXR", "PPARA")
  rownames(final_df) = c("Number of Significant Genes:", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
  write.csv(final_df, "SignificantGeneCountAndTopTenSignificantGenes.csv")

  # Histograms of FC Values
  graph_df <- as.data.frame(res_AhR) %>% 
    mutate(Group = "AhR") %>% 
      bind_rows(as.data.frame(res_CAR_PXR) %>% 
      mutate(Group = "CAR/PXR"), as.data.frame(res_PPARA) %>% 
        mutate(Group = "PPARA")) %>%
          filter(padj < 0.05)
  
  print(head(graph_df))

  graph_df %>%
    ggplot() + 
    geom_histogram(aes(x = log2FoldChange, fill = Group), bins = 50) + 
    facet_wrap(~ Group) + 
    labs(title = "Log2 Fold-Change of Significant DE Genes", 
        x = "Log2 Fold Change", y = "Count") + 
    theme(legend.position="none")

  ggsave("DEhist.png")

  graph_df %>%
    ggplot() + 
    geom_point(aes(x = log2FoldChange, y = padj, color = Group), alpha = 0.5) + 
    facet_wrap(~ Group) + 
    labs(title = "Log2 Fold-Change vs. Adjusted p-value", x = "Log2 Fold-Change", 
    y = "Adjusted p-value") + 
    theme(legend.position="none")

  ggsave("log2foldchange.png")  

  cts = counts(DESeqDataSet1, normalized=TRUE)
  write.csv(cts, 'AhR_deseq_norm_counts.csv')
  cts = counts(DESeqDataSet2, normalized=TRUE)
  write.csv(cts, 'CAR_PXR_deseq_norm_counts.csv')
  cts = counts(DESeqDataSet3, normalized=TRUE)
  write.csv(cts, 'PPARA_deseq_norm_counts.csv')

  """
}


workflow {
  sample_ids_ch = channel.of('SRR1177997', 'SRR1177999', 'SRR1178002', 'SRR1178014', 'SRR1178021', 'SRR1178047', 'SRR1177963', 'SRR1177964', 'SRR1177965')
  
  getFastqs( sample_ids_ch )

  runFastqc( getFastqs.out.fastqs.flatten() )
  
  runSTAR( getFastqs.out.fastq_pairs )

  runMultiQC_Fastqc_STAR( runFastqc.out.fastqc_dir.collect(), runSTAR.out.star_dir.collect() )

  runFeatureCounts( runSTAR.out.bams )
  
  runMultiQC_FeatCounts( runFeatureCounts.out.featureCounts_dir.collect() )

  combineFeatureCounts( runFeatureCounts.out.featureCounts_dir.collect() )
  
  makeBoxPlots( combineFeatureCounts.out.featureCounts_combined.collect() )
  
  lastPart( combineFeatureCounts.out.featureCounts_combined.collect() )
}
