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

  input:
  path fastqs
    
  output:
  path "*_fastqc.*", emit: fastqc_dir

  """
  fastqc -t 4 -o . ./*
  """
}


process run
{
  publishDir 'nf_out/bams', pattern: "*.bam", mode: 'copy'

  beforeScript 'source /etc/bashrc; module load star/2.6.0c'

  errorStrategy 'retry'
  maxRetries 3

  executor "sge"
  cpus 8
  clusterOptions "-P bf528 -pe omp 16"

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
}
