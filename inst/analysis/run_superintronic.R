#!/usr/bin/env R

# Author Sean Maden
# Run superintronic for gene of interest.

library(superintronic); library(plyranges); library(Rsamtools); library(ggplot2)

# get bam filepaths
gene.name <- "ZRSR2"; gene.dpath <- gene.name; bam.fnv <- list.files(gene.dpath)
bam.fnv <- bam.fnv[grepl("\\.sorted.bam", bam.fnv)]; head(bam.fnv)

# test bams
bam <- scanBam(file.path(gene.dpath, "ZRSR2_SRR5009377.sorted.bam"))

# get bam coverages
bamdf <- data.frame(bam=file.path(gene.dpath, list.files(gene.dpath)), 
    name = list.files(gene.dpath), stringsAsFactors=FALSE)
cvg <- superintronic::compute_coverage_long(bamdf, source = "bam")

# get features from gtf
gtf.path <- "gencode.v35.primary_assembly.annotation.gtf"
features <- gtf.path %>% 
  collect_parts() %>% 
  filter(gene_name == gene.name)

# get overlapping gene parts
cvg_over_features <- cvg %>% 
  select(-bam) %>% 
  join_parts(features)
cvg_over_features

#------------------------------
# coverage and gene model plots
#------------------------------

p <- cvg_over_features %>% 
  mutate(strand = feature_strand) %>% 
  view_coverage(score = score, colour = feature_type) + 
  scale_color_brewer(palette = "Dark2") +
  guides(colour = FALSE) +
  labs(title = "Coverage over SRM")
p

gene_track <- view_segments(unnest_parts(features), 
                            colour = feature_type)
gene_track

p / gene_track

p <- cvg_over_features %>% 
  mutate(strand = feature_strand) %>% 
  view_coverage(score = score, 
                colour = feature_type, 
                facets = vars(name)) + 
  scale_color_brewer(palette = "Dark2") +
  guides(colour = FALSE) +
  labs(title = "Coverage over SRM")
p / gene_track + patchwork::plot_layout(heights = c(3, 1))


