#!/usr/bin/Rscript
##--------------------------------------------###
##--------------------------------------------###
##--------------------------------------------###
##--------------------------------------------###
#install.packages("qtl")

setwd("C:/Users/megah/Documents/Investigacion/Proyectos Jose Luis/Retitinits/Gene_Revision/Haplotype_Flow")

library(qtl)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(dplyr)

LOD_name= "RP1.pruned_w3.QTL.Graphic.All"
LOD_File = read.table(file = LOD_name,
                      header=FALSE, sep="\t", stringsAsFactors = FALSE, fill = T,
                      colClasses = c("numeric","numeric","numeric") )
colnames(LOD_File) <- c("chr",  "pos", "lod")


qtl_plot <- function(input,              # data frame input from scanone
                     mult.pheno = FALSE, # multiple phenotypes?
                     model = "normal",   # model used in scanone
                     chrs = NA,          # chromosomes to display
                     lod = NA,           # LOD threshold
                     rug = FALSE,        # plot marker positions as rug?
                     ncol = NA,          # number of columns for facetting
                     labels = NA         # optional dataframe to plot QTL labels
) {
  
  # if we have multiple phenotypes and/or a 2part model, gather input
  if (mult.pheno & model == "2part") {
    input <- gather(input, group, lod, grep("pheno", colnames(input)))
  } else if (mult.pheno) {
    input <- gather(input, group, lod, grep("pheno", colnames(input)))
  } else if (model == "2part") {
    input <- gather(input, method, lod, lod.p.mu:lod.mu)
  }
  
  # if not all chromosomes should be displayed, subset input
  if (!is.na(chrs)[1]) {
    input <- input[as.character(input$chr) %in% chrs, ]
  }
  
  # if there is more than one LOD column, gather input
  if (!any(colnames(input) == "lod")) {
    input$lod <- input[, grep("lod", colnames(input))]
  }
  
  # if no number of columns for facetting is defined, plot all in one row
  if (is.na(ncol)) {
    ncol <- length(unique(input$chr))
  }
  
  # if labels are set and there is no name column, set from rownames
  if (!is.na(labels)[1]) {
    if (is.null(labels$name)) {
      labels$name <- rownames(labels)
    }
  }
  
  # plot input data frame position and LOD score
  plot <- ggplot(input, aes(x = pos, y = lod)) + {
    
    # if LOD threshold is given, plot as horizontal line
    if (!is.na(lod)[1] & length(lod) == 1) geom_hline(yintercept = lod, linetype = "dashed")
  } + {
    
    if (!is.na(lod)[1] & length(lod) > 1) geom_hline(data = lod, aes(yintercept = lod, linetype = group))
  } + {
    
    # plot rug on bottom, if TRUE
    if (rug) geom_rug(size = 0.1, sides = "b")
  } + {
    
    # if input has column method but not group, plot line and color by method
    if (!is.null(input$method) & is.null(input$group)) geom_line(aes(color = method), size = 2, alpha = 0.6)
  } + {
    
    # if input has column group but not method, plot line and color by group
    if (!is.null(input$group) & is.null(input$method)) geom_line(aes(color = group), size = 2, alpha = 0.6)
  } + {
    
    # if input has columns method and group, plot line and color by method & linetype by group
    if (!is.null(input$group) & !is.null(input$method)) geom_line(aes(color = method, linetype = group), size = 2, alpha = 0.6)
  } + {
    
    # set linetype, if input has columns method and group
    if (!is.null(input$group) & !is.null(input$method)) scale_linetype_manual(values = c("solid", "twodash", "dotted"))
  } + {
    
    # if input has neither columns method nor group, plot black line
    if (is.null(input$group) & is.null(input$method)) geom_line(size = 2, alpha = 0.6)
  } + {
    
    # if QTL positions are given in labels df, plot as point...
    if (!is.na(labels)[1]) geom_point(data = labels, aes(x = pos, y = lod))
  } + {
    
    # ... and plot name as text with ggrepel to avoid overlapping
    if (!is.na(labels)[1]) geom_text_repel(data = labels, aes(x = pos, y = lod, label = name), nudge_y = 0.5)
  } + 
    # facet by chromosome
    facet_wrap(~ chr, ncol = ncol, scales = "free_x") +
    # minimal plotting theme
    theme_minimal() +
    # increase strip title size
    theme(strip.text = element_text(face = "bold", size = 12)) +
    # use RcolorBrewer palette
    scale_color_brewer(palette = "Set1") +
    # Change plot labels
    labs(x = "Chromosome",
         y = "LOD",
         color = "",
         linetype = "")
  
  print(plot)
}

#Genes retrieve function
# specify the database
#ensembl = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

#d=LOD_File[,c(1,2,2,3)]
#colnames(d) <- c("chr",  "start","end", "lod")
d$start <- (d$start*1000000) 
d$end <- (d$end*1000000) 

# loop through rows, get genes, then paste with collapse,
# and finally bind back with data d.
MarkersWithGenes <- cbind(
  d,
  genes = apply(d, 1, function(i){
    x <- getBM(attributes=c("external_gene_name"), 
               filters = c("chromosome_name" , "start", "end"), 
               values = list(i[1], i[2], i[3]), 
               mart = ensembl)
    # keeping only 3 genes, as output is too long.
    # In your case remove below line
    x <- head(x, 3)
    
    # return genes, comma separated
    paste(x$external_gene_name, collapse = ",")
  })
)


##Plotting with qtl_plot function starts
qtl_plot(input = LOD_File)
qtl_plot(input = LOD_File,
         chrs = c(8), 
         rug = TRUE)

qtl_plot(input = LOD_File,
         chrs = c(3,8), 
         rug = TRUE)

  