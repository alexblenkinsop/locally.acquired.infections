library(tidyverse)
library(phyloscannerR)
library(ggtree)
library(wesanderson)
library(pals)
library(glue)

outdir <- "~/git/locally.acquired.infections/data/trees"
trsm <- 'HSX'
#trsm <- 'MSM'

if(trsm=='MSM'){
  rdas <- list.files(indir,pattern = "*AmsMSM__redacted.rda")
}else{
  rdas <- list.files(indir,pattern = "*AmsHSX__redacted.rda")
}


#walk(rdas, function(a.rda){
for(i in 1:length(rdas)){
  
	a.rda <- rdas[i]
  load(file.path(indir,a.rda))
  
  subtype <- str_match(a.rda, ".*_subtype_([A-Za-z0-9]+)+_wOutgroup_*")[,2]
  print(subtype)
  ptree <- phyloscanner.trees[[1]]
  
  tree <- ptree$tree
  
  read.counts <- ptree$read.counts[1:6487]
  
  attr(tree, 'BLACKLISTED') <- c(is.na(read.counts), rep(FALSE, tree$Nnode))
  
  read.counts <- c(read.counts, rep(1, tree$Nnode))
  read.counts[which(is.na(read.counts))] <- 1
  
  attr(tree, 'READ_COUNT') <- read.counts
  
  new.branch.colours <- attr(tree, "BRANCH_COLOURS")
  new.individual <- attr(tree, "INDIVIDUAL")
  if(trsm=='MSM'){
    new.branch.colours <- factor(new.branch.colours, levels = c("AmsMSM",  
                                    "AmsnonMSM",
                                    "NL",
                                    "WEurope",
                                    "EEuropeCentralAsia",
                                    "NorthAm",
                                    "LaAmCar",
                                    "Africa",
                                    "MENA",
                                    "Asia",
                                    "Oceania",
                                    NA))
    new.individual <- factor(new.individual, levels =  c("AmsMSM",  
                                                         "AmsnonMSM",
                                                         "NL",
                                                         "WEurope",
                                                         "EEuropeCentralAsia",
                                                         "NorthAm",
                                                         "LaAmCar",
                                                         "Africa",
                                                         "MENA",
                                                         "Asia",
                                                         "Oceania",
                                                         NA))
  }else{
  new.branch.colours <- factor(new.branch.colours, levels = c("AmsHSX",  
                                                                "AmsnonHSX",
                                                                "NL",
                                                                "WEurope",
                                                                "EEuropeCentralAsia",
                                                                "NorthAm",
                                                                "LaAmCar",
                                                                "Africa",
                                                                "MENA",
                                                                "Asia",
                                                                "Oceania",
                                                                NA))
  new.individual <- factor(new.individual, levels =  c("AmsHSX",  
                                                       "AmsnonHSX",
                                                         "NL",
                                                         "WEurope",
                                                         "EEuropeCentralAsia",
                                                         "NorthAm",
                                                         "LaAmCar",
                                                         "Africa",
                                                         "MENA",
                                                         "Asia",
                                                         "Oceania",
                                                         NA))
}
  
  huh <- glasbey(n=21)
  
  huh2 <- c(huh[c(6, 16, 7, 8, 9, 11, 12, 21, 5)], "grey50")
  if(trsm=='MSM'){
    huh2 <- c(rev(pal_npg("nrc")(8))[c(1:6,8,7)], huh[c(7,16,21,19)], "grey50") #HSX
  }else{
    huh2 <- c(pal_npg("nrc")(7), huh[c(7,16,21,19)], "grey50") #HSX
    huh2 <- c(rev(pal_npg("nrc")(8))[c(1:6,8,7)], huh[c(7,16,21,19)], "grey50") #HSX
  }
  
  tips <- tree$tip.label %>% length()
  
  if(trsm=='MSM'){
    tree.display <- ggtree(tree, aes(color=new.branch.colours), size = 0.1) +
      geom_point2(aes(color=new.individual), size=0.66) +
      geom_treescale(x=0.01, y = length(tree$tip.label)*0.5,  offset=4) +
      scale_colour_manual(values = huh2, na.value = "black", drop = F, labels = c("Amsterdam - MSM",
                                                                                  "Amsterdam - non-MSM",
                                                                                  "Netherlands",
                                                                                  "Western Europe",
                                                                                  "Eastern Europe & Central Asia",
                                                                                  "North America",
                                                                                  "Latin America & Caribbean",
                                                                                  "Sub-Saharan Africa",
                                                                                  "Middle-East and North Africa",
                                                                                  "South- and Southeast Asia",
                                                                                  "Oceania",
                                                                                  "Unassigned"), name = "Region/risk group") + 
      theme(legend.text = element_text(size=6)) +
      guides(col=guide_legend(ncol=2)) +
      ylim(-1, length(tree$tip.label) +1)
    if(subtype=='Bc4'){
      tree.display <- tree.display + geom_treescale(x=0.01, y = length(tree$tip.label)*0.5,  offset=-150)
    }
  }else{
    tree.display <- ggtree(tree, aes(color=new.branch.colours), size = 0.1) +
      geom_point2(aes(color=new.individual), size=0.66) +
      geom_treescale(x=0.01, y = length(tree$tip.label)*0.5,  offset=4) +
      scale_colour_manual(values = huh2, na.value = "black", drop = F, labels = c("Amsterdam - HSX",
                                                                                  "Amsterdam - non-HSX",
                                                                                  "Netherlands",
                                                                                  "Western Europe",
                                                                                  "Eastern Europe & Central Asia",
                                                                                  "North America",
                                                                                  "Latin America & Caribbean",
                                                                                  "Sub-Saharan Africa",
                                                                                  "Middle-East and North Africa",
                                                                                  "South- and Southeast Asia",
                                                                                  "Oceania",
                                                                                  "Unassigned"), name = "Region/risk group") + 
      theme(legend.text = element_text(size=6)) +
      guides(col=guide_legend(ncol=2)) +
      ylim(-1, length(tree$tip.label) +1)
  }
  
  ggsave(file.path(outdir,glue("type_{subtype}_tree_newlabels_{trsm}.pdf")), width = 7, height = tips*0.02, limitsize = F)
  
  if(trsm=='MSM'){
    tree.display <- ggtree(tree, aes(color=new.branch.colours), size = 0.1,layout="fan", open.angle=60) +
    geom_point2(aes(color=new.individual), size=0.66) +
      geom_treescale(x=-0.01, y = length(tree$tip.label)*0.5,  offset=-150) + # 02ag
      scale_colour_manual(values = huh2, na.value = "black", drop = F, labels = c("Amsterdam - MSM",
                                                                                "Amsterdam - non-MSM",
                                                                                "Netherlands",
                                                                                "Western Europe",
                                                                                "Eastern Europe & Central Asia",
                                                                                "North America",
                                                                                "Latin America & Caribbean",
                                                                                "Sub-Saharan Africa",
                                                                                "Middle-East and North Africa",
                                                                                "South- and Southeast Asia",
                                                                                "Oceania",
                                                                                "Unassigned"), name = "Region/risk group") + 
    theme(legend.text = element_text(size=6),
          legend.position="bottom",
          plot.margin=unit(c(-45,-20,0,-35), "mm")) +
    guides(col=guide_legend(ncol=2)) +
    ylim(0, length(tree$tip.label))
    
    if(subtype=='02AG') tree.display <- tree.display + theme(plot.margin=unit(c(-15,-20,0,-35), "mm"))
    if(subtype=='Bc2') tree.display <- tree.display + theme(plot.margin=unit(c(-5,-20,0,-35), "mm"))
    if(subtype=='Bc4') tree.display <- tree.display + theme(plot.margin=unit(c(-20,-20,0,-35), "mm"))
    if(subtype=='C') tree.display <- tree.display + theme(plot.margin=unit(c(-18,-20,0,-35), "mm"))
    if(subtype=='D') tree.display <- tree.display + geom_treescale(x=0.01, y = length(tree$tip.label)*0.5,  offset=-25) +
                                        theme(plot.margin=unit(c(-20,-20,0,-35), "mm"))
    if(subtype=='G') tree.display <- tree.display + geom_treescale(x=0.01, y = length(tree$tip.label)*0.5,  offset=-25)

  }else{
    tree.display <- ggtree(tree, aes(color=new.branch.colours), size = 0.1,layout="fan") +
      geom_point2(aes(color=new.individual), size=0.66) +
      geom_treescale(x=0.01, y = length(tree$tip.label)*0.5,  offset=-150) +
      scale_colour_manual(values = huh2, na.value = "black", drop = F, labels = c("Amsterdam - HSX",
                                                                                  "Amsterdam - non-HSX",
                                                                                  "Netherlands",
                                                                                  "Western Europe",
                                                                                  "Eastern Europe & Central Asia",
                                                                                  "North America",
                                                                                  "Latin America & Caribbean",
                                                                                  "Sub-Saharan Africa",
                                                                                  "Middle-East and North Africa",
                                                                                  "South- and Southeast Asia",
                                                                                  "Oceania",
                                                                                  "Unassigned"), name = "Region/risk group") + 
      theme(legend.text = element_text(size=6),
            legend.position="bottom") +
      guides(col=guide_legend(ncol=2)) +
      ylim(-1, length(tree$tip.label) +1)
  }
  
  ggsave(glue("type_{subtype}_tree_newlabels_fan_sm_{trsm}.pdf"), width = tips/130, height = tips/130, limitsize = T)
  if(grepl('02AG',subtype)) ggsave(glue("type_{subtype}_tree_newlabels_fan_sm_{trsm}.pdf"), width = tips/200, height = tips/200, limitsize = T)
  if(grepl('Bc1',subtype)) ggsave(glue("type_{subtype}_tree_newlabels_fan_sm_{trsm}.pdf"), width = tips/500, height = tips/500, limitsize = T)
  if(grepl('Bc2',subtype)) ggsave(glue("type_{subtype}_tree_newlabels_fan_sm_{trsm}.pdf"), width = tips/600, height = tips/600, limitsize = T)
  if(grepl('Bc3',subtype)) ggsave(glue("type_{subtype}_tree_newlabels_fan_sm_{trsm}.pdf"), width = tips/600, height = tips/600, limitsize = T)
  if(grepl('Bc',subtype)) ggsave(glue("type_{subtype}_tree_newlabels_fan_sm_{trsm}.pdf"), width = tips/300, height = tips/300, limitsize = T)
  if(grepl('C',subtype)) ggsave(glue("type_{subtype}_tree_newlabels_fan_sm_{trsm}.pdf"), width = tips/200, height = tips/200, limitsize = T)
  if(grepl('C',subtype)) ggsave(glue("type_{subtype}_tree_newlabels_fan_sm_{trsm}.pdf"), width = tips/300, height = tips/300, limitsize = T)
  if(grepl('D',subtype)) ggsave(glue("type_{subtype}_tree_newlabels_fan_sm_{trsm}.pdf"), width = tips/50, height = tips/50, limitsize = T)
  if(grepl('G',subtype)) ggsave(glue("type_{subtype}_tree_newlabels_fan_sm_{trsm}.pdf"), width = tips/80, height = tips/80, limitsize = T)
  if(subtype=='06cpx') ggsave(glue("type_{subtype}_tree_newlabels_fan_sm_{trsm}.pdf"), width = tips/20, height = tips/20, limitsize = T)
  if(subtype=='Bc4') ggsave(glue("type_{subtype}_tree_newlabels_fan_sm_{trsm}.pdf"), width = tips/20, height = tips/20, limitsize = T)
  #ggsave(glue("type_{subtype}_tree_newlabels_fan_sm_{trsm}.pdf"), width = tips/200, height = tips/200, limitsize = F)
}




if(any(read.counts!=1 & !is.na(read.counts))){
  tree.display <- tree.display + geom_point2(aes(color=INDIVIDUAL, size=READ_COUNT, shape=BLACKLISTED), alpha=0.5)
}

x.max <- ggplot_build(tree.display)$layout$panel_scales_x[[1]]$range$range[2]

# for ggplot2 version compatibility

if(is.null(x.max)){
  x.max <- ggplot_build(tree.display)$layout$panel_params[[1]]$x.range[2]
}

if(is.null(x.max)){
  x.max <- ggplot_build(tree.display)$layout$panel_ranges[[1]]$x.range[2]
}


tree.display <- tree.display + ggplot2::xlim(0, 1.1*x.max)
tree.display
file.name <- glue("type_{subtype}_tree_newlabels_smalltips.pdf")
ggsave(file.name, device="pdf", height = pdf.hm*length(tree$tip.label), width = pdf.w, limitsize = F)

