#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
library(dplyr)

# change this to the folder with the R objects
setwd("path/to/objects/")

# read in object
# change for other time points
upInKO_object <- readRDS("~/Downloads/makeNetworkPlots/D5_KOvWT_byTime_upInKO_clusterProfResult_forJesse.rds")

# look at top terms
data.frame(upInKO_object) %>% filter(qvalue < 0.05) %>% nrow()

head(upInKO_object,15)

# set number of clusters
number_of_clusters <- 2

# make figure with with label
set.seed(1001)
p_temp <- emapplot(upInKO_object, 
                   cex.params = list(category_label = 0.6), 
                   cluster.params = list(cluster = TRUE, 
                                         n = number_of_clusters),
                   repel = T,
                   #node_label = "none", # uncomment this to remove labels in fugure
                   showCategory = 15 # number of terms to show
                   ) 

minx <- min(p_temp$data$x) * .7
maxx <- max(p_temp$data$x) * 1.3
miny <- min(p_temp$data$y) * .7
maxy <- max(p_temp$data$y) * 1.2
p_temp$coordinates$limits$x <- c(minx, maxx)
p_temp$coordinates$limits$y <- c(miny, maxy)
p_temp

# change name of the output file  depending on which time point is used
pdf(paste0("D5_KOvWT_byTime_upInKO_networkPlot_labels", number_of_clusters, "clusters.pdf"), width = 7, height = 5)
print(p_temp)
dev.off()


# change for other time points
upInKO_object <- readRDS("~/Downloads/forJesse_makeNetworkPlots/D3_KOvWT_byTime_upInKO_clusterProfResult_forJesse.rds")

# look at top terms
data.frame(upInKO_object) %>% filter(qvalue < 0.05) %>% nrow()

head(upInKO_object,15)

# set number of clusters
number_of_clusters <- 3

# make figure with with label
set.seed(1001)
p_temp <- emapplot(upInKO_object, 
                   cex.params = list(category_label = 0.6), 
                   cluster.params = list(cluster = TRUE, 
                                         n = number_of_clusters),
                   repel = T,
                   node_label = "none", # uncomment this to remove labels in fugure
                   showCategory = 15 # number of terms to show
) 

minx <- min(p_temp$data$x) * 1.3
maxx <- max(p_temp$data$x) * 1.3
miny <- min(p_temp$data$y) * 1.3
maxy <- max(p_temp$data$y) * 1.3
p_temp$coordinates$limits$x <- c(minx, maxx)
p_temp$coordinates$limits$y <- c(miny, maxy)
p_temp

# change name of the output file  depending on which time point is used
pdf(paste0("D3_KOvWT_byTime_upInKO_networkPlot", number_of_clusters, "clusters.pdf"), width = 7, height = 5)
print(p_temp)
dev.off()

######################  up in WT ###################### 

# change this to the folder with the R objects
setwd("path/to/objects/")

# read in object
# change for other time points
upInWT_object.d3 <- readRDS("/Users/jessenovak/Dropbox (Dropbox @RU)/20230918_JN_singleCell (3)/rmPhgdhInCKO/DEG_GSEA_rmPhgdhInCKO_db/DEG_GSEA_KOvsWT_byTime/GSVA_top_DEpathways/forJesse_makeNetworkPlots/D3_KOvWT_byTime_upInWT_clusterProfResult_forJesse.rds")
# look at top terms
data.frame(upInWT_object.d3) %>% filter(qvalue < 0.05) %>% nrow()
upInWT_object.d3
head(upInWT_object,15)

# set number of clusters
number_of_clusters <- 3

# make figure with with label
set.seed(1001)
p_temp <- emapplot(upInWT_object.d3, 
                   cex.params = list(category_label = 0.6), 
                   cluster.params = list(cluster = TRUE, 
                                         n = number_of_clusters),
                   repel = T,
                   node_label = "none", # uncomment this to remove labels in fugure
                   showCategory = 15 # number of terms to show
) 

minx <- min(p_temp$data$x) * 1.3
maxx <- max(p_temp$data$x) * 1.3
miny <- min(p_temp$data$y) * 1.3
maxy <- max(p_temp$data$y) * 1.3
p_temp$coordinates$limits$x <- c(minx, maxx)
p_temp$coordinates$limits$y <- c(miny, maxy)
p_temp

upInWT_object.d5 <- readRDS("/Users/jessenovak/Dropbox (Dropbox @RU)/20230918_JN_singleCell (3)/rmPhgdhInCKO/DEG_GSEA_rmPhgdhInCKO_db/DEG_GSEA_KOvsWT_byTime/GSVA_top_DEpathways/forJesse_makeNetworkPlots/D5_KOvWT_byTime_upInWT_clusterProfResult_forJesse.rds")
# look at top terms
data.frame(upInWT_object.d5) %>% filter(qvalue < 0.05) %>% nrow()
upInWT_object.d5
head(upInWT_object,15)

# set number of clusters
number_of_clusters <- 4

# make figure with with label
set.seed(1001)
p_temp <- emapplot(upInWT_object.d5, 
                   cex.params = list(category_label = 0.6), 
                   cluster.params = list(cluster = TRUE, 
                                         n = number_of_clusters),
                   repel = T,
                   #node_label = "none", # uncomment this to remove labels in fugure
                   showCategory = 15 # number of terms to show
) 

minx <- min(p_temp$data$x) * 1.3
maxx <- max(p_temp$data$x) * 1.3
miny <- min(p_temp$data$y) * 1.3
maxy <- max(p_temp$data$y) * 1.3
p_temp$coordinates$limits$x <- c(minx, maxx)
p_temp$coordinates$limits$y <- c(miny, maxy)
p_temp


# change name of the output file  depending on which time point is used
pdf(paste0("D5_KOvWT_byTime_upInKO_networkPlot_labels", number_of_clusters, "clusters.pdf"), width = 7, height = 5)
print(p_temp)
dev.off()


# change for other time points
upInWT_object <- readRDS("~/Downloads/forJesse_makeNetworkPlots/D3_KOvWT_byTime_upInKO_clusterProfResult_forJesse.rds")

# look at top terms
data.frame(upInWT_object) %>% filter(qvalue < 0.05) %>% nrow()

head(upInWT_object,15)

# set number of clusters
number_of_clusters <- 3

# make figure with with label
set.seed(1001)
p_temp <- emapplot(upInWT_object, 
                   cex.params = list(category_label = 0.6), 
                   cluster.params = list(cluster = TRUE, 
                                         n = number_of_clusters),
                   repel = T,
                   node_label = "none", # uncomment this to remove labels in fugure
                   showCategory = 15 # number of terms to show
) 

minx <- min(p_temp$data$x) * 1.3
maxx <- max(p_temp$data$x) * 1.3
miny <- min(p_temp$data$y) * 1.3
maxy <- max(p_temp$data$y) * 1.3
p_temp$coordinates$limits$x <- c(minx, maxx)
p_temp$coordinates$limits$y <- c(miny, maxy)
p_temp

# change name of the output file  depending on which time point is used
pdf(paste0("D3_KOvWT_byTime_upInKO_networkPlot", number_of_clusters, "clusters.pdf"), width = 7, height = 5)
print(p_temp)
dev.off()




