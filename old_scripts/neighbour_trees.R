# ploting neighbour joining tree (sketch)

diss <- snpgdsIBS(genofile_all, sample.id = pca_samples, num.thread = 3, maf= 0.05, missing.rate = 0.05)
diss_mat <- diss$ibs
row.names(diss_mat) <- diss$sample.id
colnames(diss_mat) <- diss$sample.id

cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

distCor <- function(x) as.dist(1-cor(t(x)))
hclustAvg <- function(x) hclust(x, method="average")

heatmap.3(diss_mat, trace="none", zlim=c(-3,3), reorder=FALSE,
          distfun=distCor, hclustfun=hclustAvg, col=rev(cols), symbreak=FALSE) 


clust <- snpgdsHCluster(diss)
tmp <- snpgdsDrawTree(clust, yaxis.height = FALSE, outlier.n = 3, outlier.col = "red")

snpgdsDrawTree(clust, main="HapMap Phase II",
               edgePar=list(col=rgb(0.5,0.5,0.5, 0.75), t.col="black"))

# males 2012 and 2014

male_df <- meta_df %>%
  filter(sex == "M") %>%
  #filter(region !="CB") %>%
  filter(!(grepl(bad_samples, id))) 

male_ids <- male_df %>%  
  select(id) %>% 
  unlist %>% 
  as.character

#male_ids <- sample(male_ids, 50)

diss_males <- snpgdsDiss(genofile_all, sample.id = male_ids, num.thread = 3, maf= 0.01, missing.rate = 0)

diss_males$sample.id <- diss_males$sample.id %>% gsub("whtstbk_gbs_2012_brds_", "2012_", .)
clust <- snpgdsHCluster(diss_males)
#plot(clust$hclust, hang = -1, cex = 0.6)


# Default plot

nodePar <- list(lab.cex = 0.0, pch = c(NA, 19), 
                cex = 0.5)

edgePar <- list(col = "grey", lwd = 1)

plot(hcd,  xlab = "Height", nodePar = nodePar, 
     edgePar = edgePar, horiz = TRUE, leaflab = "none")

node_cols <- (grepl("2012", labels(hcd)) %>% as.numeric) +1
branch_cols <- as.facotr() %>% as.numeric) +1


hcd <- as.dendrogram(clust$hclust, hang = 0.1)
hcd <- hcd %>% 
  #set("labels", c("")) %>% 
  set("leaves_col", node_cols) %>%
  set("branches_k_color", k=3) %>%
  set("labels_cex", 0.5) %>%
  set("leaves_pch", 19)

d1=color_branches(hcd,5, col = c(3,1,1,4,1))
plot(d1) # selective coloring of branches :)

ggd1 <- as.ggdend(hcd)
ggplot(ggd1) 

ggplot(ggd1, labels = FALSE) + 
  scale_y_reverse(expand = c(0.2, 0)) +
  coord_polar(theta="x")