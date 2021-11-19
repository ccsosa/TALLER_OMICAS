########################################################################################
#cargar librerias
require(gprofiler2);library(biomaRt);library(topGO);
require(clusterProfiler);require(GOsummaries)

########################################################################################
#MULTICATEGORÍA
########################################################################################
url_file = "https://raw.githubusercontent.com/ccsosa/R_Examples/master/Hallmarks_of_Cancer_AT.csv"
x_group <- read.csv(url_file)
#leyendo archivo
x_group[,1] <- NULL
#definiendo categorias
CH <- c("AID","AIM","DCE","ERI","EGS","GIM","IA","RCD","SPS","TPI")

#preparación
x_Hsap <- lapply(seq_len(length(CH)), function(i){
  x_unique <- unique(na.omit(x_group[,i]))
  x_unique <- x_unique[which(x_unique!="")]
  x_unique <- as.list(x_unique)
  return(x_unique)
})

names(x_Hsap) <- CH
#GO
x_s_group <-  gprofiler2::gost(query = x_Hsap,
                               organism = "hsapiens", ordered_query = FALSE,
                               multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                               measure_underrepresentation = FALSE, evcodes = T,
                               user_threshold = 0.05, correction_method = "fdr",
                               domain_scope = "annotated", custom_bg = NULL,
                               numeric_ns = "", sources = "GO:BP", as_short_link = FALSE)

res_group <- x_s_group$result





x_Hsap_2 <- list()
for(i in 1:length(x_Hsap)){
  x_Hsap_2[[i]] <- clusterProfiler::bitr(as.character(unlist(x_Hsap[[i]])),
                                         fromType = "SYMBOL",
                                         toType = c("ENTREZID"),
                                         OrgDb = "org.Hs.eg.db")[,2]
}
names(x_Hsap_2) <-CH
x_compare <- compareCluster(geneClusters=x_Hsap_2,enrichGO, OrgDb = org.Hs.eg.db)

dotplot(x_compare)
dotplot(x_compare, x="group") + facet_grid(~Cluster)
cnetplot(x_compare)

# Modificar archivo de resultados para correr enrichplot
gp_mod =              res_group[,c("query",
                                   "source",
                                   "term_id",
                                   "term_name",
                                   "p_value",
                                   "query_size",
                                   "intersection_size",
                                   "term_size",
                                   "effective_domain_size",
                                   "intersection")]

gp_mod$GeneRatio = paste0(gp_mod$intersection_size, "/", gp_mod$query_size)
gp_mod$BgRatio = paste0(gp_mod$term_size, "/", gp_mod$effective_domain_size)
names(gp_mod) = c("Cluster", "Category", "ID", "Description", "p.adjust",
                  "query_size", "Count", "term_size", "effective_domain_size",
                  "geneID", "GeneRatio", "BgRatio")
gp_mod$geneID = gsub(",", "/", gp_mod$geneID)
gp_mod$Cluster <- factor(gp_mod$Cluster)


gp_mod_cluster = new("compareClusterResult", compareClusterResult = gp_mod)
gp_mod_enrich = new("enrichResult", result = gp_mod)


#dotplot
enrichplot::dotplot(gp_mod_cluster)
#ca
p2 <- enrichplot::cnetplot(gp_mod_cluster, categorySize="pvalue")
p2

edox2 <- pairwise_termsim(gp_mod_enrich)
p4 <- emapplot(edox2, cex_category=1.5,layout="kk")

########################################################################################
results_genes = gconvert(x[,1], organism = "hsapiens",
                         target = "ENTREZGENE_ACC", filter_na = FALSE)


ggo <- groupGO(gene     = results_genes,
               OrgDb    = "org.Hs.eg.db",
               ont      = "BP",
               level    = 3,
               readable = TRUE)




head(ggo)


gene.df <- clusterProfiler::bitr(x[,1], fromType = "SYMBOL",
                                 toType = c("ENTREZID"),
                                 OrgDb = org.Hs.eg.db)



gs = gosummaries(x_Hsap_2)



#gostplot(gp_up, interactive = RUE)
#https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html
ego <- enrichGO(gene          = gene.df$ENTREZID,
                OrgDb         = "org.Hs.eg.db",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)

#Grafo aciclico
goplot(ego)
barplot(ego, showCategory=20)
heatplot(ego, showCategory=10)

edox2 <- pairwise_termsim(ego)
p1 <- emapplot(edox2)
p1
#Running function to get graph of a list of features and GO terms

x <- GOCompare::graphGOspecies(df=x_s$result,
                               GOterm_field="term_name",
                               option = 2,
                               numCores=6,
                               saveGraph=FALSE,
                               outdir = NULL)


#################################################

