########################################################################################
#cargar librerias
require(gprofiler2);library(biomaRt);library(topGO);
require(clusterProfiler);require(GOsummaries)

###############################################################################################
#gprofiler2
###############################################################################################

url_file = "https://raw.githubusercontent.com/ccsosa/TALLER_OMICAS/master/Hallmarks_of_Cancer_AT.csv"
x_group <- read.csv(url_file)
#leyendo archivo
x_group[,1] <- NULL

#definiendo categorias
CH <- c("AID","AIM","DCE","ERI","EGS","GIM","IA","RCD","SPS","TPI")

#preparacion de los datos
x_Hsap <- lapply(seq_len(length(CH)), function(i){
  x_unique <- unique(na.omit(x_group[,i]))
  x_unique <- x_unique[which(x_unique!="")]
  x_unique <- as.list(x_unique)
  return(x_unique)
})

names(x_Hsap) <- CH
#Correr enriquecimiento funcional para todas las listas de genes
x_s_group <-  gprofiler2::gost(query = x_Hsap,
                               organism = "hsapiens", ordered_query = FALSE,
                               multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                               measure_underrepresentation = FALSE, evcodes = T,
                               user_threshold = 0.05, correction_method = "fdr",
                               domain_scope = "annotated", custom_bg = NULL,
                               numeric_ns = "", sources = "GO:BP", as_short_link = FALSE)

res_group <- x_s_group$result

#Obtener numero de GO enriquecidos por lista de genes
tapply(res_group$query,res_group$query,length)

#Obtener top tres GO por categoria
x_res <- list()
for(i in 1:length(CH)){
  x_res[[i]] <- res_group[which(res_group$query==CH[[i]]),]
  x_res[[i]] <- x_res[[i]][c(1:3),]
}
x_res <- do.call(rbind, x_res)


#obtener la tabla de frecuencias para graficar
x_res$query <- factor(x_res$query,levels = CH)
other_table <- table( x_res$query,x_res$term_name)

par(mar = c(13, 4, 11, 1) + 0.2) #add room for the rotated labels
#par(mar=c(1,1,1,1))

bar <- barplot(other_table,
        main = "",
        xlab = "", ylab = "Frecuencia",
        col = rainbow(10),
        xaxt="n",
        axes=T,
       # legend.text = rownames(other_table),
        beside = F,
        las=2,horiz = F,
        space=0,cex.names = 0.8)

labs <- paste(colnames(other_table))
text(cex=1, x=bar-1, y=-.52,labs, xpd=TRUE, srt=45)
#axis(2, at = 0:5, labels = 0:5)
legend("top", rownames(other_table), fill = rainbow(10), bty = "n",horiz = T,inset = c(0,-0.5),xpd = T,
       cex = 0.5)


###############################################################################################
#clusterProfiler
###############################################################################################

#Organizar la lista de genes para usarse en clusterProfiler
x_Hsap_2 <- list()
for(i in 1:length(x_Hsap)){
  x_Hsap_2[[i]] <- clusterProfiler::bitr(as.character(unlist(x_Hsap[[i]])),
                                         fromType = "SYMBOL",
                                         toType = c("ENTREZID"),
                                         OrgDb = "org.Hs.eg.db")[,2]
}
names(x_Hsap_2) <-CH


x_compare <- clusterProfiler::compareCluster(geneClusters=x_Hsap_2,enrichGO, OrgDb = org.Hs.eg.db)

clust_results <- x_compare@compareClusterResult
tapply(clust_results$Cluster,clust_results$Cluster,length)


#Graficar en un dotplot
enrichplot::dotplot(x_compare)

#Graficar una red 
enrichplot::cnetplot(x_compare)

#
###############################################################################################
#ver los resultados de gprofiler en clusterprofiler
###############################################################################################
#Modificar archivo de resultados para correr enrichplot
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

# Convertir de gprofiler a clusterprofiler
gp_mod_cluster = new("compareClusterResult", compareClusterResult = gp_mod)
gp_mod_enrich = new("enrichResult", result = gp_mod)


#dotplot
enrichplot::dotplot(gp_mod_cluster)
#red de interacciones
p2 <- enrichplot::cnetplot(gp_mod_cluster)
p2


###############################################################################################
#GOsummaries
###############################################################################################
#Ajustar el formato a GOSummaries
x_Hsap3 <- lapply(seq_len(length(CH)), function(i){
  x_unique <- as.character(x_Hsap[[i]])
  return(x_unique)
})

names(x_Hsap3) <- CH

#Correr el analisis y graficar
gs1 = gosummaries(x_Hsap3)
plot(gs1[1:5])
plot(gs1[6:10])


###############################################################################################
#GOCompare
###############################################################################################
#Cargar datos de ejemplo (cuatro cancer hallmarks)
data(H_sapiens_compress)
data(A_thaliana_compress)

#Definir la columna que tiene la info de los GO enriquecidos
GOterm_field <- "Functional_Category"

#Nombrar las especies
species1 <- "H. sapiens"
species2 <- "A. thaliana"

x <- compareGOspecies(df1=H_sapiens_compress,
                      df2=A_thaliana_compress,
                      GOterm_field=GOterm_field,
                      species1=species1,
                      species2=species2)

#graficar el PCoA
x$graphics
  
#cluster con las distancias
plot(hclust(x$distance,"ward.D"))

#Extraer pesos para GO enriquecidos en una sola especie 
x_graph <- graphGOspecies(df=H_sapiens_compress,
                    GOterm_field=GOterm_field,
                    option = "GO",
                    numCores=2,
                    saveGraph=FALSE,
                    outdir = NULL)


#Extraer pesos para GO enriquecidos entre dos especies y categorias
x_graph_two <- graph_two_GOspecies(x=x,
                               species1=species1,
                               species2=species2,
                               GOterm_field=GOterm_field,
                               numCores=2,
                               saveGraph = FALSE,
                               option= "GO",
                               outdir = NULL)


