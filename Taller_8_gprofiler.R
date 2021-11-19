#cargar librerias
  require(gprofiler2);library(biomaRt);library(topGO);
require(clusterProfiler);require(GOsummaries)

#cargar ejemplo desde CSV
  url_file = "https://raw.githubusercontent.com/ccsosa/R_Examples/master/GIM.csv"
  x <- read.csv(url_file,header = T)

#alternativa usando gosummaries
  # gl = list(GIM = x[,1]) # Two lists per component
  # gs = gosummaries(gl)
  # plot(gs, fontsize = 8)



###############################################################################################
#gprofiler2
###############################################################################################
#Functional enrichment analysis
  x_s <-  gprofiler2::gost(query = x[,1],
                           organism = "hsapiens", ordered_query = FALSE,
                           multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                           measure_underrepresentation = FALSE, evcodes = FALSE,
                           user_threshold = 0.05, correction_method = "false_discovery_rate",
                           domain_scope = "annotated", custom_bg = NULL,
                           numeric_ns = "", as_short_link = FALSE,
                           sources="GO:BP")
#GO PLOT
  p <- gostplot(x_s, capped = F, interactive = FALSE)
  p


#GO PLOT with the top five of GO
  publish_gostplot(p, highlight_terms = x_s$result$term_id[1:10])

###############################################################################################
#topGO
###############################################################################################



# create GO db for genes to be used using biomaRt - please note that this takes a while
db= useMart('ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl', host="www.ensembl.org")
go_ids= getBM(attributes=c('go_id', 'external_gene_name', 'namespace_1003'), filters='external_gene_name', values=x[,1], mart=db)
gene_2_GO=unstack(go_ids[,c(1,2)])

# remove any candidate genes without GO annotation
keep = x[,1] %in% go_ids[,2]
keep =which(keep==TRUE)
candidate_list=x[,1][keep]
geneList=factor(as.integer(x[,1] %in% candidate_list),levels = c(0,1))
names(geneList)= x[,1]

# Create the class topGOdata
GOdata=new('topGOdata', ontology='BP', allGenes = geneList, annot = annFUN.gene2GO, gene2GO = gene_2_GO)

#classic_fisher_result=runTest(GOdata, algorithm='classic', statistic='fisher')
#https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf

#Run stastical tests
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

#summarize in a table
allRes <- GenTable(GOdata, classicFisher = resultFisher,
                    classicKS = resultKS, elimKS = resultKS.elim,
                    orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)



#diversas formas de obtener p valores
pValue.classic <- score(resultKS)
pValue.elim <- score(resultKS.elim)[names(pValue.classic)]
gstat <- termStat(GOdata, names(pValue.classic))

#add colors for graphics
 colMap <- function(x) {
   .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
   return(.col[match(1:length(x), order(x))])
 }

 gSize <- gstat$Annotated / max(gstat$Annotated) * 4
 gCol <- colMap(gstat$Significant)


#plots
  #plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim",
      #pch = 19, cex = gSize, col = gCol)


showSigOfNodes(GOdata, score(resultKS.elim), firstSigNodes = 5,reverse = T,showEdges = T,useInfo = "np")
