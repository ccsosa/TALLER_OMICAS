#cargar librerias
library(gprofiler2);library(biomaRt);library(topGO);
library(clusterProfiler);library(GOsummaries)

#cargar ejemplo desde CSV
  url_file = "https://raw.githubusercontent.com/ccsosa/R_Examples/master/GIM.csv"
  x <- read.csv(url_file,header = T)


###############################################################################################
#gprofiler2
###############################################################################################
#Analisis de enriquecimiento funcional
  x_s <-  gprofiler2::gost(query = x[,1],
                           organism = "hsapiens", ordered_query = FALSE,
                           multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                           measure_underrepresentation = FALSE, evcodes = FALSE,
                           user_threshold = 0.05, correction_method = "false_discovery_rate",
                           domain_scope = "annotated", custom_bg = NULL,
                           numeric_ns = "", as_short_link = FALSE,
                           sources="GO:BP")
#GO PLOT
  p <- gprofiler2::gostplot(x_s, capped = F, interactive = FALSE)
  p


#GO PLOT con el top diez
  gprofiler2::publish_gostplot(p, highlight_terms = x_s$result$term_id[1:10])

###############################################################################################
#topGO
###############################################################################################


#Obtener los GO desde ENSEMBL
db= biomaRt::useMart('ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl', host="www.ensembl.org")
go_ids= biomaRt::getBM(attributes=c('go_id', 'external_gene_name', 'namespace_1003'),
              filters='external_gene_name', 
              values=x[,1], 
              mart=db)

#cargar en caso de que no funcione desde biomaRt
#write.csv(go_ids,"D:/REPO_GITHUB/TALLER_OMICAS/go_ids.csv",na = "",row.names = F,quote=F)
go_ids <- read.csv("https://raw.githubusercontent.com/ccsosa/TALLER_OMICAS/master/go_ids.csv",header = T)

gene_2_GO=unstack(go_ids[,c(1,2)])


#remover genes sin anotacion
keep = x[,1] %in% go_ids[,2]
keep =which(keep==TRUE)
candidate_list=x[,1][keep]
geneList=factor(as.integer(x[,1] %in% candidate_list),levels = c(0,1))
names(geneList)= x[,1]

#crear un objeto topGOdata
GOdata=new('topGOdata', ontology='BP', allGenes = geneList, annot = topGO::annFUN.gene2GO,
           gene2GO = gene_2_GO)

#classic_fisher_result=runTest(GOdata, algorithm='classic', statistic='fisher')
#https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf

#Run stastical tests
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

#summarize in a table
allRes <- topGO::GenTable(GOdata, classicFisher = resultFisher,
                    classicKS = resultKS, elimKS = resultKS.elim,
                    orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)



#diversas formas de obtener p valores
pValue.classic <- score(resultKS)
pValue.elim <- score(resultKS.elim)[names(pValue.classic)]
gstat <- topGO::termStat(GOdata, names(pValue.classic))

par(cex = 1)
showSigOfNodes(GOdata, score(resultKS.elim), firstSigNodes = 3)



###############################################################################################
#GOsummaries
###############################################################################################
#alternativa usando gosummaries
gl = list(GIM = x[,1]) # Two lists per component
gs = GOsummaries::gosummaries(gl)
plot(gs, fontsize = 8)

