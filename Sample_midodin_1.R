

library(BioStudies)
f= getBio("E-MTAB-4632",path ="data")
countTable <- read.table( "data/counts_clincal_rnaseq_updated_July_16.txt",
                         header=TRUE, as.is=TRUE, row.names=1, sep="\t") 
treatment <- factor(rep(c("after", "before"), 7))
colnames((countTable))
individual <- factor(rep(1:7, each=2))
countTable[, c("P_207.S", "P_207.T")] <- countTable[, c("P_207.T", "P_207.S")]
View(countTable)

######################################################################
#miodin
library(miodin)
mp= MiodinProject(name = " trial",author ="deepak",path=".")
sampleTable =data.frame(SampleName =colnames(countTable),
                        samplingPoint =rep("sp1",14),
                        treatment= treatment,
                        individual =individual)
assayTable= data.frame(SampleName =colnames(countTable),
                       DataFile="data/counts_clincal_rnaseq_updated_July_16.txt",
                       DataColumn =colnames(countTable))
ms = studyDesignCaseControl(
  studyName = "trial",
  factorName = "treatment",
  caseName = "after",
  controlName = "before",
  contrastName= "treatment",
  numCase = 7,
  numControl = 7,
  sampleTable = sampleTable,
  assayTable = assayTable,
  assayTableName = "RNAseq",
  paired = "individual"
)
insert(ms,mp)
#work flow object 
mw = MiodinWorkflow("trial")
mw=mw +
  importProcessedData(
    name = "RNA-seq importer",
    experiment = "sequencing",
    dataType = "rna",
    studyName = "trial",
    assayName = "RNAseq",
    datasetName = "E-MTAB-4632",
    contrastName = "treatment"
    )%>%
  processSequencingData(
    name = "RNA-seq processor",
    contrastName = "treatment",
    filterLowCount = TRUE
  )
mw=insert(mw,mp)
mw = execute(mw)
saveDataFile(mp)
export(mp ,"dataset","E-MTAB-4632")
test <- read.table(
  "data/counts_clincal_rnaseq_updated_July_16.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1
)

dim(test)
head(test[,1:3])
