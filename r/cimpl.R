
args <- commandArgs(TRUE)
if(length(args) < 2) {
  print('Not enough arguments provided!')
  q(save='no')
}

insFile = args[1]
outputPrefix = args[2]
 
library(cimpl)
library(reshape)
library(BSgenome.Mmusculus.UCSC.mm10)

parseScales <- function(scaleCol) {
  parseScale <- function(x) { 
    as.numeric(gsub('e.', 'e+', gsub('X', '', as.character(x))))
  }
  sapply(scaleCol, parseScale)
}

splitIds <- function(cisTable, col.id) {
  splitIds = sapply(as.character(cisTable[,col.id]), function(x) { strsplit(x, '|', fixed=T) })
  names(splitIds) = rownames(cisTable)
  
  tableIds = rep(names(splitIds), as.numeric(sapply(splitIds, length)))
  mergeTable = data.frame(ids=tableIds)
  mergeTable[col.id] = unlist(splitIds, use.names=F)
  
  cisTable.noCIS = cisTable[!(colnames(cisTable) %in% c(col.id))]
  merged = merge(cisTable.noCIS, mergeTable, by.x='row.names', by.y='ids')
  merged[!(colnames(merged) %in% c('Row.names'))]
}

# Load data
data = read.csv(insFile, sep='\t', stringsAsFactors=F)

# Do CIMPL analysis
res = doCimplAnalysis(data, BSgenome=BSgenome.Mmusculus.UCSC.mm10, system='SB',
                      lhc.method="exclude", n_iterations=1000, scales = c(5000, (1:15) * 10000))

save(data, res, file=paste(outputPrefix, '.rdump.Rdata', sep=''))

mart = useMart('ENSEMBL_MART_ENSEMBL', dataset='mmusculus_gene_ensembl', host='www.ensembl.org')
genes = getEnsemblGenes(res, mart=mart)

ciss <- getCISs(res, genes=genes, alpha=0.05, mul.test=TRUE)
cisMatrix <- getCISMatrix(res, ciss)

# Reformat output
cols.ids = c('ins_id')
cols.scales = colnames(cisMatrix)[grepl('^X.+', colnames(cisMatrix))]
cisMatrix.scales = cisMatrix[,c(cols.ids, cols.scales)]

cisTable = melt(cisMatrix.scales, id=cols.ids, variable_name='scale')
cisTable = cisTable[cisTable$value != '',]
cisTable = rename(cisTable, c('value'='cis_id'))
cisTable = splitIds(cisTable, 'cis_id')
cisTable$scale = NULL

colnames(ciss) <- c('chromosome', 'peakLocation', 'peakHeight', 'start', 'end', 'width',
                    'nInsertions', 'pValue', 'scale', 'associatedEnsemblGeneId', 'otherEnsemblGeneId',
                    'associatedExternalGeneId', 'otherExternalGeneId')
ciss$id = rownames(ciss)

# Write output
options(scipen=999)

write.table(ciss, file=paste(outputPrefix, '.cis-list.tsv', sep=''), sep='\t', row.names=F)
write.table(cisTable, file=paste(outputPrefix, '.cis-ins-mapping.tsv', sep=''), sep='\t', row.names=F)
#export.bed(ciss, file=paste(outputPrefix, '.cis-list.bed', sep=''))
