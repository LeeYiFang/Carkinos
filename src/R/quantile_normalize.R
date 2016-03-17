# Install the affy package by
#   source("http://bioconductor.org/biocLite.R")
#   biocLite("affy")
library(affy)

# To download the raw data as TAR file
#   curl -o GSE32474_RAW.tar "http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE32474&format=file"
# Extract under folder GSE32474_RAW
#   tar -C GSE32474_RAW -xf GSE32474_RAW.tar

raw_affy <- ReadAffy(
    celfile.path = "Dataset/GSE32474_RAW",
    compress = TRUE,
    verbose = TRUE
)

eset <- expresso(
    raw_affy,
    bgcorrect.method="none",
    normalize.method="quantiles",
    pmcorrect.method="pmonly",
    summary.method="medianpolish"
)
# Get all probes
# probe_names <- featureNames(eset)
write.exprs(eset, file="GSE32474.quantile_normalized.tsv")
