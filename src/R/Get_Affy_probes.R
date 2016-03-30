# To install hgu133a.db, run
# ## try http:// if https:// URLs are not supported
# source("https://bioconductor.org/biocLite.R")
# biocLite("hgu133a.db")
library("hgu133plus2.db")  # "hgu133b.db")  # "hgu133a.db")
affy_plat <- hgu133plus2.db
columns(affy_plat)
keytypes(affy_plat)
all_probes <- keys(affy_plat, keytype = 'PROBEID')
plat_df <- select(
    affy_plat,
    keys = all_probes,
    keytype = 'PROBEID',
    columns = c("SYMBOL", "ENTREZID", "GENENAME")
)

setwd('~/code/dj_carkinos/data/')
write.csv(
    plat_df,
    # "Affy_U133A_probe_info.csv",
    # "Affy_U113B_probe_info.csv",
    "Affy_U133plus2_probe_info.csv",
    row.names = FALSE
)
