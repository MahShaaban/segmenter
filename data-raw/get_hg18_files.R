# load required libraries
library(rtracklayer)
library(GenomeInfoDb)

# start a session
session <- browserSession()

# assign genome name
genome(session) <- "hg18"

# query the session
query <- ucscTableQuery(session, "refGene")

tableName(query) <- "refGene"

# retrieve the track in a data.frame
tab <- getTable(query)

# write table into tsv
write.table(tab, 'inst/extdata/hg18/genome_table.tsv')

# check the table has the correct columns
names(tab)

# get chromosome lengths file
chrom_length <- getChromInfoFromUCSC('hg18')

# write to file
write.table(chrom_length[, c('chrom', 'size')],
            'inst/extdata/hg18/chromsize.tsv',
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE)
