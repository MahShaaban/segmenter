## code to generate `inst/extdata/cell_mark_table.tsv`

df <- cbind(
  paste0('cell', rep(c(1,2), each = 2)),
  paste0('mark', rep(c(1,2), times = 2)),
  rep('randomBam.bam', 4)
)

write.table(df,
            'inst/extdata/cell_mark_table.tsv',
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE)
