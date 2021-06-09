## code to generate `test_obj`

# locate input and output files
inputdir <- system.file('extdata/SAMPLEDATA_HG18',
                        package = 'segmenter')
coordsdir <- system.file('extdata/COORDS',
                         package = 'chromhmmData')
anchorsdir <- system.file('extdata/ANCHORFILES',
                          package = 'chromhmmData')
chromsizefile <- system.file('extdata/CHROMSIZES',
                             'hg18.txt',
                             package = 'chromhmmData')

# run command
test_obj <- segmenter::learn_model(
    inputdir = inputdir,
    coordsdir = coordsdir,
    anchorsdir = anchorsdir,
    chromsizefile = chromsizefile,
    numstates = 3,
    assembly = 'hg18',
    cells = c('K562', 'GM12878'),
    annotation = 'RefSeq',
    binsize = 200
)

# write object to data
usethis::use_data(test_obj, overwrite = TRUE)
