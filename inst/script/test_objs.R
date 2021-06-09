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

# run command with multiple states
test_objs <- lapply(3:8,
                    function(x) {
                      segmenter::learn_model(inputdir = inputdir,
                                             coordsdir = coordsdir,
                                             anchorsdir = anchorsdir,
                                             chromsizefile = chromsizefile,
                                             numstates = x,
                                             assembly = 'hg18',
                                             cells = c('K562', 'GM12878'),
                                             annotation = 'RefSeq',
                                             binsize = 200)
                    })

# remove the segments slot to reduce the size of the object
test_objs <- lapply(test_objs, function(x) {
  x@segment <- list()
  x
})

# write object to data
usethis::use_data(test_objs, overwrite = TRUE)
