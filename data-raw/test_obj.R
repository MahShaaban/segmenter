## code to prepare `test_obj` dataset goes here

# locate input and output files
inputdir <- system.file('extdata/ChromHMM/SAMPLEDATA_HG18',
                        package = 'segmentr')
outputdir <- system.file('extdata/output',
                         package = 'segmentr')
coordsdir <- system.file('extdata/ChromHMM/COORDS',
                         package = 'segmentr')
anchorsdir <- system.file('extdata/ChromHMM/ANCHORFILES',
                          package = 'segmentr')
chromsizefile <- system.file('extdata/ChromHMM/CHROMSIZES',
                             'hg18.txt',
                             package = 'segmentr')
# run command
test_obj <- segmentr::learn_model(
    inputdir = inputdir,
    outputdir = outputdir,
    coordsdir = coordsdir,
    anchorsdir = anchorsdir,
    chromsizefile = chromsizefile,
    numstates = 3,
    assembly = 'hg18',
    cells = c('K562', 'GM12878'),
    annotation = 'RefSeq',
    binsize = 200,
    read_only = TRUE
)

# write object to data
usethis::use_data(test_obj, overwrite = TRUE)

# run command with multiple states
test_objs <- lapply(3:8,
                    function(x) {
                        segmentr::learn_model(inputdir = inputdir,
                                              coordsdir = coordsdir,
                                              anchorsdir = anchorsdir,
                                              chromsizefile = chromsizefile,
                                              numstates = x,
                                              assembly = 'hg18',
                                              cells = c('K562', 'GM12878'),
                                              annotation = 'RefSeq',
                                              binsize = 200)
                    })

# write object to data
usethis::use_data(test_objs, overwrite = TRUE)
