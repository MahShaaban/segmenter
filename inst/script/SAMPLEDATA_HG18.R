## code to download sample data from 
## source: https://github.com/jernst98/ChromHMM
## to: inst/extdata/SAMPLEDATA_HG18

# create a directory
dir.create('inst/extdata/SAMPLEDATA_HG18')

# create a temporary directory
temp <- tempdir()

# download the github rep as a zip archive
download.file('https://github.com/jernst98/ChromHMM/archive/master.zip',
              destfile = file.path(temp, 'master.zip'))

# extract the sample files to destination
unzip(file.path(temp, 'master.zip'),
      files = c('ChromHMM-master/SAMPLEDATA_HG18/K562_chr11_binary.txt.gz',
                'ChromHMM-master/SAMPLEDATA_HG18/GM12878_chr11_binary.txt.gz'),
      junkpaths = TRUE,
      exdir = 'inst/extdata/SAMPLEDATA_HG18')
