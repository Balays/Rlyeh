#### !!! IMPORTANT !!! RENAMING BAMFILES !!!
###
rename.bamfiles <- F
if(rename.bamfiles) {
  from <- list.files(bamdir, pattern, recursive = T, full.names = T)
  to   <- gsub('PRV_MdBio_', '', from)
  to   <- gsub('_virus_only_filtered', '', to)
  file.rename(from, to)
  
  bamfiles <- grep('.bai', 
                   list.files(bamdir, pattern, recursive = T, full.names = T), 
                   invert = T, value = T)
}
#### !!!
##
