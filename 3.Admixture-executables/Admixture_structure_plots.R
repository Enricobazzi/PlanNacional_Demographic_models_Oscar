######################################
##### Import and Prepare Dataset #####
######################################
# Define the input Directory
wd_input <- "/Users/enricobazzicalupo/admixture_structure/"

# Create a list of the .Q output file names (the ADMIXTURE output files)
Q_sample_files <- list.files(wd_input, pattern="WholeGenome_plink_Humchr_r2_0.3_4lp.*.Q$")


for (i in 1:length(Q_sample_files)){
  # Import the sample file table
  tbl <- read.table(paste0(wd_input,Q_sample_files[[i]]))
  
  barplot(t(as.matrix(tbl)), col=rainbow(i+2),
          xlab="Individual #", ylab="Ancestry", main = paste0("K = ",i + 2), border=NA)
}

###

## Plot K = n for all 10 repetitions ##
# I do this to check that all repetitions give same result for same K. If not it's a case of a K with conflicting signals.

wd_input <- "/Users/enricobazzicalupo/admixture_structure"
rep_dirs <- dir(path = wd_input, pattern = "rep_", full.names = T)
par(mfrow=c(5,2))

for (k in 3:8){
 for (i in 1:length(rep_dirs)){
   tbl <- read.table(paste0(rep_dirs[[i]],"/WholeGenome_plink_Humchr_r2_0.3_4lp.",k,".Q"))
   barplot(t(as.matrix(tbl)), col=rainbow(k),
          ylab="Ancestry", main = paste0("K = ",k), border=NA)
 }
}

## Plot all K = n for one repetition ##
# I do this because I observed no differences in any K value for all the repetitions. rep_dirs[[n]] determines the repetition

wd_input <- "/Users/enricobazzicalupo/admixture_structure"
rep_dirs <- dir(path = wd_input, pattern = "rep_", full.names = T)
par(mfrow=c(3,2))

for (k in 3:8){
    tbl <- read.table(paste0(rep_dirs[[1]],"/WholeGenome_plink_Humchr_r2_0.3_4lp.",k,".Q"))
    barplot(t(as.matrix(tbl)), col=rainbow(k),
            ylab="Ancestry", main = paste0("K = ",k), border=NA)
}
