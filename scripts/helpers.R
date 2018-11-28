library(data.table)

getExtension <- function(file){
  ex <- strsplit(basename(file), split="\\.")[[1]]
  return(ex[2])
}

get_df <- function (f){

  extension = getExtension(f)
  if (extension == "gtf"){
    cols = "1,4,5,7"
  } else if (extension == "bed"){
    cols = "1-3,6"
  } else {
    stop(paste0("Bad extension: ", extension, " in ", f))
  }


  if (Sys.info()["sysname"] == "Darwin"){
    cmd = paste0("gzcat ", f, " | cut -f ", cols)
  } else {
    cmd = paste0("zcat ", f, " | cut -f ", cols)
  }
  print(cmd)
  df = fread(cmd, header=FALSE, col.names=c("Chromosome", "Start", "End", "Strand"), stringsAsFactors=TRUE)
  return(df)
}

file_to_grange <- function(f){
  df <- get_df(f)
  return(GRanges(seqnames = df$Chromosome, ranges = IRanges(start = df$Start, end = df$End), strand = df$Strand))
}
