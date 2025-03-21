library(data.table)
library(stringr)
library(msa)
library(Biostrings)

#' Title
#'
#' @param fasta_path 
#'
#' @return
#' @export
#'
#' @examples
fasta_header = function(fasta_path){
  
  dt = read_lines(fasta_path)
  dt = dt[grep('^>.*', dt)]
  dt = str_replace_all(dt,"(ave_reads|posterior)=",";")
  dt = gsub("^>haplotype[0-9]\\s+|\\s*\\|\\s*", "", dt)
  dt = gsub("\\s+", "", dt) %>% as.data.table()
  dt = dt[, c("hap", "posterior","ave_reads") := tstrsplit(., ";", fixed=TRUE)]
  dt = dt[,.:=NULL]
  dt[,"barcode"] = str_extract(fasta_path,"barcode[0-9]*")
  dt[,"window"] = str_extract_all(fasta_path,"-[0-9]*-[0-9]*")
  dt[,"window"] = str_remove(dt[,window],"^-")
  dt = dt[, ave_reads:=as.numeric(ave_reads)]
  dt = dt[, posterior:=as.numeric(posterior)]
  
  
  return(dt)
  
}

#' Title
#'
#' @param fasta_path 
#'
#' @return
#' @export
#'
#' @examples
#' 
fasta_seq = function(fasta_path){
  
  
  f = readDNAStringSet(fasta_path)
  
  #metadata for filtering later needed
  hap = str_extract(names(f),"hap_[0-9]*")
  ave_reads = str_extract(names(f),"ave_reads=[0-9]*")%>% 
    gsub("ave_reads=", "", .)
  posterior = str_extract(names(f),"posterior=[0-9]\\.[0-9]*")%>% 
    gsub("posterior=", "", .)
  barcode= str_extract(fasta_path,"barcode[0-9]*")
  window=str_extract_all(fasta_path,"-[0-9]*-[0-9]*") %>% 
    gsub("^-", "", .)
  
  meta = DataFrame(
    hap = hap,
    ave_reads = ave_reads,
    posterior = posterior,
    barcode = barcode,
    window = window
  )
  
  meta$posterior = as.numeric(meta$posterior)
  meta$ave_reads = as.numeric(meta$ave_reads)
  
  mcols(f) = meta
  names(f) = paste0(meta$barcode,"_",meta$hap)
  
  return(f)
  
}

#' Title
#'
#' @param seqs sequences red in with fasta_seq()
#' @param window string format "start-end" 
#'
#' @return tree ready to be plotted
#' @export
#'
#' @examples
#' 
filter_seq = function(seqs,
                      window,
                      posterior=0.99,
                      ave_reads=7,
                      aa_conversion = FALSE){
  
  seq_window = seqs[(mcols(seqs)$window == window) &
                      (mcols(seqs)$posterior >= posterior) &
                      (mcols(seqs)$ave_reads >= ave_reads)]
  
  if (aa_conversion == TRUE){
    
    seq_window = Biostrings::translate(seq_window)
    
  }else{
    
    seq_window = seq_window
    
  }
  
  return(seq_window)
  
}


#' Title
#'
#' @param gisaid_fasta_file 
#' @param gisaid_metadata 
#' @param collection_start 
#' @param collection_end 
#' @param target_location 
#'
#' @return
#' @export
#'
#' @examples
#' 
making_M_consensus <- function(gisaid_fasta_file,
                               gisaid_metadata,
                               collection_start = '2023-12-01',
                               collection_end = '2024-03-31',
                               target_location = 'Switzerland'){
  
  
  meta_dt = suppressWarnings({fread(paste("zstdcat",gisaid_metadata))})
  meta_dt[["Collection_Date"]] = strftime(as.Date(meta_dt$Collection_Date,
                                                  format="%Y-%m-%d"))
  
  meta_filt_dt = meta_dt[Collection_Date >= as.Date(collection_start) &
                           Collection_Date <= as.Date(collection_end) &
                           grepl(target_location, Location)]
  
  meta_filt_dt = meta_filt_dt[, c("Collection_Date", "Isolate_Id","Clade")]
  
  
  s = suppressWarnings({readDNAStringSet(gisaid_fasta_file)})
  seq_name = names(s)
  sequence = paste(s)
  
  seq_dt <- as.data.table(data.frame(seq_name,sequence))
  
  seq_dt[, c("EPI_ID") := str_extract(seq_name,"EPI_ISL_[0-9]+") ]
  #unique segment names found, not sure if they are all correct, but needed ones are there
  seq_dt[, c("segment") := str_extract(seq_name,"NA|PB2|PB1|MP|NP|NS|PA|HA|HE|P3") ]
  
  #only MP needed
  dt_typing = seq_dt[segment == 'MP', ]
  
  dt_typing = merge(dt_typing,meta_filt_dt, by.x="EPI_ID",by.y="Isolate_Id")
  dt_typing[["week"]] = strftime(as.Date(dt_typing$Collection_Date, format="%Y-%m-%d"),
                                 format = "%V")
  
  M_seq = DNAStringSet(dt_typing$sequence)
  
  M_msa =  msa(M_seq)
  
  aligned_DNAStringSet <- as(M_msa, "DNAStringSet")
  
  consensus_mat <- consensusMatrix(aligned_DNAStringSet)
  
  consensus_seq <- apply(consensus_mat, 2, function(x) {
    names(x)[which.max(x)]
  })
  consensus_sequence <- DNAString(paste(consensus_seq, collapse = ""))
  
  return(consensus_sequence)
  
}



#' Title
#'
#' @param AAStringSet 
#'
#' @return
#' @export
#'
#' @examples
aa_consensus <- function(AAStringSet){
  
  msa =  msa(AAStringSet)
  
  aligned_AAStringSet <- as(msa, "AAStringSet")
  
  consensus_mat <- consensusMatrix(aligned_AAStringSet)
  
  consensus_seq <- apply(consensus_mat, 2, function(x) {
    names(x)[which.max(x)]
  })
  
  consensus_sequence <- AAString(paste(consensus_seq, collapse = ""))
  
  return(consensus_sequence)
  
}




