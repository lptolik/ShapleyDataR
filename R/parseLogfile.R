plot.log <- function(fp.log){
  library(ggplot2)
  con1 <- file(fp.log, "r")
  proc.log1 <- readLines(con1)
  save.starts <- proc.log1[which(stringi::stri_detect_fixed(proc.log1, pattern = "tol="))]
  save.ends <- proc.log1[which((stringi::stri_detect_fixed(proc.log1, pattern = "Save is complete")))]
  save.starts.dt <- lapply(stri_split_fixed(save.starts, pattern = " ", n=3, tokens_only = TRUE), function(x) strptime(x[[3]], format="%X"))
  save.ends.dt <- lapply(stri_split_fixed(save.ends, pattern = " ", n=3, tokens_only = TRUE), function(x) strptime(x[[3]], format="%X"))
  saving.times <- sapply(1:length(save.ends.dt), function(x) difftime(save.ends.dt[[x]], save.starts.dt[[x]], units = "secs"))
  st.df <- as.data.frame(saving.times)
  ggplot(st.df, aes(x = 1:length(saving.times), y = saving.times)) + geom_line() + geom_smooth(method = "loess")
  return(st.df)

}