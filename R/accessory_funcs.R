non_vect_code = function(files,verbose = FALSE,col_idx) {
  
  for (i in 1:length(files)) {
  
    if (i == 1) {
      
      data <- if(exists(data)) cbind(data,data_temp) else data_temp
      
    }
    
    
    # optimize the input window
    y <- optimal rows
    x <- optimal cols
    
    for (i in 1:(ncols(data)/y)) {  
      
      data_temp <- fread(input = files[i],
                       select = (x*(i-1)+1):min(x*i,cols), # min() catches % remainder
                       header = TRUE,
                       sep = "\t",
                       stringsAsFactors = FALSE,
                       #colClasses = rep("numeric", 970),
                       col.names = col_idx,
                       data.table = TRUE)
    
      data <- data_temp %>% as.matrix() %>% dropNA() %>% cbind(data, .)
    
    }
  
  }
}