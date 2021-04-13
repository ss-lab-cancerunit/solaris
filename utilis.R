require(RColorBrewer)
# # Colour palettes
divergingPal <- c(brewer.pal(11, "Spectral"), "grey")
divergingPal_long <- c("grey",
                       brewer.pal(9, "YlOrRd")[c(2,3,5,7,9)],
                       brewer.pal(9, "YlGnBu")[c(2,3,4,5,6,7,8,9)],
                       brewer.pal(9, "YlGn")[c(8,7,5)],
                       brewer.pal(9, "PuRd")[c(3,5,6,7)],
                       brewer.pal(9, "BuPu")[c(9,8,7,6,5,4,3,2)])

divergingPalFull <- divergingPal_long[c(5,4,3,2,6,7,16:14,13:8,17:30)]
divergingPalMed <- divergingPal_long[c(5,4,3,2,13:8,7,16:14,17:28)]

# # Vectrorized strsplit  
splitvec <- function(vector, split, select, merge = "_"){
  processed <- sapply(vector, function(x){
    separated <- unlist(strsplit(x, split = split))[select]
    if (length(separated) > 1){
      return(paste(separated, collapse = merge))
    } else
      return(separated)
    })
  processed <- unname(processed)
  return(processed)
}