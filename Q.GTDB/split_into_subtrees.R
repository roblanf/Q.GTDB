# Split an input tree into subtrees of maximum size N
library(castor)

split_tree <- function(tree, Nmax = 1000, Nmin = 4){
  
  to_split = list(tree)
  to_keep = list() # we'll put trees to keep here
  
  # we don't need to do anything if the input tree is already small enough  
  if(Ntip(tree)<=N){
    return(c(tree))
  }
  
  while(length(to_split)>0){
    
    cat("to split: ", length(to_split), "    to keep: ", length(to_keep), "\n")
    
    # poor man's pop: get the first element and delete it
    t = to_split[[1]]
    to_split[[1]] <- NULL
    
    # split the tree (this is effectively just splitting it at the root node of the tree)
    splits <- split_tree_at_height(t, 0.00000000000000000000000000000001)
    # extract the subtrees from the castor object
    subtrees = sapply(splits$subtrees, function(x) x$tree, simplify = FALSE)
    
    for(subtree in subtrees){
      
      if( Ntip(subtree) > Nmax ){ # split it again...
        to_split[[length(to_split) + 1]] <- subtree
      } else if ( Ntip(subtree) > (Nmin - 1)){ # it's between 4 and N, which is what we want
        to_keep[[length(to_keep) + 1]] <- subtree
      } 
      # otherwise we ditch the tree because it has <=3 tips
    }
  }
  
  class(to_keep) <- 'multiPhylo'
  return(to_keep)
}
