
  timer$start("git")
  
  gitMessage <- paste0(runReason," AOI for this run: ",vec_to_sentence(aoiName))

  if(commitToGit) {
    
    gitadd()
    
    gitcommit(msg = gitMessage)
    
    gitpush()
    
  }

  timer$stop("git", comment = paste0("commit to git = ",commitToGit,". ",gitMessage))
  