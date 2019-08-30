adaptive_jump = function(current_jump, jump.mat,wei){
  # take 50 last jump probabilities
  jump_prob <- jump.mat[(wei-48):wei]
  new_jump  <- ifelse(mean(jump_prob)<0.40, max(current_jump - min(0.05,1/(wei %/% 50)),0.01),
                      current_jump + min(0.05,1/(wei %/% 50)))
  return(new_jump)
}
