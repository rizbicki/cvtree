#' CV Tree score
#'
#' @param tree an object of the class "phylo" which the function will evaluate
#' @param chars_test a matrix in which each row name is
#' a tip of "tree" and each column is a hold-out character (typically not used in the construction of tree)
#' @param root the tip of "tree" that should be used as a root
#' @param nsim number of Monte Carlo simulations to be used by simmap. Defaults to 10
#' @param cl Number of cores to be used. Defaults to 1
#' @param normalize Whether edge lengths should be normalized to have the same value. Defaults to TRUE
#'
#' @return A list with two components: mean_score, the estimated predictive score (smaller is better),
#'  and sd_score, its standard error
#' @export
#'
#' @examples
#'
#' set.seed(0)
#' my_tree<- ape::rtree(n=26,tip.label=LETTERS)
#' chars_test <- replicate(2,ape::rTraitDisc(my_tree,rate = 0.5))
#' cvtree(my_tree,chars_test,root="A")
cvtree <- function(tree,chars_test,root,nsim=10,
                    cl=1,normalize=TRUE)
{
  if(normalize)
  {
    # Normalize edge lenghts so that all
    tree$edge.length <- rep(1,length(tree$edge.length))
  }

  cl <- parallel::makeCluster(cl)
  doParallel::registerDoParallel(cl)
  tree <- ape::root(tree,root)
  score <- foreach::foreach(ii=1:ncol(chars_test)) %dopar%
    {
      chars <- names(table(chars_test[,ii]))
      if(length(chars)==1)
        return(list(score_mean=0,score_se=0))
      complete_sequence <- phytools::to.matrix(chars_test[,ii],chars)
      prob_chars <- matrix(NA,nrow(chars_test),length(chars))
      for(jj in 1:nrow(chars_test))
      {
        incomplete_sequence <- complete_sequence
        incomplete_sequence[jj,] <- 1/length(incomplete_sequence[jj,])
        reconstruction<-phytools::make.simmap(tree,incomplete_sequence,
                                              nsim=nsim,model="ER")
        reconstruction_results<-phytools::describe.simmap(reconstruction,plot=FALSE)
        if(nrow(reconstruction_results$tips)==1)
        {
          prob_chars[jj,] <- incomplete_sequence[jj,]
        } else {
          which_row <- which(rownames(reconstruction_results$tips)==rownames(incomplete_sequence)[jj])
          prob_chars[jj,] <- reconstruction_results$tips[which_row,]
        }
      }

      score_mean <- mean(rowMeans((complete_sequence-prob_chars)^2))
      score_se <- sqrt(var(rowMeans((complete_sequence-prob_chars)^2))/nrow(prob_chars))
      return(list(score_mean=score_mean,
                  score_se=score_se))
    }
  score_means <- sapply(score,function(x)x$score_mean)
  score_se <- sapply(score,function(x)x$score_se)
  mean_score <- mean(score_means)
  sd_score <- sqrt(sum(score_se^2))/length(score_se)
  parallel::stopCluster(cl)
  return(list(mean_score=mean_score,
              sd_score=sd_score))
}



