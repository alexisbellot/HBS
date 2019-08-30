options(java.parameters = "-Xmx5000m")
library(bartMachine)
library(ICEbox)
library('xlsx')
library(ggplot2)
library('gtools')

### Variable selection for Hierarchical Bayesian Survival Model


### --------------------------------------------------------------------------------------

# Assess feature importance by calculating proportion of times a variable is selected for split
# in determining location estimate of survival function

feature_importance = function(surv.output){
  y = as.numeric(colMeans(surv.output$beta.train))
  bart = bartMachine(X=X[1:500,], y=y, verbose=FALSE)
  return(investigate_var_importance(bart))}

### --------------------------------------------------------------------------------------

# Partial dependence plot to assess influence on covariate

partial_dependence = function(surv.output, feature){
  y = as.numeric(colMeans(surv.output$beta.train))
  bart = bartMachine(X=X[1:500,], y=y, verbose=FALSE)
  return(pd_plot(bart, j = feature))
}
### --------------------------------------------------------------------------------------

# feature selection procedure by permuting the response to create null distribution

feature_selection = function(surv.output){
  y = as.numeric(colMeans(surv.output$beta.train))
  bart = bartMachine(X=X, y=y, verbose=FALSE)
  return(var_selection_by_permute(bart))
}


### --------------------------------------------------------------------------------------

# variable interaction investigation by calculating proportion of times both variables are selected in a given tree

interaction = function(surv.output){
  y = as.numeric(colMeans(surv.output$beta.train))
  bart = bartMachine(X=X[1:500,], y=y, verbose=FALSE)
  return(interaction_investigator(bart))
}

### --------------------------------------------------------------------------------------


# non-parametric covariate importance test by permuting the response to create null distribution

feature_importance_test = function(surv.output,features, num_perm){
  y = as.numeric(colMeans(surv.output$beta.train))
  bart = bartMachine(X=X[1:500,], y=y, verbose=FALSE)
  return(cov_importance_test(bart,covariates = features,num_permutation_samples = num_perm))
}

# individual conditional expectation plot to assess influence of feature on single observations

ice_plot = function(surv.output, feature,centered){
  y = c(as.numeric(colMeans(surv.output$beta.train)), as.numeric(colMeans(surv.output$beta.test)))
  bart = bartMachine(X=X, y=y, verbose=FALSE)
  surv.ice = ice(object = bart,X=X,y=y,predictor=feature)
  #return(surv.ice)
  return(plot(surv.ice, x_quantile = TRUE, plot_pdp = TRUE,centered = centered))
}

## Rank test for individual significance of a variable

observed_effect = function(model){
  X = subset(dat,select=-c(Survival,Status))
  y = as.numeric(colMeans(model$beta.train))
  bart = bartMachine(X=X[1:500,], y=y, verbose=FALSE)
  surv.ice = ice(object = bart,X=X,y=y,predictor="Variable_Name",verbose=FALSE)
  return(surv.ice$ice_curves[,2]-surv.ice$ice_curves[,1])
}
null_effect = function(n){
  ## Need to specify Variable name
  null_effects <- matrix(model,nrow=500,ncol=n)
  for (i in 1:n){
    X$Variable_Name <- permute(X$Variable_Name)
    y = as.numeric(colMeans(model$beta.train))
    bart = bartMachine(X=X, y=y, verbose=FALSE)
    surv.ice = ice(object = bart,X=X,y=y,predictor="Variable_Name",verbose=FALSE)
    ind_ice <- surv.ice$ice_curves[,2]-surv.ice$ice_curves[,1]
    null_effects[,i] <- ind_ice
    print(sprintf("Interation: %d", i))
  }
  return(null_effects)
}

test = function(model,n){
  observed = observed_effect()
  null = rank_test(n)

  pval = numeric(n)
  for (i in 1:n){
    pval[i] = mean(abs(observed[i]) <= abs(null[i,])) # p value
  }
  return(p_val)
}

