options(java.parameters = "-Xmx5000m")
library(bartMachine)
library(ICEbox)
library('xlsx')
library(ggplot2)
library('gtools')

### Variable selection for Bayesian Non-Linear Survival Model

### Need bartmachine package
### Same bart settings as Survival model

### --------------------------------------------------------------------------------------

# Assess feature importance by calculating proportion of times a variable is selected for split
# in determining location estimate of survival function

feature_importance = function(surv.output){
  y = as.numeric(colMeans(surv.output$parameters$beta.train))
  bart = bartMachine(X=X[1:500,], y=y, verbose=FALSE)
  return(investigate_var_importance(bart))}

### --------------------------------------------------------------------------------------

# Partial dependence plot to assess influence on covariate

partial_dependence = function(surv.output, feature){
  y = as.numeric(colMeans(surv.output$parameters$beta.train))
  bart = bartMachine(X=X[1:500,], y=y, verbose=FALSE)
  return(pd_plot(bart, j = feature))
}
### --------------------------------------------------------------------------------------

# feature selection procedure by permuting the response to create null distribution

feature_selection = function(surv.output){
  y = as.numeric(colMeans(surv.output$parameters$beta.train))
  bart = bartMachine(X=X[1:500,], y=y, verbose=FALSE)
  return(var_selection_by_permute(bart))
}
feat_imp = feature_selection(ggITE2)
imp = data.frame(Observed = feat_imp$var_true_props_avg, Null = apply(feat_imp$permute_mat,2,function(x)quantile(x,probs=0.90)),
                 Names = names(feat_imp$var_true_props_avg))
imp$Names = factor(imp$Names, levels=imp$Names)
ggplot(data=imp[1:17,]) + geom_bar(aes(x=Names,y=Observed,colour="Observed"),stat='identity',fill='white') + 
  geom_point(aes(x=Names,y=Null,colour='Null'))+
  scale_colour_manual(name="",values=c('Observed'="black", 'Null'="red"))+
  theme(legend.title = element_blank(),legend.position = c(0.9,0.9),
        axis.text=element_text(size=12),
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(colour = "gray90"),
        panel.background = element_rect(fill='white'),
        axis.text.x = element_text(angle = 40, hjust = 1))+labs(x = "", y="Inclusion Proportions")


### --------------------------------------------------------------------------------------

# variable interaction investigation by calculating proportion of times both variables are selected in a given tree

interaction = function(surv.output){
  y = as.numeric(colMeans(surv.output$parameters$beta.train))
  bart = bartMachine(X=X[1:500,], y=y, verbose=FALSE)
  return(interaction_investigator(bart))
}

### --------------------------------------------------------------------------------------


# non-parametric covariate importance test by permuting the response to create null distribution

feature_importance_test = function(surv.output,features, num_perm){
  y = as.numeric(colMeans(surv.output$parameters$beta.train))
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

## Rank test for individual significance

observed_effect = function(){
  X = subset(dat,select=-c(Survival,Status))
  y = as.numeric(colMeans(ggITE2$parameters$beta.train))
  bart = bartMachine(X=X[1:500,], y=y, verbose=FALSE)
  surv.ice = ice(object = bart,X=X,y=y,predictor="Diabetes",verbose=FALSE)
  return(surv.ice$ice_curves[,2]-surv.ice$ice_curves[,1])
}
rank_test = function(n){
  
  null_effects <- matrix(nrow=500,ncol=n)
  for (i in 1:n){
    X$Diabetes <- permute(X$Diabetes)
    y = as.numeric(colMeans(ggITE2$parameters$beta.train))
    bart = bartMachine(X=X[1:500,], y=y, verbose=FALSE)
    surv.ice = ice(object = bart,X=X[1:500,],y=y,predictor="Diabetes",verbose=FALSE)
    ind_ice <- surv.ice$ice_curves[,2]-surv.ice$ice_curves[,1]
    null_effects[,i] <- ind_ice
    print(sprintf("Interation: %d", i))
  }
  return(null_effects)
}
observed = observed_effect()
null = rank_test(500)

pval = numeric(500)
for (i in 1:500){
  pval[i] = mean(abs(observed[i]) <= abs(null[i,])) # p value
}

# average reduction in median survival times for diabetic patients related to diabetes is
exp(mean(observed))-1

mean(pval[which(dat$Diabetes[1:500]==1)]<0.05)
which(pval[which(dat$Diabetes[1:500]==1)]<0.05)
mean(dat$Weight[which(pval[which(dat$Diabetes[1:500]==1)]<0.05)])
mean(dat$Weight[which(pval[which(dat$Diabetes[1:500]==1)]>=0.05)])
mean(dat$Body.Mass.Index[which(pval[which(dat$Diabetes[1:500]==1)]<0.05)])
mean(dat$Body.Mass.Index[which(pval[which(dat$Diabetes[1:500]==1)]>=0.05)])


#test42 = data.frame(null=exp(null42$x)-1)
#test43 = data.frame(null=exp(null43$x)-1)
write.csv(null, 'LVAD_null.csv')
write.csv(truth, 'LVAD_truth.csv')
