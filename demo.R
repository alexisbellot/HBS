""
"Example performance estimation on the publicly available GBSG2 cancer data "
""

## Required packages

library('survival')
library('pec')
library('ggplot2')


## Load model and utility function for performance computation

source('/utils.R')
source('/HBS.R')
source('/Variable Importance and tests for HBS model.R')

## Load and preprocess data ---------------------------------------------

data(GBSG2)
data = GBSG2
colnames(data)[10]<- 'Status'
data$horTh= as.numeric(data$horTh) -1
data$menostat = as.numeric(data$menostat) -1
data$tgrade1 = as.numeric(data$tgrade=="I")
data$tgrade2 = as.numeric(data$tgrade=="II")
#data$tgrade3 = as.numeric(data$tgrade=="III")
data$time = data$time/7
data = subset(data, select=-c(tgrade))
colnames(data)[1] = "group"
data$group = data$group+1
colnames(data)[8] = "Survival"
Y = subset(data,select=c(Survival,Status))
X = subset(data,select=-c(Survival,Status))

## Run algorithm and plot predictions ----------------------------------------------

hbs = HBS(X.train=X,Y.train=Y,X.test=X,n.iter=2000,ntree=50,burn.in=500,thinning=5,n.group=2)

# Create survival matrix and estimate probabilities
predictions = matrix(NA, nrow(X),100)
for (j in 1:2){
  for (i in which(X$group == j)){
    predictions[i,] = 1 - pgengamma(seq(1,100,length.out = 100),mu=median(hbs$beta.test[,i]),
                                    sigma=exp(median(hbs$logsig.mat[,j])),
                                    Q=median(hbs$loglam.mat[,j])) }}

predictions_frame = data.frame(times = seq(1,100,length.out = 100), patient_1=predictions[9,1:100],
                               patient_2=predictions[10,1:100])

# Plot predictions for the first 2 patients
ggplot(predictions_frame, aes(times))+ xlim(-1,100) + ylim(0.5,1)+
  geom_line(aes(y=patient_1,col="1"),size=1)+
  geom_line(aes(y=patient_2, col="2"),size=1)+
  scale_colour_manual(name="",values=c("1"="orange", "2"="blue"))+
  theme(legend.position="",legend.text=element_text(size=10),
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(colour = "gray90"),
        panel.background = element_rect(fill='white'))+labs(x = "Time (Weeks)", y="Survival probability")


## Run algorithm and compute performance (C index and brier score) -------------------------------------------

perf = perf.est(data,k=3,n.group=2)

## Feature importance example --------------------------------------------------------------

feat_imp = feature_selection(hbs)

imp = data.frame(Observed = feat_imp$var_true_props_avg, Null = apply(feat_imp$permute_mat,2,function(x)quantile(x,probs=0.90)),
                 Names = names(feat_imp$var_true_props_avg))
imp$Names = factor(imp$Names, levels=imp$Names)

ggplot(data=imp[1:8,]) + geom_bar(aes(x=Names,y=Observed,colour="Observed"),stat='identity',fill='white') + 
  geom_point(aes(x=Names,y=Null,colour='Null'))+
  scale_colour_manual(name="",values=c('Observed'="black", 'Null'="red"))+
  theme(legend.title = element_blank(),legend.position = c(0.9,0.9),
        axis.text=element_text(size=12),
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(colour = "gray90"),
        panel.background = element_rect(fill='white'),
        axis.text.x = element_text(angle = 40, hjust = 1))+labs(x = "", y="Inclusion Proportions")
