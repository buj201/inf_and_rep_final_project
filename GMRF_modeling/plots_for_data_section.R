setwd("~/Desktop/Other/MS_Courses/Inf_and_rep/project/nycSchoolPredictions/GMRF_modeling")

source('get_data.R')

library(reshape)
library(ggplot2)

melt_train = melt(train[,!(colnames(train) %in% c('CT2010','School.Year'))])

png('figs/grade_distribution.png',width=800,height=700, res=144)
ggplot(data=melt_train, aes(x=value, color=variable)) +
  geom_density() +
  ggtitle('Distributions of Counts by Grade') +
  xlab('Number of Students') +
  ylab('Density')
dev.off()

melt_train_2 = melt(train[,!(colnames(train) %in% c('CT2010'))],id.vars=c('School.Year'))
melt_train_2$School.Year = as.factor(melt_train_2$School.Year)

png('figs/year_distribution.png',width=800,height=700, res=144)
ggplot(data=melt_train_2, aes(x=value, color=School.Year)) +
  geom_density() +
  ggtitle('Distributions of Counts by School Year') +
  xlab('Number of Students') +
  ylab('Density')
dev.off()

library(matrixStats)

#Row/col/diagonal SDs

train_split = split(train, train$CT2010)
for (i in 1:length(train_split)){
  train_split[[i]] = train_split[[i]][,!(names(train_split[[i]]) %in% c('CT2010','School.Year'))]
  rownames(train_split[[i]]) = c('20012002','20022003','20032004','20042005','20052006','20062007','20072008','20082009')
  train_split[[i]] = t(train_split[[i]])
  train_split[[i]] = train_split[[i]][c(6,1,2,3,4,5),]
  train_split[[i]] = matrix(train_split[[i]],nrow(train_split[[i]]), ncol(train_split[[i]]))
}

grade_sds = list()
year_sds = list()
cohort_sds = list()
for (i in 1:length(train_split)){
  row_sds = rowSds(train_split[[i]])/rowMeans(train_split[[i]])
  col_sds = colSds(train_split[[i]])/colMeans(train_split[[i]])
  diag_sds = rep(0,3)
  for (j in 1:3){
    temp = rep(0,6)
    for (counter in 0:5){
      temp[counter+1] = train_split[[i]][1+counter, j+counter]
    }
    diag_sds[j] = sd(temp)/mean(temp)
  }
  grade_sds[[i]] = row_sds
  year_sds[[i]] = col_sds
  cohort_sds[[i]] = diag_sds
}

grade_sds = unlist(grade_sds)
year_sds = unlist(year_sds)
cohort_sds = unlist(cohort_sds)

var_plot_df = cbind(grade_sds,year_sds,cohort_sds)
colnames(var_plot_df) = c('Grade','Year','Cohort')
var_plot_df = melt(var_plot_df)
colnames(var_plot_df) = c('__','Orientation','value')


png('figs/directions_of_variance.png',width=800,height=700, res=144)
ggplot(data=var_plot_df, aes(x=value, color=Orientation)) +
  geom_density() +
  ggtitle('Relative Standard Deviation for Different Orientations') +
  xlab('Relative SD (SD/Mean)') +
  ylab('Density')
dev.off()

library(gplots)
png('figs/heat_map.png',width=800,height=700, res=144)
heatmap_data = train_split[[4]]
rownames(heatmap_data) = c('K','1','2','3','4','5')
colnames(heatmap_data) = c('01-02','02-03','03-04','04-05','05-06','06-07','07-08','08-09')
heatmap.2(heatmap_data,Rowv=FALSE, Colv=FALSE,dendrogram='none',scale="row",
          density.info="none", trace="none",
          margins=c(6,0), # ("margin.Y", "margin.X")
          key.par=list(mar=c(3.5,0,3,0)),
          heatmap.par(mar=c(4,0,1,0)),
          # lmat -- added 2 lattice sections (5 and 6) for padding
          lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(2.5, 5), lwid=c(1, 10, 1))
dev.off()