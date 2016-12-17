setwd("~/Desktop/Other/MS_Courses/Inf_and_rep/project/nycSchoolPredictions/GMRF_modeling")

exp_corr_error = read.csv('./stan_models/exponential_change_corr_RMSE_error.csv')
exp_no_corr_error = read.csv('./stan_models/exponential_change_no_corr_RMSE_error.csv')
lin_corr_error = read.csv('./stan_models/linear_change_corr_RMSE_error.csv')
lin_no_corr_error = read.csv('./stan_models/linear_change_no_corr_RMSE_error.csv')
no_change_error = read.csv('./stan_models/no_change_RMSE_error.csv')

library(gplots)
library(RColorBrewer)
cols_red_yellow_green = colorRampPalette(c("green", "yellow", "red"))(n = 299)

make_plot = function(input_df, title, save_file_path) {
  png(save_file_path,width=800,height=700, res=144)
  heatmap_data = round(data.matrix(input_df[1:5,2:6]),2)
  rownames(heatmap_data) = c('Given K','Given 1','Given 2','Given 3','Given 4')
  colnames(heatmap_data) = c('Predict 1','Predict 2','Predict 3','Predict 4','Predict 5')
  heatmap.2(heatmap_data,Rowv=FALSE, Colv=FALSE,dendrogram='none',key=FALSE,
            density.info="none", trace="none", cellnote=heatmap_data,
            notecex=1.0, notecol="black", na.color='grey',
            margins=c(2,0), # ("margin.Y", "margin.X")
            heatmap.par(mar=c(0,0,0,0)),
            col=cols_red_yellow_green,
            # lmat -- added 2 lattice sections (5 and 6) for padding
            lmat=rbind(c(7,8,9),c(5, 1, 2), c(6, 4, 3)), lhei=c(0.5,5,1.5), lwid=c(1, 10, 4))
  dev.off()
}

make_plot(exp_corr_error, 'Exponential Growth with GMRF Prior','./figs/exp_corr_error.png')
make_plot(exp_no_corr_error, 'Exponential Growth without GMRF Prior', './figs/exp_no_corr_error.png')
make_plot(lin_corr_error, 'Linear Growth with GMRF Prior','./figs/lin_corr_error.png')
make_plot(lin_no_corr_error, 'Linear Growth without GMRF Prior','./figs/lin_no_corr_error.png')
make_plot(no_change_error, 'Baseline Model- 0 Growth', './figs/no_change_error.png')