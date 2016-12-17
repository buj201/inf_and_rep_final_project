source('get_data.R')

#Now merge wide2 dataset with the SpatialPolygosDataFrame

CTs_for_D20@data = merge(CTs_for_D20@data, wide2, by='CT2010')

#Make maps showing cohort

library(sp)

non_counts = c('CT2010','CTLabel','BoroCode','BoroName','BoroCT2010',
               'CDEligibil','NTACode','NTAName','PUMA','Shape_Leng','Shape_Area')
breakpoints = seq(0, max(CTs_for_D20@data[,!(names(CTs_for_D20) %in% non_counts)])+5, length.out = 11)

png('figs/cohort_map.png',width=800,height=700, res=144)
spplot(CTs_for_D20,
       c("Count.of.Students.K.20012002",
         "Count.of.Students.1.20022003",
         "Count.of.Students.2.20032004",
         "Count.of.Students.3.20042005",
         "Count.of.Students.4.20052006",
         "Count.of.Students.5.20062007"),
       names.attr = c("K, 01/02","1, 02/03","2, 03/04", "3, 04/05", '4, 05/06','5, 06/07'),
       col = "transparent",
       colorkey=list(space="bottom"),
       scales = list(draw = FALSE),
       at = breakpoints,
       col.regions = terrain.colors(n = length(breakpoints)-1),
       main = "Tracking a single cohort, 01/02-06/07",
       sub = 'Number of students')
dev.off()

png('figs/year_map.png',width=800,height=700, res=144)
spplot(CTs_for_D20,
       c("Count.of.Students.K.20012002",
         "Count.of.Students.K.20022003",
         "Count.of.Students.K.20032004",
         "Count.of.Students.K.20042005",
         "Count.of.Students.K.20052006",
         "Count.of.Students.K.20062007"),
       names.attr = c("K, 01/02","K, 02/03","K, 03/04", "K, 04/05", 'K, 05/06','K, 06/07'),
       col = "transparent",
       colorkey=list(space="bottom"),
       scales = list(draw = FALSE),
       at = breakpoints,
       col.regions = terrain.colors(n = length(breakpoints)-1),
       main = "Across Years for a Single Grade (K)",
       sub = 'Number of students')
dev.off()

png('figs/grade_map.png',width=800,height=700, res=144)
spplot(CTs_for_D20,
       c("Count.of.Students.K.20012002",
         "Count.of.Students.1.20012002",
         "Count.of.Students.2.20012002",
         "Count.of.Students.3.20012002",
         "Count.of.Students.4.20012002",
         "Count.of.Students.5.20012002"),
       names.attr = c("K, 01/02","1, 01/02","2, 01/02", "3, 01/02", '4, 01/02','5, 01/02'),
       col = "transparent",
       colorkey=list(space="bottom"),
       scales = list(draw = FALSE),
       at = breakpoints,
       col.regions = terrain.colors(n = length(breakpoints)-1),
       main = "Across Grades for a Single Year (01-02)",
       sub = 'Number of students')
dev.off()
