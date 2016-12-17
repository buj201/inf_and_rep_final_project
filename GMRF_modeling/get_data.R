setwd("~/Desktop/Other/MS_Courses/Inf_and_rep/project/nycSchoolPredictions/GMRF_modeling")

library(rgdal)
library(rgeos)

#Read in shapefiles
CTs = readOGR('./data/nyct2010_16d',"nyct2010")
Districts = readOGR('./data/nysd_16d',"nysd")

#Select District 20 and Brooklyn for efficiency
D20 = Districts[Districts$SchoolDist == 20, ]
Brooklyn = CTs[CTs$BoroName == 'Brooklyn', ]

#Find overlapping census tracts in D20
overlaps = gIntersects(D20, Brooklyn, byid=TRUE)
CTs_for_D20 = Brooklyn[which(overlaps), ]

#Make census tract adjacency matrix
adjacency_matrix = gTouches(CTs_for_D20, byid=TRUE)
rownames(adjacency_matrix) = CTs_for_D20$CT2010
colnames(adjacency_matrix) = CTs_for_D20$CT2010
adjacency_matrix = adjacency_matrix[order(rownames(adjacency_matrix)),order(colnames(adjacency_matrix))]
adj = matrix(as.numeric(adjacency_matrix), nrow=dim(adjacency_matrix)[1], ncol=dim(adjacency_matrix)[2])

#Now parse census tracts to make short codes
short_codes = list()
index = 1
for (ct in CTs_for_D20$CT2010){
  short_codes[index] = sub('0{2}$', '', sub('^0+','', ct))
  index = index + 1
}

#Make dataframe mapping full  census tracts to short codes
DT = data.frame(CT2010 = CTs_for_D20$CT2010, X2010.Census.Tract = unlist(short_codes))

#Read in raw data and replace short codes with full census tract names
data = read.csv('data/nyc_doe.csv')
correct_CTs = merge(data, DT, by='X2010.Census.Tract')
correct_CTs = correct_CTs[ , !(names(correct_CTs)=='X2010.Census.Tract')]

#Fill missing values (missing rows correspond to 0 counts)
years = unique(correct_CTs$School.Year)
grades = unique(correct_CTs$Grade.Level)
all_cts = union(unique(correct_CTs$CT2010), CTs_for_D20$CT2010)
full_coverage = expand.grid(School.Year = years, Grade.Level = grades, CT2010 = all_cts)
full_coverage = merge(full_coverage, correct_CTs, by=c("CT2010","Grade.Level", "School.Year"), all.x=TRUE)
full_coverage[is.na(full_coverage)] = 0

library(reshape)
wide = reshape(full_coverage, idvar=c('CT2010','School.Year'),timevar='Grade.Level', direction='wide')
wide = wide[order(wide$School.Year, as.character(wide$CT2010)),]
test = wide[(wide$School.Year %in% c(20092010,20102011)),]
train = wide[!(wide$School.Year %in% c(20092010,20102011)),]
wide2 = reshape(wide, idvar='CT2010', timevar='School.Year', direction='wide')

#Potentially drop tracks with 0 counts- note methodology seems questionable.
#total_student_count = rowSums(wide2[,!(names(wide2) %in% c('CT2010'))])
#wide$CT2010[which(total_student_count == 0)]
