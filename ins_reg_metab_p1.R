#insulin_regulated_metabolites_project
#aim: To find metabolites asoociated with insulin regulated proteins (insulin regulted metabolites)
#data: using phosphoproteomic data from Humphrey et al, 2013 supp table S2 at 'http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3690479/?report=classic'
#first made a new spreadsheet with only relevent data from the above. This includes uniprot protein IDs and data columns relating to insulin/basal at 20 min insulin stimulation (9 data columns, 3 expts, 3 reps). Note that values are log2.
#exported spreadsheet as 'humphrey_20_min_data.csv'
#there are a lot of missing values denoted by 'NaN'.
humphrey_20_min_data=read.csv("insulin_regulated_metabolites_project1/humphrey_20_min_data.csv", sep=",", header=TRUE,na.string="NaN")
fix(humphrey_20_min_data) #fix displays the data table
dim(humphrey_20_min_data) #dim tell how many rows and columns there are
#remove rows with all missing values for data ie column 2 and onwards (column 1 is unprot ID)
#convert data to matrix so can remove rows
matrix_data=as.matrix(humphrey_20_min_data)
#collect rows for which col 2-10 are not Na
nonullrows_matrix_data=matrix_data[!apply(matrix_data[,2:10],1,function(x) {all(is.na(x))} )]
#convert matrix back to .csv
nonullrows_humphrey_20_min_data=write.csv(nonullrows_matrix_data, file="insulin_regulated_metabolites_project1/nonullrows_humphrey_20_min_data.csv")