#insulin_regulated_metabolites_project
#aim: To find metabolites asoociated with insulin regulated proteins (insulin regulted metabolites)
#data: using phosphoproteomic data from Humphrey et al, 2013 supp table S2 at 'http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3690479/?report=classic'
#first made a new spreadsheet with only relevent data from the above. This includes uniprot protein IDs and data columns relating to insulin/basal at 20 min insulin stimulation (9 data columns, 3 expts, 3 reps). Note that values are log2.
#exported spreadsheet as 'humphrey_20_min_data.csv'

## cleaned up data
#there are a lot of missing values denoted by 'NaN'.
humphrey_20_min_data=read.csv("insulin_regulated_metabolites_project1/humphrey_20_min_data.csv", sep=",", header=TRUE,na.string="NaN")
fix(humphrey_20_min_data) #fix displays the data table
dim(humphrey_20_min_data) #dim tell how many rows and columns there are: 37248  10
#remove rows with all missing values for data ie column 2 and onwards (column 1 is unprot ID)
#convert data to matrix so can remove rows
matrix_data=as.matrix(humphrey_20_min_data)
# first label rows for which col 2-10 are all Na
# the '1,' indicates to pass the function over rows (2 would indicate to pass it over a column direction, which would be one cell at a time).
row_entries_all_na = apply(matrix_data[,2:10], 1, function(x) {all(is.na(x))} ) 
# then collect all rows that are not all Na
# comma at the end makes it select all columns, not just the first eg if you type ">matrix_data[3]" you will get the 3rd row col 1 only same as "matrix_data[3,1]" but if you type ">matrix_data[3,]" you will get the entire 3rd row.
no_null_rows = matrix_data[!row_entries_all_na,]
#check how many rows remain
dim(no_null_rows) #  24756  10
#convert matrix back to .csv
fileout = write.csv(no_null_rows, file="insulin_regulated_metabolites_project1/no_null_rows.csv")
#check how many rows we have now
no_null_rows_data=read.csv("insulin_regulated_metabolites_project1/no_null_rows.csv", sep=",", header=TRUE,na.string="NaN")
fix(no_null_rows_data)
dim(no_null_rows_data) # rows 24756 cols 11 - an extra column[1] 'X' was added with the former row numbers, in .csv conversion?
#made a new matrix without the extra column
clean_data = no_null_rows_data[,2:11]
dim(clean_data) #  24756    10

##2 Average and visualise data
