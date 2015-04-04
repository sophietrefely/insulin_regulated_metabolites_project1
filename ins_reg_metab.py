"""insulin_regulated_metabolites_project

aim: To find metabolites asoociated with insulin regulated proteins (insulin regulted metabolites)
data: using phosphoproteomic data from Humphrey et al, 2013 supp table S2 at 'http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3690479/?report=classic'
first made a new spreadsheet with only relevent data from the above. This includes uniprot protein IDs and data columns relating to insulin/basal at 20 min insulin stimulation (9 data columns, 3 expts, 3 reps). Note that values are log2.
exported spreadsheet as 'humphrey_20_min_data.csv'
"""
##1 cleaned up data
#there are a lot of missing values denoted by 'NaN'.

original_data= open("humphrey_20_min_data.csv", 'rU') 

# split creates a list out of the line [<first>, <second>] field
# [0] is gene uniprot ID, [1] first expt log2 ins/bas....[9] 9th expt. log2 ins/bas.

original_data_list = []
for line in original_data:
    stripped_line = line.rstrip()   # remove the newline from the line, it has one by default
    line_split = stripped_line.split(',') #fields are comma seperated
    original_data_list.append(line_split)

#remove lines that contain only missing data ie lines with 'NaN' for [1]...[9]
#make new data_list

data_list = []
for line in original_data_list:
    #call all function 'all' () on list comprehension []
    all_items_NaN = all([item == 'NaN' for item in line[1:9]])
    if not all_items_NaN:
        data_list.append(line)

#Remove NaNs from data_list
clean_data_list = []
for line in data_list:
    new_line = [] #make a new list of values without NaNs
    for item in line:
        if not item == 'NaN':
            new_line.append(item)
    clean_data_list.append(new_line)
print(clean_data_list[1])

##2 create z-factors for each line of data

#First unlog the data ie convert each integer to 2**(integer)
unlogged_data_list = []
for line in clean_data_list:
    unlogged_line = []
    for item in line[1:]:  #not uniprot id, which is line[0]
        integer_item = float(item)
        unlogged_integer_item = 2**(integer_item)
        unlogged_line.append(unlogged_integer_item)
    unlogged_data_list.append(unlogged_line)

print(unlogged_data_list[1])

#Plot data as mean for each line with standard devation
