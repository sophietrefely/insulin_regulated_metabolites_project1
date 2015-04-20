"""insulin_regulated_metabolites_project

aim: To find metabolites associated with insulin regulated proteins (insulin regulted metabolites)
data: using phosphoproteomic data from Humphrey et al, 2013 supp table S2 at 'http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3690479/?report=classic'
first made a new spreadsheet with only relevent data from the above. This includes uniprot protein IDs and data columns relating to insulin/basal at 20 min insulin stimulation (9 data columns, 3 expts, 3 reps). Note that values are log2 (this allows fold change in either direction to be compared).
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

print(clean_data_list[1:5])

##2 create z-factors for each line of data

#Find mean for each line (phosphosite) and assign to uniprot ID, ie column [0]. Add standard deviation (stdev)
import statistics
stats_data_list = []
for line in clean_data_list:
    stats_data_line = [line[0]]
    data = line[1:]  #need to convert strings to floats here to get mean etc
    float_data = [float(item) for item in data]
    mean = statistics.mean(float_data)
    n = len(data)
    #can't do stdev for 1 value. Therefore set stdev to 0
    if n==1:
        stdev = 0
    else:
        stdev = statistics.stdev(float_data)
    #extend instead of append twice so [0] = uniprot id, [1] = mean, [2] = stdev, [3] = n.
    stats_data_line.extend([mean,stdev,n])
    stats_data_list.append(stats_data_line)

print(stats_data_list[1:5])
# Find min, max and median
import numpy as np
means1 = []
for item in stats_data_list:
    mean = item[1]
    means1.append(mean)   
minimum = np.amin(means1)
print('Stats for All Data')
print('minimum is', minimum)
maximum = np.amax(means1)
print('maximum is', maximum)
median = np.median(means1)
print('median is', median)
sites = len(stats_data_list)
print(sites, 'phosphosites')

#plot the data to look at close to 0 vs regulated.
import matplotlib.pyplot as plt
plt.hist([item[1] for item in stats_data_list if item[1]<10], bins=100, color = '0.75', label='All data')
plt.title('All Data')
plt.xlabel('Mean log2(insulin/basal) phosphorylation')
plt.ylabel('# of phosphorylation sites')

#plot control sites as dot plot that can overlay onto the histogram.
#Used controls used as examples for validation in Humphrey et al:
# positive control (increase w insulin): AKT1/2 T308, TSC2 S981, TSC2 T1462, PRAS40 T246, GSK3a S21,
# negative control (decrease w insulin):

poscontrol = [1, 2.0, 3.0, 3.5, 5.0, 6.0]
yaxisas1= [1000,1000,1000,1000,1000,1000]
negcontrol = [-1, -2, -4, -2.5, - 5, -2.25]
plt.plot(poscontrol, yaxisas1, marker='o', linestyle='None', color='b', label='Positive controls')
plt.plot(negcontrol, yaxisas1, marker='o', linestyle='None', color='r', label='Negative controls')
plt.legend()
plt.show()

#I only want high confidence proteins, therefore removed sites where n=1 from analysis
#some proteins will have regulated and unregulated sites. Therefore, if a protein has one regulated site it will be designated regulated.
# Another way to find 'non-insulinregulated proteins' would be to use the proteome list generated in this study and sutract reulated proteins from it.
stringent_data_list = []
for item in stats_data_list:
    if item[3]>1: #ie n>1
        stringent_data_list.append(item)

print(stringent_data_list[1:5])

stringent_means = []
for item in stringent_data_list:
    mean = item[1]
    stringent_means.append(mean)   

minimum = np.amin(stringent_means)
print('Stats for n>1 Data')
print('minimum is', minimum)
maximum = np.amax(stringent_means)
print('maximum is', maximum)
median = np.median(stringent_means)
print('median is', median)
sites = len(stringent_data_list)
print(sites, 'phosphosites')                
    
#plot the n>1 data 'stringent_data_list'
plt.hist([item[1] for item in stringent_data_list], bins=100, color = '0.75', label='Stringent (n>1) data')
plt.title('Stringent (n>1) data')
plt.xlabel('Mean log2(insulin/basal) phosphorylation')
plt.ylabel('# of phosphorylation sites')

#plot control sites as dot plot that can overlay onto the histogram.
#Used controls used as examples for validation in Humphrey et al:
# positive control (increase w insulin): AKT1/2 T308, TSC2 S981, TSC2 T1462, PRAS40 T246, GSK3a S21,
# negative control (decrease w insulin):

poscontrol = [1, 2.0, 3.0, 3.5, 5.0, 6.0]
yaxisas1= [1000,1000,1000,1000,1000,1000]
negcontrol = [-1, -2, -4, -2.5, - 5, -2.25]
plt.plot(poscontrol, yaxisas1, marker='o', linestyle='None', color='b', label='Positive controls')
plt.plot(negcontrol, yaxisas1, marker='o', linestyle='None', color='r', label='Negative controls')
plt.legend()
plt.show()

#Calculate Median Absolute Deviation (MAD) +/-2.5
from numpy import median, absolute

def mad(data, axis=None):
    return median(absolute(data - median(data, axis)), axis)
stringent_data_mad = mad(stringent_means)
print('MAD:', stringent_data_mad)

2.5_MAD = stringent_data_mad*2.5
range_pos_MAD = median + 2.5_MAD
range_neg_MAD = median - 2.5_MAD
