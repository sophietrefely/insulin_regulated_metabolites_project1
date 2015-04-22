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

# Find min, max, median etc
import numpy as np
means1 = []
for item in stats_data_list:
    mean = item[1]
    means1.append(mean)
    
minimum = np.amin(means1)
print('Stats for All Data:')
print('minimum is', minimum)
maximum = np.amax(means1)
print('maximum is', maximum)
all_median = float(np.median(means1))
print('median is', all_median)
mean = np.mean(means1)
print('mean is', mean)
stddev= np.std(means1)
print('stddev is', stddev)
sites = len(stats_data_list)
print(sites, 'phosphosites')

#plot std dev for All data
import matplotlib.pyplot as plt
plt.hist([item[2] for item in stats_data_list], bins=100, color = '0.75', label='std dev')
plt.title('Std Dev for All Data')
plt.xlabel('Std Dev')
plt.ylabel('# of phosphorylation sites')
plt.show()

# Only want high confidence proteins, therefore removed sites where n=1 from analysis, ie non-repeated results
#Reminder: [0] = uniprot id, [1] = mean, [2] = stdev, [3] = n.
stringent_data_list = []
for item in stats_data_list:
    if item[3]>1: #ie n>1
        stringent_data_list.append(item)
print(stringent_data_list[1:5])

#Determined min, max, mean etc of data
stringent_means = []
for item in stringent_data_list:
    mean = item[1]
    stringent_means.append(mean)

print('Stats for n>1 Data:')
stringent_minimum = np.amin(stringent_means)
print('minimum is', stringent_minimum)
stringent_maximum = np.amax(stringent_means)
print('maximum is', stringent_maximum)
stringent_median = float(np.median(stringent_means))
print('median is', stringent_median)
stringent_mean = np.mean(stringent_means)
print('mean is', stringent_mean)
stringent_stddev= np.std(stringent_means)
print('stddev is', stringent_stddev)
sites = len(stringent_data_list)
print(sites, 'phosphosites')

###convert 'Mean log2(insulin/basal) phosphorylation' into z-scores. For each mean -item[1]- (z=(item[1]-overallmean)/overallstddev))
##for item in stringent_data_list:
##    z = (item[1]-stringent_mean)/stringent_stddev
##    item[1] = z
##print(stringent_data_list[1:5])

#Calculate Median Absolute Deviation (MAD) +/-2.5 from stringent_data_list means
from numpy import median, absolute

def mad(data, axis=None):
    return median(absolute(data - median(data, axis)), axis)
##z_score_list = []
##for item in stringent_data_list:
##    z=item[1]
##    z_score_list.append(z)
##    
##z_score_mad = float(mad(z_score_list))
##print('MAD:', z_score_mad)
##
##MAD2_5 = 2.5*(z_score_mad)
##print('2.5*MAD:', MAD2_5)
MAD = float(mad(stringent_means))
print('Stringent MAD', MAD)

MAD2_5 = 2.5*(MAD)
print('2.5*Stringent MAD:', MAD2_5)

all_data_mad = float(mad(means1))
print('All data MAD:', all_data_mad)

all_data_mad2_5 = 2.5*(all_data_mad)
print('2.5*All data MAD:', all_data_mad2_5)

#plot data to compare original and stringent and show cutoff lines for MAD and positive control proteins
import matplotlib.pyplot as plt
plt.hist([item[1] for item in stats_data_list if item[1]<10], bins=100, color = '0.75', label='All data')
plt.title('Data Processing')
plt.xlabel('Mean log2(insulin/basal) phosphorylation')
plt.ylabel('# of phosphorylation sites')
#plot the n>1 data means in 'stringent_data_list'
plt.hist([item[1] for item in stringent_data_list], bins=100, color = 'b', label='Stringent (n>1) data')
plt.xlabel('Mean log2(insulin/basal) phosphorylation')
plt.ylabel('# of phosphorylation sites')
#overlay MAD*2.5 barrier lines
plt.plot([stringent_median+MAD2_5, stringent_median+MAD2_5],[0, 5000], 'r--')
plt.plot([stringent_median-MAD2_5, stringent_median-MAD2_5],[0, 5000], 'r--', label='Stringent data +/-2.5*MAD')
#overlay MAD*1 barrier lines
plt.plot([all_median+MAD, all_median+MAD],[0, 5000], 'm--')
plt.plot([all_median-MAD, all_median-MAD],[0, 5000], 'm--', label='Stringent data +/-MAD')
#overlay All data MAD*2.5 barrier lines
plt.plot([all_data_mad2_5, all_data_mad2_5],[0, 5000], 'r--')
plt.plot([-all_data_mad2_5, -all_data_mad2_5],[0, 5000], 'r--', label='All data +/-2.5*MAD')
#overlay All data MAD*1 barrier lines
plt.plot([all_data_mad, all_data_mad],[0, 5000], 'm--')
plt.plot([-all_data_mad, -all_data_mad],[0, 5000], 'm--', label='All data +/-MAD')
#overlay control phosphosites:
#Used controls used as examples for validation in Humphrey et al:
# positive control (known Akt pathway sites increase w insulin): AKT1/2 T308, TSC2 S981, TSC2 T1462, PRAS40 T246, GSK3a S21, GSK3b S9, ENOS S1177, PFKFB2 S466, AS160 T588, AS160 S642, AS160 S318, AS250 T715
poscontrol = [1, 2.0, 3.0, 3.5, 5.0, 6.0]
yaxisas1= [100,100,100,100,100,100]
plt.plot(poscontrol, yaxisas1, marker='o', linestyle='None', color='r', label='Positive controls')
#annotate divisions for insulin-regulated/non-regulated lists
plt.annotate('non-regulated sites', xy=(0, 2500), xytext=(-5, 2600),
            arrowprops=dict(facecolor='k', shrink=0.05),
            )
plt.text(-7, 2000, 'negatively regulated sites')
plt.text(3, 2000, 'positively regulated sites')
plt.legend()
plt.show()

#From the graph, we are better off using 'All data' than 'Stringent data' because both data sets behave very similarly and you get more information by including all.
#Therefore, for 'All data', create lists of positive and negatively changed sites with mean log2(insulin/basal) beyond 2.5*MAD from 0
#Reminder: [0] = uniprot id, [1] = mean, [2] = stdev, [3] = n.
all_pos_reg = []
for item in stats_data_list:
    if item[1] > all_median+all_data_mad2_5:
        all_pos_reg.append(item)
all_pos_uniprots = []
for item in all_pos_reg:
    uniprot = item [0]
    all_pos_uniprots.append(uniprot)
print('all_pos_reg_sites:', len(all_pos_reg)) 

all_neg_reg = []
for item in stats_data_list:
    if item[1] < all_median-all_data_mad2_5:
        all_neg_reg.append(item)
all_neg_uniprots = []
for item in all_neg_reg:
    uniprot = item [0]
    all_neg_uniprots.append(uniprot)
print('all_neg_reg_sites:', len(all_neg_reg))     

#Create list of non-regulated sites within 1*MAD from 0
all_non_reg = []
for item in stats_data_list:
    if item[1] < all_median+all_data_mad and item[1] > all_median-all_data_mad:
        all_non_reg.append(item)
print('all_non_reg_sites:', len(all_non_reg))       

#some proteins will have regulated and unregulated sites. If a protein has one regulated site it will be designated regulated. Therefore remove from non-reg list
true_non_reg = []
for item in all_non_reg:
    if item[0] not in all_pos_uniprots and item[0] not in all_neg_uniprots:
      true_non_reg.append(item)
print('true_non_reg_sites:', len(true_non_reg))
true_non_uniprots = []
for item in true_non_reg:
    uniprot = item[0]
    true_non_uniprots.append(uniprot)

#proteins will be represented within lists multple time (as they can contain several regulated/non-regulatred sites). Therefore must remove duplication within lists.
pos_reg_proteins = []
for item in all_pos_uniprots:
    if item not in pos_reg_proteins:
        pos_reg_proteins.append(item)
print('pos_reg_proteins:', len(pos_reg_proteins))

neg_reg_proteins = []
for item in all_neg_uniprots:
    if item not in neg_reg_proteins:
        neg_reg_proteins.append(item)
print('neg_reg_proteins:', len(neg_reg_proteins))

non_reg_proteins = []
for item in true_non_uniprots:
    if item not in non_reg_proteins:
        non_reg_proteins.append(item)
print('non_reg_proteins:', len(non_reg_proteins))

#If a protein is in both the positive and negative list, positive trumps negative. Therefore remove from negative list
neg_reg_proteins_only = []
for item in neg_reg_proteins:
    if item not in pos_reg_proteins:
        neg_reg_proteins_only.append(item)
print('neg_reg_proteins_only:', len(neg_reg_proteins_only))
        
