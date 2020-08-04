#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as st
import numpy as np
import seaborn as sns


# In[2]:


mouse_metadata_path = os.path.join(r"/Users/Munit/Desktop/MatplotlibHW/Resources/Mouse_metadata.csv")
study_results_path = os.path.join(r"/Users/Munit/Desktop/MatplotlibHW/Resources/Study_results.csv")


# In[3]:


mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)
mouse_metadata.head()


# In[4]:


study_results.head()


# In[5]:


metadata_results = pd.DataFrame.merge(study_results, mouse_metadata, on="Mouse ID", how="outer")
metadata_results


# In[6]:


#checking number of mice
len(metadata_results["Mouse ID"].value_counts())


# In[7]:


#Getting the duplicate mice by ID number that shows up for Mouse ID and Timepoint
metadata_results[metadata_results.duplicated(['Mouse ID', 'Timepoint'], keep=False)]


# In[9]:


# Optional: Get all the data for the duplicate mouse ID. 
duplicate = metadata_results[metadata_results.duplicated(['Mouse ID', 'Timepoint'], keep=False)]
duplicate.info("Mouse ID")


# In[10]:


# Create a clean DataFrame by dropping the duplicate mouse by its ID.
metadata_results_df = metadata_results.drop_duplicates(subset="Mouse ID")
metadata_results_df


# In[11]:


# Checking the number of mice in the clean DataFrame.
len(metadata_results["Mouse ID"].value_counts())


# In[12]:


metadata_results_sort = metadata_results.sort_values(["Tumor Volume (mm3)"], ascending=True)
metadata_results_sort.head()


# In[14]:


#Summary Statistics
# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen
drug_regimen_grp = metadata_results_sort.groupby(["Drug Regimen"])
drug_regimen_grp
regimen_median = drug_regimen_grp['Tumor Volume (mm3)'].median()
regimen_median
regimen_var = drug_regimen_grp['Tumor Volume (mm3)'].var()
regimen_std = drug_regimen_grp['Tumor Volume (mm3)'].std()
regimen_sem = drug_regimen_grp['Tumor Volume (mm3)'].sem()
regimen_mean = drug_regimen_grp['Tumor Volume (mm3)'].mean()


# In[15]:


summary_stats = pd.DataFrame({"Mean Tumor Volume": regimen_mean, "Median Tumor Volume": regimen_median, 
                              "Tumor Volume Variance":regimen_var, "Tumor Volume Std. Dev": regimen_std, 
                              "Tumor Volume Std. Err": regimen_sem})
summary_stats


# In[16]:


regimen_data_points = metadata_results.groupby(["Drug Regimen"]).count()["Mouse ID"]
regimen_data_points


# In[17]:


#Bar and Pie Charts
# Generate a bar plot showing the total number of mice for each treatment throughout the course of the study using pandas.
regimen_data_points.plot(kind='bar', figsize = (9,6))
plt.xlabel("Drug Regimen")
plt.ylabel("Number of Data Points")
plt.show()
plt.tight_layout() 


# In[18]:


# Generate a bar plot showing the total number of mice for each treatment throughout the course of the study using pyplot.
x_axis = np.arange(len(regimen_data_points))
tick_locations = [x for x in x_axis]
dataframe = pd.DataFrame(regimen_data_points)
dataframe.plot.bar(legend=False,rot='vertical',figsize = (9,6))
plt.ylabel("Number of Data Points")
plt.show()


# In[19]:


# Generate a pie plot showing the distribution of female versus male mice using pandas
gender_grp = metadata_results.groupby(["Sex"])
gender_grp['Sex'].value_counts()
gender_grp['Sex'].value_counts().sum()
female_grp = (gender_grp['Sex'].value_counts()/gender_grp['Sex'].value_counts().sum())*100
female_grp
data = {'Sex': [49, 51]}
gender_grp_df = pd.DataFrame(data,columns=['Sex'], index=['Female', 'Male'])
gender_grp_df        
gender_grp_df.plot.pie(y='Sex', figsize=(5, 5), autopct="%1.1f%%")
plt.show()


# In[20]:


# Generate a pie plot showing the distribution of female versus male mice using pyplot
gender = ["Female", "Male"]
percentages = ['49', '51']
plt.pie(percentages, labels=gender,autopct="%1.1f%%")
plt.pie


# In[21]:


#Quartiles, Outliers and Boxplots
# Calculate the final tumor volume of each mouse across four of the treatment regimens: Capomulin, Ramicane, Infubinol, and Ceftamin

capomulin = metadata_results.loc[metadata_results["Drug Regimen"] == "Capomulin",:]
ramicane = metadata_results.loc[metadata_results["Drug Regimen"] == "Ramicane", :]
infubinol = metadata_results.loc[metadata_results["Drug Regimen"] == "Infubinol", :]
ceftamin = metadata_results.loc[metadata_results["Drug Regimen"] == "Ceftamin", :]


# In[22]:


#group data by drug regimen and mouse id to capture last tumor volume
cap_last = capomulin.groupby('Mouse ID').max()['Timepoint']
cap_last_vol = pd.DataFrame(cap_last)
cap_merge = pd.merge(cap_last_vol, metadata_results, on=("Mouse ID","Timepoint"),how="left")

tumor_vol1 = cap_merge["Tumor Volume (mm3)"]
# Calculate the IQR and quantitatively determine if there are any potential outliers. 
quartiles = tumor_vol1.quantile([.25,.5,.75])
lowerq = quartiles[0.25]
upperq = quartiles[0.75]
iqr = upperq-lowerq
lower_bound = lowerq - (1.5*iqr)
upper_bound = upperq + (1.5*iqr)
#print(f"Capomulin potential outliers: below {lower_bound} and above {upper_bound}")


# In[23]:


#group data by drug regimen and mouse id to capture last tumor volume
ram_last = ramicane.groupby('Mouse ID').max()['Timepoint']
ram_last_vol = pd.DataFrame(ram_last)
ram_merge = pd.merge(ram_last_vol, metadata_results, on=("Mouse ID","Timepoint"),how="left")

tumor_vol2 = ram_merge["Tumor Volume (mm3)"]
# Calculate the IQR and quantitatively determine if there are any potential outliers. 
quartiles2 = tumor_vol2.quantile([.25,.5,.75])
lowerq2 = quartiles2[0.25]
upperq2 = quartiles2[0.75]
iqr2 = upperq2-lowerq2
lower_bound2 = lowerq2 - (1.5*iqr2)
upper_bound2 = upperq2 + (1.5*iqr2)
#print(f"Ramicane potential outliers: below {lower_bound2} and above {upper_bound2}")


# In[24]:


#group data by drug regimen and mouse id to capture last tumor volume
inf_last = infubinol.groupby('Mouse ID').max()['Timepoint']
inf_last_vol = pd.DataFrame(inf_last)
inf_merge = pd.merge(inf_last_vol, metadata_results, on=("Mouse ID","Timepoint"),how="left")
inf_merge.head(10)

tumor_vol3 = inf_merge["Tumor Volume (mm3)"]
#Calculate the IQR and quantitatively determine if there are any potential outliers. 
quartiles3 = tumor_vol3.quantile([.25,.5,.75])
lowerq3 = quartiles3[0.25]
upperq3 = quartiles3[0.75]
iqr3 = upperq3-lowerq3
lower_bound3 = lowerq3 - (1.5*iqr3)
upper_bound3 = upperq3 + (1.5*iqr3)
#print(f"Infubinol potential outliers: below {lower_bound3} and above {upper_bound3}")


# In[25]:


#group data by drug regimen and mouse id to capture last tumor volume
cef_last = ceftamin.groupby('Mouse ID').max()['Timepoint']
cef_last_vol = pd.DataFrame(cef_last)
cef_merge = pd.merge(cef_last_vol, metadata_results, on=("Mouse ID","Timepoint"),how="left")
cef_merge.head(10)

tumor_vol4 = cef_merge["Tumor Volume (mm3)"]
# Calculate the IQR and quantitatively determine if there are any potential outliers. 
quartiles4 = tumor_vol4.quantile([.25,.5,.75])
lowerq4 = quartiles4[0.25]
upperq4 = quartiles4[0.75]
iqr4 = upperq4-lowerq4
lower_bound4 = lowerq4 - (1.5*iqr4)
upper_bound4 = upperq4 + (1.5*iqr4)
#print(f"Ceftamin potential outliers: values below {lower_bound4} and above {upper_bound4}")


# In[26]:


print(f"Capomulin potential outliers: below {lower_bound} and above {upper_bound}")
print(f"Ramicane potential outliers: below {lower_bound2} and above {upper_bound2}")
print(f"Infubinol potential outliers: below {lower_bound3} and above {upper_bound3}")
print(f"Ceftamin potential outliers: values below {lower_bound4} and above {upper_bound4}")


# In[27]:


# Generate a box plot of the final tumor volume of each mouse across four regimens of interest

plot_data = [tumor_vol1, tumor_vol2, tumor_vol3, tumor_vol4]

fig1, ax1 = plt.subplots()
ax1.set_ylabel('Final Tumor Volume (mm3)')
ax1.set_xlabel('Drug Regimen')

ax1.boxplot(plot_data, labels=["Capomulin","Ramicane","Infubinol","Ceftamin",])
plt.show()
#plt.savefig('boxplot')

# top_regimen = metadata_results[metadata_results['Drug Regimen'].isin(['Capomulin', 'Ramicane'
#                                                 , 'Infubinol', 'Ceftamin'])]
# top_regimen = top_regimen.sort_values(['Timepoint'], ascending=True)
# top_regimen_df = top_regimen[['Drug Regimen', 'Mouse ID', 'Timepoint', 'Tumor Volume (mm3)']]


# In[28]:


#Line and Scatter Plots
# Generate a line plot of time point versus tumor volume for a mouse treated with Capomulin

time_vs_tumor = metadata_results[metadata_results["Mouse ID"].isin(["g316"])]
time_vs_tumor_df = time_vs_tumor[['Mouse ID', 'Timepoint', 'Tumor Volume (mm3)']]
time_vs_tumor_df
line_plot = time_vs_tumor_df.reset_index()
line_plot
line_plot2 = line_plot[['Mouse ID', 'Tumor Volume (mm3)']]
line_plot2.plot.line()
plt.title('Capomulin Treatment of mouse g316')
plt.xlabel('Timepoint(days)')
plt.ylabel('Tumor Volume (mm3)')
plt.show()


# In[30]:


# Generate a scatter plot of mouse weight versus average tumor volume for the Capomulin regimen

cap_average = capomulin.groupby(['Mouse ID']).mean()
plt.scatter(cap_average['Weight (g)'],cap_average['Tumor Volume (mm3)'])
plt.xlabel('Weight (g)')
plt.ylabel('Average Tumor Volume (mm3)')

#plt.savefig('scatterplot')
plt.show()


# In[31]:


#Correlation and Regression
# Calculate the correlation coefficient and linear regression model for mouse weight and average tumor volume for the Capomulin regimen
x = cap_average ["Weight (g)"]
y = cap_average ["Tumor Volume (mm3)"]
slope, intercept, rvalue, pvalue, stderr = st.linregress(x, y)
regress_values = x * slope + intercept
plt.scatter(x, y)
plt.plot(x,regress_values,'r-')
plt.xlabel("Weight (g)")
plt.ylabel("Average Tumor Volume (mm3)")
plt.show()


# In[32]:


pearson_coef, p_value = st.pearsonr(cap_average["Weight (g)"], cap_average["Tumor Volume (mm3)"])
print(f"The correlation between mouse weight and the average tumor volume is {pearson_coef}")


# In[ ]:




