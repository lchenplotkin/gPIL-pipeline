'''
CSV FORMAT: 
Region pixels Number of pixels representing the region.
Region area Area representing the region (region pixels x pixel scale).
Area unit Region area unit (unit of the pixel scale).
Object count Number of objects located in the region (N/A if object splitting “ON”).
Object pixels Number of pixels representing objects in this region.
Object area Area representing objects in this region (object pixels x pixel scale).
Load Ratio of Object pixels to Region pixels (Object pixels / Region pixels).
'''


#import libraries

import os
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from functools import cmp_to_key
import matplotlib as mpl
import csv
import numpy as np

leftpath='1/leftcustom'
rightpath='1/rightcustom'
#specifying filepaths for the left and right brains

plt.rc('xtick', labelsize=5)#setting x-axis labels to prefered size
namelist=[]#names of regions containing pathology will be appended to this list
artifactlist=['Arcuate hypothalamic nucleus', 'Ventromedial hypothalamic nucleus', 'Median eminence', 'third ventricle']
artifactlist=[]

def compare(a,b): 
    if a[2]=='M' and b[2]=='F':
        return -1
    elif a[2]!=b[2]:
        return 1
    elif a[0]==b[0]:
        if a[3]<b[3]:
            return -1
        else:
            return 1
    elif a[0]=='d':
        return -1
    else:
        return 1
#a sorting key to sort the order of the presented sections. (males on left, with d2 genotype before g+)

def prepare(filepath,indexlist):
    filedict={}#A dictionary, the keys are the names of the brains (eg. D2F1, G+M2), and the values are the files containing the respective data
    braindict={}#A dictionary that will be populated with the data, (once again, keys are the names of the brains)
    maxload=0 #to scale for colormap

    i=0
    for subfolder in os.listdir(filepath):
        filedict[subfolder]=[]
        for entry in os.listdir(filepath+'/'+subfolder):
            filedict[subfolder].append(filepath+'/'+subfolder+'/'+entry)
        i+=1
    #populating file dictionary

    for brainfolder in filedict:
        braindict[brainfolder]=[]
        for sectionfile in filedict[brainfolder]:
            with open(sectionfile,'r') as f:
                templist=list(csv.reader(f,delimiter=';'))
                templist.pop(0)
                #reading in file

                for line in templist:
                    if templist.index(line) not in indexlist:
                        indexlist.append(templist.index(line))
                #taking the indices of the regions that have a non-negligible amount of pathology and appending them to a list of such indices

    indexlist.sort()
    for brain in braindict:
        totalarea=[]#this list will be used to take the total area of each region across all csvs in each brain folder, which will be used to take a weighted average
        for index in indexlist:
            braindict[brain].append(0)
            totalarea.append(0)
        #appending to lists so that they are the right length (the number of indices)

        for section in filedict[brain]: #ie. "For each csv in each folder"
            with open(section,'r') as f: 
                templist=list(csv.reader(f,delimiter=';'))
                templist.pop(0)
                #opens the data from the csv into a list

                i=0
                for index in indexlist:
                    totalarea[i]+=float(templist[index][2])
                    i+=1
                #adding to total area

                while 0 in totalarea:
                    totalarea[totalarea.index(0)]+=0.000000000001
                #this ensures that python does not divide by 0 when taking an average

        for section in filedict[brain]:
            with open(section,'r') as f:
                templist=list(csv.reader(f,delimiter=';'))
                templist.pop(0)
                indexlist.sort()
                #same as above

                i=0
                for index in indexlist:
                    name=templist[index][0]
                    if name not in namelist:
                        namelist.insert(i,name)
                    modifier=float(templist[index][2])/totalarea[i]
                    #some images have more or less of different regions present, so we take a weighted average

                    braindict[brain][i]+=float(templist[index][8])*modifier #applying weighted average
                    if float(templist[index][8])*modifier>maxload:
                        maxload=float(templist[index][8])*modifier#the maxload is used to scale on the colorbar
                    i+=1

                    


    sorted_keys=sorted(list(braindict.keys()),key=cmp_to_key(compare))
    sorted_vals=[]
    for key in sorted_keys:
        sorted_vals.append(braindict[key])
    #using the sorting function defined in the beginning to order the labels in the right way for the plot (males on left, with G+ genotype on right)

    return([sorted_keys,sorted_vals,maxload,indexlist])#returning a list: #the keys(names of brains), values(big lists with each brain's data), the maxload (to be used for scaling the colorplot), the list of indices corresponding to regions that have non-negligible amounts of pathology in one or more brains

def plot_test(cmaps,pylistdat,pylistnames,samplenames,filename,maxlist):
    datalist=[]

    i=0
    for data in pylistdat:
        datalist.append(np.swapaxes(np.asarray(data),0,1))
        i+=1

    names=np.asarray(pylistnames)
    splnames=np.asarray(samplenames)

    n = len(datalist)

    fig, axs = plt.subplots(1,n, figsize=(20,35), sharey=True)
   
    plt.xticks([])
    plt.yticks([])

    i=0
    for [ax, data] in zip(axs.flat, datalist):
        tempmax=maxlist[i]

        psm = ax.pcolormesh(data, cmap=cmaps[i], rasterized=False, vmin=0, vmax=tempmax) 
        fig.colorbar(psm, ax=ax)
        #creating the plot

        ytickslist=np.arange(len(names))
        xtickslist=np.arange(len(samplenames))
        ax.set_xticks([float(n)+0.5 for n in xtickslist],splnames,minor=True)
        ax.set_yticks([float(n)+0.5 for n in ytickslist],names,minor=False)
        #setting up ticks and formating ticks to be in the middle of each square, rather than on the left
        i+=1
    
    fig.subplots_adjust(left=0.3)
    plt.suptitle('Pathology load across 10 different mice')
    plt.savefig(filename)
    plt.show()


leftrun1=prepare(leftcustom,[])
leftkeys=leftrun1[0]
leftvals=leftrun1[1]
leftmax=leftrun1[2]
indlist=leftrun1[3]

'''
rightrun=prepare(rightcustom,indlist)
rightkeys=rightrun[0]
rightvals=rightrun[1]
rightmax=rightrun[2]
indlist=rightrun[3]
'''

viridis=mpl.colormaps['viridis'].resampled(256)
inferno=mpl.colormaps['inferno'].resampled(256)
magma=mpl.colormaps['magma'].resampled(256)
#choose whichever schemes look good, can refer to matplotlib.pycolormesh documentation

plot_test([viridis,magma],[leftvals,leftvals],namelist,leftkeys,"customleft.png",[leftmax,leftmax])
#exporting image
