#script to track where mitosis occurs per generations. Construct 3D graph (in 2 dimensions, 3rd dimension depicted as pointcolour. plot position of mitosis events in 2D space, colour code with time. Using shades of grey for ease of viewing.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ast
from matplotlib import cm
import math
from numpy.polynomial.polynomial import polyfit
import csv

tLENHIGH=100
tLENLOW=20

mitosisData = pd.read_csv("5embtrack_list.csv")
'''
mitosisSUB=mitosisData[mitosisData["GEN"].str.match("\[3")] #take your generation of interest. reduce dataframe to contain only nuclei of yoru genreation of interest, and those that are registered as dividing (i.e the ones which have a generation + 1 encoded in their track.
mitosisSUB=mitosisSUB[mitosisSUB['GEN'].str.contains('4')] #find if it is tracked to mitosis. this is just going to be generation of interest+1 at the end of the string/list.  these 2 values are important as they are modifiable parameters which decide which generation you are finding the time point for mitosos for.
'''

mitosisSUB=mitosisData

mitTimeList=[]
for a in range(0,len(mitosisSUB)):
	valueInList=mitosisSUB["GEN"].iloc[a]
	genToList=ast.literal_eval(valueInList)
	xLOC=mitosisSUB["X"].iloc[a]
	xVAL=ast.literal_eval(xLOC)
	yLOC=mitosisSUB["Y"].iloc[a]
	yVAL=ast.literal_eval(yLOC)
	zLOC=mitosisSUB["Z"].iloc[a]
	zVAL=ast.literal_eval(zLOC)
	tLOC=mitosisSUB["T"].iloc[a]
	tVAL=ast.literal_eval(tLOC)
	nucDIAMLOC=mitosisSUB["nuc_diam"].iloc[a]
	nucDIAMVAL=ast.literal_eval(nucDIAMLOC)
	cellDIAMLOC=mitosisSUB["cell_diam"].iloc[a]
	cellDIAMVAL=ast.literal_eval(cellDIAMLOC)	
	ncDIAMLOC=mitosisSUB["NC_ratio"].iloc[a]
	ncDIAMVAL=ast.literal_eval(ncDIAMLOC)
	famIDLOC=mitosisSUB["famID"].iloc[a]
	famIDVAL=ast.literal_eval(famIDLOC)
	datafileLOC=mitosisSUB["datafile"].iloc[a]
	genLOC=mitosisSUB["GEN"].iloc[a]	
	genVAL=ast.literal_eval(genLOC)
	#datafileVAL=ast.literal_eval(datafileLOC)	
	genLEN=len(genToList)
	print("{},{},{},{}".format(xVAL[-1],yVAL[-1],zVAL[-1],tVAL[-1],genLEN))
	tLEN=tVAL[-1]-tVAL[1]+1
	max_nuc_diam_ = int(max(nucDIAMVAL[:]))
	max_cell_diam_ = int(max(cellDIAMVAL[:]))*2
	if max_cell_diam_ == 0:
		maxNCR=0
	else:
		maxNCR=max_nuc_diam_/max_cell_diam_
	
	plotcrd=(xVAL[-1],yVAL[-1],zVAL[-1],int(tVAL[-1]),genLEN,tLEN,max_nuc_diam_ ,max_cell_diam_,maxNCR,datafileLOC,famIDVAL[1],int(max(genVAL))) #create tuple with all the coordinates we could need, x,y,z,T,gen
	mitTimeList.append(plotcrd)
print(mitTimeList)

range_finder = []
for b in mitTimeList:
	if b[5] > 50 or b[5] < tLENHIGH:
		range_finder.append(b[3])
		#range_finder[i]=range_finder[i]
	else:
		continue
rangemin=min(range_finder)
rangemax=max(range_finder)

def rangeNorm(timevalue,rangemin,rangemax):
	#normVal=((timevalue-(rangemin-1))/(rangemax-(rangemin-1)))/254
	normVal=((timevalue-rangemin)/(rangemax-rangemin))
	#normVal = timevalue-rangemin
	print("normVal",normVal)
	return normVal
'''
for b in mitTimeList:
	if b[5] < tLENLOW or b[5] > tLENHIGH:
		continue
	print("b= ", b)
	plt.scatter(b[0],b[2],c=cm.gist_yarg(b[3]),s=(b[6]**2)/3) #need to multiply the mitosis time value by 15 as the colour pallette doesn't register the small time differences. scaling up puts distance between time points and these register better on colour scale. make points nice and big to emphasise colours.

'''
x_coord_list=[]
y_coord_list=[]
with open("scatter_outs_ALL.csv", "w", newline="") as scatterO:
	plotwrite = csv.writer(scatterO)
	plotwrite.writerow(["x","y","z","mit_time","genLen","CC_LEN","nucdiam","celldiam","ncratio","datafile", "id","generation"])
	for b in mitTimeList:
		#if b[5] < 20 or b[5] > tLENHIGH:
		if b[8]<0.01 or b[8]>0.85 or b[5] < 20 or b[5] > tLENHIGH:
			continue
		print("b= ", b)
		plt.scatter(b[5],b[8],c='black',s=40)
		x_coord_list.append(b[0])
		y_coord_list.append(b[8])
		plotwrite.writerow([b[0],b[1],b[2],b[3],b[4],b[5],b[6],b[7],b[8],b[9],b[10],b[11]])
		
		#plt.scatter(b[0],b[1],c=cm.binary(rangeNorm(b[3],rangemin,rangemax)-0.3),s=(b[5]*b[5])/10) 
		#need to multiply the mitosis time value by 15 as the colour pallette doesn't register the small time differences. scaling up puts distance between time points and these register better on colour scale. make points nice and big to emphasise colours.
#b,m=polyfit(np.array(x_coord_list),np.array(y_coord_list),1)
#plt.plot(np.array(x_coord_list),m*np.array(y_coord_list) + b)	

	
plt.title("Mitosis points in the transition from 512->1000cell stage (XY plane, dataset: 20210310)")	
plt.xlim(0,100)
plt.ylim(0,1)
plt.ylabel("cell cycle length (frames)")
plt.xlabel("x coordinate (microns)")
#plt.legend()
#plt.grid(True)	
#plt.colorbar()
plt.show()
