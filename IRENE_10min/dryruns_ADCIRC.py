# -*- coding: utf-8 -*-
"""
Created on Thu May  5 12:37:01 2022

@author: jfgorski
"""
import numpy as np
from netCDF4 import Dataset
from scipy import spatial
from operator import itemgetter
import shutil
import os
from geopy import distance
import time


def GetNodes(fort_file,stat_list):
    """Finds ADCIRC nodes, use Alireza's kdtree method """
    ncf=Dataset(fort_file,'r')
    elem  = ncf.variables['element'][:]
    xnode = ncf.variables['x'][:]
    ynode = ncf.variables['y'][:]
    nodes=np.column_stack((xnode,ynode))


    ### nearest k nodes to each station in stat_list 
    CKD=list(spatial.cKDTree(nodes).query(stat_list,k=4)[1])

    triangles=[]
    for i in range(len(stat_list)): #for each station
        tri=[]    
        ps=CKD[i]                #nearest k nodes 
        for node in ps:
            tri.extend(list(np.where(elem == node+1)[0]))  # find nearest elements/triangles
        triangles.append(set(tri))    

    GAMMA=[]
    ADJ=[]
    for i in range(len(stat_list)):
        xp,yp=stat_list[i]
        for j in triangles[i]:
            adjc=elem[j]-1
            XX  = [xp,xnode[adjc[0]],xnode[adjc[1]],xnode[adjc[2]]]
            YY  = [yp,ynode[adjc[0]],ynode[adjc[1]],ynode[adjc[2]]]

            def Area(X,Y):
                return round(abs((X[1]*Y[2]-X[2]*Y[1])-(X[0]*Y[2]-X[2]*Y[0])+(X[0]*Y[1]-X[1]*Y[0])),12)
            def Sub_Area(a,b,c):
                return Area(list(itemgetter(a,b,c)(XX)),list(itemgetter(a,b,c)(YY)))

            Tot_A  = Area(XX[1:],YY[1:])
            G0 = Sub_Area(0,2,3)/Tot_A
            G1 = Sub_Area(1,0,3)/Tot_A
            G2 = Sub_Area(1,2,0)/Tot_A
            
            if round(G0+G1+G2,5) == 1.0 :
                print(i,[G0,G1,G2])
                GAMMA.append([G0,G1,G2])          
                ADJ.append(adjc)
                break
                 
    return(GAMMA,ADJ)
    
def PullT(file,dt,GAMMA,ADJ,stat_list):
    """use GAMMA and ADJ to get data"""

    ncf=Dataset(file,'r')
    var = ncf.variables[dt][:]
    nt=var.shape[0]
              
    Stat_val=[]
    ## output data
    for s in range(len(stat_list)):
        G0,G1,G2=GAMMA[s]
        n0,n1,n2=ADJ[s]
        S=[]
        for t in range(nt):
            mval= G0*var[t][n0]+G1*var[t][n1]+G2*var[t][n2]
            try:
                mval= round(mval,4)
                S.append(float(mval))
            except:
                mval = 0.01
                S.append(float(mval))
        Stat_val.append(S)
    return (Stat_val)              
    
def params_file(Transect,length,fi_deg,tstop):
    """creates param file for exporting"""
    paramfile=open(path+'params.txt','w')
    
    nx=length-1
    
    thetamin=180-fi_deg
    thetamax=360-fi_deg

    
    textblock='''%%%%%%   This is for Transect {}
%%%% GRID %%%%
gridform     = xbeach
depfile      = y_file_tnsct_{}.txt
posdwn       = -1
nx           = {}
ny           = 0
vardx        = 1
alfa         = {}
xfile        = x_file_tnsct_{}.txt
thetamin     = {}
thetamax     = {}
dtheta       = 20
thetanaut    = 1


%%%% WAVE %%%%
wbcversion= 3
instat    = jons_table
bcfile    = wv_{}.txt
random    = 0
leftwave  = neumann
rightwave = neumann
back = wall

%%%% WATER LVL %%%%
tideloc = 1
zs0file = wl_{}.txt
zs0=0
left       = neumann_v
right      = neumann_v


%%%% PAPER PARAM %%%%
break=roelvink1
gamma=0.42
alpha=2
beta=0.1
facAs=0.1
hmin=0.2
eps=0.01
bdslpeffmag=roelvink_bed
lws=0
fw=0.02
wetslp=0.3

%%%%%no need to chng
ngd=1
nd=3
D50=0.0003
D90=0.0005
morfac = 10

%%%% Output %%%
morstart= 3600
tstart = 0
tintg  = 1800
tstop  = {}

outputformat = netcdf
nglobalvar = 3
zs
zb
H

'''.format(Transect,Transect, nx, fi_deg, Transect, thetamin, thetamax, Transect, Transect, tstop)  

    paramfile.write(textblock)
    paramfile.close()


#                            USER DEFINED PATHS AND VARS
#______________________________________________________________________________

wv_file= 'swan_HS.63.nc'
tp_file= 'swan_TPS.63.nc'
dir_file= 'swan_DIR.63.nc'
fort_file= 'fort.63.nc'
   
off_station_file = '../Final_profiles/offshore_stations.txt'

###############################################################################
# ask for user input for landfall location and radius
t0=time.time()

landfall_lon = input("please eneter the longitutde of the storm landfall (ex: -76.596):")
landfall_lat = input("please enter the latitude of the storm landfall (ex: 34.664): ")
radius = input("please eneter the radius of interest (kilometers):")

#                   END TRANSECT GENERATION, START WV & WL
#_______________________________________________________________________________
stat_list=[]
sim_list=[]
with open(off_station_file,'r') as f_st:
    f_st.readline()
    for line in f_st:
        m=line.split()
        # calc distance between offshore location and point
        landfall = (float(landfall_lat), float(landfall_lon))
        new = (float(m[2]),float(m[1]))
        dist = distance.distance(landfall,new).km
        if dist <= float(radius):
            sim_list.append(int(m[0]))
            stat_list.append([float(m[1]),float(m[2])])
            
gamma,adj=GetNodes(fort_file,stat_list)

t1 = time.time()
print("Time to find nodes: {} s".format(round(t1-t0,3)))

HS=PullT(wv_file,str('swan_HS'),gamma,adj,stat_list)
TPS=PullT(tp_file,str('swan_TPS'),gamma,adj,stat_list)
DIR=PullT(dir_file,str('swan_DIR'),gamma,adj,stat_list)
WL=PullT(fort_file,str('zeta'),gamma,adj,stat_list)   

t2 = time.time()
print("Time to PullT: {} s".format(round(t2-t1,3)))

G=0
for T in sim_list:
    source_dir = "/gpfs_common/share01/jcdietri/jfgorski/grad_project/Final_profiles/T{}/".format(T)
    dirr = os.getcwd()
    path = dirr+"/T{}/".format(T)
    try:
        shutil.copytree(source_dir, path)
    except:
        shutil.rmtree(path)
        shutil.copytree(source_dir, path)
    
    filename_wv=path+str('wv_{}.txt'.format(T))
    filename_wl=path+str('wl_{}.txt'.format(T)) 
    
    # every 10 minutes 
    e = 0
    while HS[G][e]<=0.01:
        e+=1
    
    with open(filename_wv,'w') as wv:
        for i in range(e,len(HS[G])):
            wvline='{} {} {} {} {} {} {}'.format(HS[G][i],TPS[G][i],(270-DIR[G][i]),3.3,20,600,1)
            wv.write(str(wvline))
            wv.write('\n')
    k=0
    with open(filename_wl,'w') as wl:
        for j in range(e,len(HS[G])):
            wlline='{} {}'.format(k*600,WL[G][j])
            wl.write(str(wlline))
            wl.write('\n')
            k+=1
            
    last_t=(k-1)*600
    #get nx and fi from info file
    g=[]
    info_file = path+'T{}_info.txt'.format(T)
    with open(info_file,'r') as data:
        for line in data:
            o = line.strip().split()
            g.append([o[0],o[1]])
            
    nx=int(g[2][1])
    fi_deg=float(g[1][1])
    params_file(T,nx,fi_deg,last_t)
    G+=1
    
t3 = time.time()
print("Time to prepare all files: {} s".format(round(t3-t0,3)))