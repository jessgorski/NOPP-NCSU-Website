# -*- coding: utf-8 -*-
"""
Created on Mon May  25 10:39:02 2022

@author: jfgorski

This code plots the 1D results for simulations. Includes: 
    - lat/lon
    - Dune Crest Elev change
    - Volume change (crest to shoreline)
    - "Sallenger scale" (collision vs overwash/inundation)
    
"""

import netCDF4 as nc
import matplotlib.pyplot as plt
import os
import time
import numpy as np
from mpl_toolkits.basemap import Basemap


def Dune_Crest_Elev_change(crest_index, T_pre, T_post):
    """This function takes two lists of elevation data,
    Finds the difference between the max elevation of pre and post DEM """

    init = T_pre[crest_index]
    final = T_post[crest_index]
    if final != None:
        diff = init - final
    else:
        diff = None
    return diff


def Vol_change(shore_index, crest_index, X, T_pre, T_post):
    """This function computes the volume change between the shoreline and the 
    dune crest of the pre-storm elevation """
    
    vol=0
    for a in range(shore_index,(crest_index)):
        dist_last = X[a] - X[a-1]
        dist_next = X[a+1] - X[a]
        X_dist= 0.5*(dist_last + dist_next)
        Y_diff = T_pre[a] - T_post[a]
        if abs(Y_diff) <=100:
            vol += (X_dist * Y_diff)    
    return vol 


def Regime(crest_ind, xb_out):
    """Returns either 0 = collision or 1 = inundation/overwash if there is 
    water at the crest"""       
    
    TIME = xb_out.variables['globaltime']
    
    WL = []
    count =0
    for t in range(len(TIME)):
        wl = xb_out.variables['zs'][t][0][crest_ind]
        if wl !=0:
            count+=1
        WL.append(wl)    
    if count >=1:
        regime = "Inundation/Overwash"
    else: 
        regime = "Collision"
    return regime


#------------------------------------------------------------------------------------------------------------------------------
t0 = time.time()
path = os.getcwd()
transects = os.listdir(path)
print (path, transects)

image_link = 'https://raw.github.ncsu.edu/jfgorski/NOPP-NCSU-XBeach/main/IRENE_10min_plots/'

try:
    os.makedirs(path+"/plots/")
except OSError:
    print ("/%s already exists, updating content" % path)
else:
    print ("Successfully created the directory /%s" % path)

lon_0 = -76.596
lat_0 = 34.664

max_lat = lat_0 + 2.0
min_lat = lat_0 - 2.0
max_lon = lon_0 + 5.0
min_lon = lon_0 - 5.0

plt.figure(figsize=(10,8))   
my_map = Basemap(projection='ortho', lon_0=lon_0, lat_0=lat_0, resolution='i')

#my_map.bluemarble()
#my_map.shadedrelief()
my_map.drawcoastlines()
my_map.drawcountries()

xmin, ymin = my_map(min_lon, min_lat)
xmax, ymax = my_map(max_lon, max_lat)

google_csv = []
google_csv.append(['T','lon','lat','Regime','Link'])

for u in transects:    
    try:
        T = int(u[1:])
        #get lat lon from info file
        g=[]
        info_file = path+'/T{}/T{}_info.txt'.format(T,T)
        with open(info_file,'r') as data:
            for line in data:
                o = line.strip().split()
                g.append([o[0],o[1]])
        lon = round(float(g[6][1]),3)
        lat = round(float(g[7][1]),3)
        crest_elev = float(g[9][1])
        
        file=path+"/T{}/xboutput.nc".format(T)
        xb_out = nc.Dataset(file,'r')
        X = xb_out.variables['globalx'][0]
        init = xb_out.variables['zb'][0]
        final = xb_out.variables['zb'][-1]
        
        X=abs(X)
        
        q=[]
        for r in init[0]:
            if r != None:
                q.append(r)
            else:
                q.append(-9999)

        crest_ind = q.index(crest_elev)
        zero_crossings = np.where(np.diff(np.sign(init[0])))[0]
        shore_index=zero_crossings[0]
        crest_x = X[crest_ind]-X[shore_index]
        dist_arr = np.array(X[shore_index:])
        shore_dist = dist_arr - X[shore_index]
        
        D = Dune_Crest_Elev_change(crest_ind,init[0],final[0])
        V = Vol_change(shore_index,crest_ind,X,init[0],final[0])
        S = Regime(crest_ind,xb_out)
    
        g=plt.figure()
        plt.plot(crest_x,init[0][crest_ind],'r*',label="Dune Crest")
        plt.plot(shore_dist,init[0][shore_index:],'grey',label='init')
        plt.plot(shore_dist,final[0][shore_index:],'b--',label='final')
        plt.xlabel("Distance from mean shoreline (meters)")
        plt.ylabel("Elevation (meters)")
        plt.title('Shoreline Lon/Lat: {} {}'.format(lon,lat),fontdict={'fontsize': 6})
        plt.suptitle('Dune Crest Elev Change = {}m\n Volume Change = {}m3\n Regime = {}'.format(round(D,2),round(V,3),S),fontsize=6)
        plt.legend()
        filen=path+'/plots/T{}_results.png'.format(T)
        g.savefig(filen)
        plt.close()
        
        loc = str(image_link)+'T{}_results.png'.format(T)
        google_csv.append([T,lon,lat,S,loc])

        # Map (long, lat) to (x, y) for plotting
        xs=[]
        ys=[]
        xy_file = path+'/T{}/T{}_XYdata.txt'.format(T,T)
        with open(xy_file,'r') as f:
            first = f.readline()
            f_line = first.strip().split()
            xs.append(float(f_line[0]))
            ys.append(float(f_line[1]))
            for line in f:
                n=line.strip().split()
            l_line = n
        xs.append(float(l_line[0]))
        ys.append(float(l_line[1]))
        x,y = my_map(xs, ys)
        if S=="Collision":
            my_map.plot(x, y, 'g-')
        elif S=="Inundation/Overwash":
            my_map.plot(x , y, 'r-')
        print("Added T{} to plot".format(T))
    except:
        print('No Transect {}'.format(T))


ax = plt.gca()         # set the axes limits
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)

output_file = path+'/plots/FinalMap.png'
plt.savefig( output_file, bbox_inches='tight', dpi=300)        # save image

plt.close('all')

with open(path+"/plots/Shoreline_locations_w_Regime.csv", "w") as f:
    for row in google_csv:
        f.write("%s\n" % ','.join(str(col) for col in row))

t1 = time.time()
print ("Time to make plots: {} s".format(round(t1-t0,3))) 

   