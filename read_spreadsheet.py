# -*- coding: utf-8 -*-
"""
Created on Sat Oct 15 22:40:48 2016

@author: lj
"""
import numpy as np
import matplotlib.dates as md
import matplotlib.pylab as plt
import datetime,dateutil
def read_archive_spread(filename):
    infile=open(filename,'r')
#    print infile.readline()
    idnames=(infile.readline()).split('\t')
    print idnames
    dates=[]
    paras=[]
    dict1={}
    for line in infile:
        columns=line.split('\t')
        date = columns[0]
        para = columns[1:]
        para = [float(x) for x in para]
        dates.append(date)
        paras.append(para)
    infile.close()
    dict1[idnames[0]]=dates
    paras=np.array(paras)
    for i in range(len(idnames)):
        if i>0:
            dict1[idnames[i]]=paras[:,i-1]        
    return dict1
tz=read_archive_spread('CGIExport2.txt')
print tz
xxxx=tz.keys()
print xxxx
datestr=tz[xxxx[6]]
datex=[dateutil.parser.parse(s) for s in datestr]
print datex
fig,ax=plt.subplots()
datesFmt = md.DateFormatter('%m/%d/%Y %H:%M:%S.%f')
ax.xaxis.set_major_formatter(datesFmt)
ax.plot_date(datex,tz[xxxx[1]])
plt.xticks( rotation=25 )


plt.show()
#datecol= [datetime.datetime.strptime(elem, '%m/%d/%Y %H:%M:%S.%f') for elem in tz[xxxx[6]] ]
#(fig,ax)=plt.subplots(1,1)
#ax.plot(datecol,tz[xxxx[1]])
#plt.show()
#timed,raw=np.loadtxt('CGIExport2.txt',skiprows=1,delimiter='\t',unpack=True,dtype=[('date', '|S10'), ('floatmio', float)])#converters={ 0: md.strpdate2num('%m/%d/%Y %H:%M:%S.%f')})


#convertfunc = lambda x: dateutil.parser.parse(x)
#timed=np.loadtxt('CGIExport2.txt',skiprows=1,usecols=[0],delimiter='\t')
