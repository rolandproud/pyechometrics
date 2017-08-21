# -*- coding: utf-8 -*-
"""
Summarise Sound Scattering Layers (SSLs)

@author: Roland Proud
"""

## import packages
import matplotlib.pyplot as plt
import gzip
import pickle
import numpy as np
from pyechoplot.plotting import plot_pseudo_SSL, save_png_plot, plot_Sv

## import pyechometrics modules
from pyechometrics.metrics import stats, dims, nasc

## get Sv data and mask
def get_obj(filepath):
    f   = gzip.open(filepath,'rb')
    obj = pickle.load(f,encoding = 'bytes')
    f.close()
    return obj

## noise_level
noise_level = -999

## read Sv
Sv18 = get_obj('./data/PS_Sv18.pklz')

## get SSL mask - see 'ident_SSLs' example in pyechomask
Sv18mask = get_obj('./data/SSL_flag_mask_18.pklz')

## plot
plt.figure(1)
plt.subplot(211)
plot_Sv(Sv18)
plt.subplot(212)
plot_Sv(Sv18,mask = Sv18mask)
plt.title('SSL identification - 18 kHz echosounder data')
plt.show()

## sample interval in meters for this echogram
sample_int  = 0.2 ## in meters

## calculate NASC (include all SSLs)
NASC = nasc(Sv18, sample_int, mask = Sv18mask)

## plot NASC by ping
plt.plot(NASC)
plt.xlabel('ping')
plt.ylabel(r'NASC $m^2nmi^{-2}$')
plt.title('NASC values for SSLs')
plt.show()

## save plot   
#save_png_plot('./','NASCexampleWiki')

## make binary mask for a single sound scattering layer (SSL) (Sv18mask == 2)
SSLmask                = np.zeros(Sv18mask.shape)
SSLmask[Sv18mask == 2] = 1

## get SSL stats and dimensions
SSL_mean, SSL_median, SSL_std, n             = stats(Sv18, mask = SSLmask)
mean_row, mean_height, mean_col, mean_length = dims(Sv18, mask = SSLmask)

## change row to depth
mean_depth  = mean_row * sample_int
mean_height = mean_height * sample_int

## plot a pseudo SSL using metrics
## *assume single normal distribution
plot_pseudo_SSL(SSL_mean,SSL_std,mean_height,mean_depth)
plt.ylabel('depth (m)')  
plt.xlabel('pings')
plt.title('pseudo DSL produced using summary metrics',fontsize = 16)
plt.show()

## save plot   
#save_png_plot('./','exampleWiki')






