# -*- coding: utf-8 -*-
"""
.. :module:: metrics
    :synopsis: summarise signal

| Developed by: Roland Proud (RP) <rp43@st-andrews.ac.uk> 
|               Pelagic Ecology Research Group, University of St Andrews
| Contributors:
|
| Maintained by:
| Modification History:      
|
"""

import numpy as np

def abc(Sv,sample_int, mask = None,noise_level = -999):
    '''
    :param Sv: gridded Sv values (dB re 1m^-1)
    :type  Sv: 2D numpy.array
    
    :param sample_int: sample interval in meters
    :type  sample_int: float
    
    :param mask: binary mask (0 - noise; 1 - signal)
    :type  mask: 2D numpy.array
    
    :param noise_level: level of background noise (db re 1m^-1)
    :type  noise_level: float
    
    
    return:
    :param ABC: ABC values (m^2m^-2)
    :type  ABC: 1D numpy.array
    
    desc: Calculate Area backscattering coefficient (ABC) values
    
    defined by RP
    
    status: dev
    
    
    '''
    Sv = 10**(Sv/10.)  ##  log to linear
    if mask is not None:
        Sv[mask == 0] = 10**(noise_level/10.)
    Sv  = np.ma.masked_invalid(Sv)
    ABC = np.ma.mean(Sv,axis = 0) * np.ma.count(Sv,axis = 0) * sample_int

    return ABC

def nasc(Sv,sample_int, mask = None,noise_level = -999):
    '''
    :param Sv: gridded Sv values (dB re 1m^-1)
    :type  Sv: 2D numpy.array
    
    :param sample_int: sample interval in meters
    :type  sample_int: float
    
    :param mask: binary mask (0 - noise; 1 - signal)
    :type  mask: 2D numpy.array
    
    :param noise_level: level of background noise (db re 1m^-1)
    :type  noise_level: float
    
    
    return:
    :param NASC: NASC values (m^2nmi^-2)
    :type  NASC: 1D numpy.array
    
    desc: Calculate Nautical area scattering coefficient (NASC) values
    
    defined by RP
    
    status: dev
    
    '''
    return abc(Sv,sample_int, mask = mask, noise_level = noise_level)\
                 * 1852**2 * 4 * np.pi

def idx_feature(Sv, mask = None,noise_level = -999,linear = False):
    '''
    :param Sv: gridded Sv values (dB re 1m^-1)
    :type  Sv: 2D numpy.array
        
    :param mask: binary mask (0 - noise; 1 - signal)
    :type  mask: 2D numpy.array
    
    :param noise_level: level of background noise (db re 1m^-1)
    :type  noise_level: float
    
    :param linear: return linear of decibel values (True = linear)
    :type  linear: boolean
    
    
    return:
    :param Sv: gridded masked Sv values (dB re 1m^-1)
    :type  Sv: 2D numpy.masked.array
    
    :param row_idx: row index of signal values
    :type  row_idx: 1D numpy.array
    
    :param col_idx: column index of signal values
    :type  col_idx: 2D numpy.array
    
    desc: mask Sv, change to linear if required and index data values
    
    defined by RP
    
    status: dev
    
    '''
    ## linear or log?
    if linear:
        Sv          = 10**(Sv/10.)
        noise_level = 10**(noise_level/10.)
        
    ## mask noise and invalid
    Sv = np.ma.masked_where(Sv <= noise_level,Sv)
    Sv = np.ma.masked_invalid(Sv)
    
    ## select mask
    if mask is not None:
        row_idx,col_idx  = np.where(mask == 1)
    else:
        row_idx,col_idx  = np.where(Sv.mask == False)
        
    return Sv,row_idx,col_idx
        
def stats(Sv, mask = None,noise_level = -999):
    '''
    :param Sv: gridded Sv values (dB re 1m^-1)
    :type  Sv: 2D numpy.array
    
    :param mask: binary mask (0 - noise; 1 - signal)
    :type  mask: 2D numpy.array
    
    :param noise_level: level of background noise (db re 1m^-1)
    :type  noise_level: float
    
    return:
    
    desc:  return common statistics for a masked feature
    
    defined by RP
    
    status: dev
    
    '''
    ## mask noise and get data idx
    Sv,row_idx,col_idx = idx_feature(Sv, mask = mask,noise_level = noise_level,\
                                     linear = True)
    ## get stats
    Sv_median         = np.ma.median(Sv[row_idx,col_idx])
    Sv_mean           = np.ma.mean(Sv[row_idx,col_idx])
    Sv_std            = np.ma.std(Sv[row_idx,col_idx])
    n                 = len(row_idx)
    
    return Sv_mean, Sv_median, Sv_std, n

def dims(Sv, mask = None,noise_level = -999):
    ''' 
    :param Sv: gridded Sv values (dB re 1m^-1)
    :type  Sv: 2D numpy.array
    
    :param mask: binary mask (0 - noise; 1 - signal)
    :type  mask: 2D numpy.array
    
    :param noise_level: level of background noise (db re 1m^-1)
    :type  noise_level: float
    
    return:
    
    desc:  return row/col dimensions for a masked feature
    
    defined by RP
    
    status: dev
    '''
    ## mask noise and get data idx
    Sv,row_idx,col_idx = idx_feature(Sv, mask = mask,noise_level = noise_level)
        
    ## get dims
    mean_row         = np.mean(row_idx)
    mean_height      = len(row_idx)/len(np.unique(col_idx))  
    mean_col         = np.mean(col_idx)
    mean_length      = len(col_idx)/len(np.unique(row_idx)) 


    return mean_row, mean_height, mean_col, mean_length







