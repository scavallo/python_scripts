#!/usr/bin/python

import numpy as np

def pstats(m):
    """
    
    p_Mean, p_Max, p_Min, p_Median, p_StDev = pstats(m)
    
    Computes the basic statistics of input array m
    
    Steven Cavallo
    University of Oklahoma
    March 2014
        
    """
    eps = 2.2204e-16

    MatrixSize = np.size(m);
    MatrixShape = np.shape(m);
    try:
        nans = np.isnan(m);
    except:
        m = np.array( m, dtype=float)
        nans = np.isnan(m);    
    
    nans_vec = np.reshape(nans, np.prod(np.size(nans)), 1);   
    mm = m[~np.isnan(m)]

    mp = np.reshape(m, np.prod(np.size(m)), 1);     
    Nnans = np.sum(nans_vec);

    NElements = np.prod(MatrixSize);
    NAnalyzedElements = NElements - Nnans; 

    p_Mean = np.mean(mm); 
    p_Max = np.max(mm[np.isfinite(mm)]);
    p_Min = np.min(mm[np.isfinite(mm)]);
    p_Range = p_Max-p_Min; 
    p_Median = np.median(mm[np.isfinite(mm)]);
    p_StDev = np.std(mm[np.isfinite(mm)]);
    p_absmp = np.abs(mm[np.isfinite(mm)]);
    p_MeanAbs = np.mean(p_absmp);      
    p_MinAbs = np.min(p_absmp[p_absmp>eps]);       
    p_FracZero = float(len(mp[mp==0]))/float(NAnalyzedElements)            
    p_FracNaN = float(Nnans)/float(NElements);

    print(" ");
    print("MatrixSize = ", MatrixSize);
    print("MatrixShape = ", MatrixShape);
    print("NElements = ", NElements);
    print("NAnalyzedElements = ", NAnalyzedElements); 
    print("Mean = ", p_Mean);
    print("Median = ", p_Median);
    print("Max = ", p_Max);
    print("Min = ", p_Min);
    print("Range = ", p_Range);   
    print("StDev = ", p_StDev);
    print("MeanAbs = ", p_MeanAbs);
    print("MinAbs = ", p_MinAbs);
    print("FracZero = ", p_FracZero);
    print("FracNaN = ", p_FracNaN);

    return p_Mean, p_Max, p_Min, p_Median, p_StDev
