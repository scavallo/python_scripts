#!/home/scavallo/anaconda2/bin/python

import numpy as np

def mstats(m):
   eps = 2.2204e-16
   
   MatrixSize = np.size(m);
   MatrixShape = np.shape(m);
   try:
       nans = np.isnan(m);

       nans_vec = np.reshape(nans, np.prod(np.size(nans)), 1);   
       mm = m[~np.isnan(m)]
   except:
       m = np.array( m, dtype=float)
       nans = np.isnan(m);
       nans_vec = np.reshape(nans, np.prod(np.size(nans)), 1);   
       mm = m[~np.isnan(m)]      
   
   mp = np.reshape(m, np.prod(np.size(m)), 1);     
   Nnans = np.sum(nans_vec);
   
   NElements = np.prod(MatrixSize);
   NAnalyzedElements = NElements - Nnans; 

   Mean = np.mean(mm); 
   Max = np.max(mm[np.isfinite(mm)]);
   Min = np.min(mm[np.isfinite(mm)]);
   Range = Max-Min; 
   Median = np.median(mm[np.isfinite(mm)]);
   StDev = np.std(mm[np.isfinite(mm)]);
   absmp = np.abs(mm[np.isfinite(mm)]);
   MeanAbs = np.mean(absmp);      
   MinAbs = np.min(absmp[absmp>eps]);       
   FracZero = float(len(mp[mp==0]))/float(NAnalyzedElements)            
   FracNaN = float(Nnans)/float(NElements);
   
   print(" ");
   print("MatrixSize = ", MatrixSize);
   print("MatrixShape = ", MatrixShape);
   print("NElements = ", NElements);
   print("NAnalyzedElements = ", NAnalyzedElements); 
   print("Mean = ", Mean);
   print("Median = ", Median);
   print("Max = ", Max);
   print("Min = ", Min);
   print("Range = ", Range);   
   print("StDev = ", StDev);
   print("MeanAbs = ", MeanAbs);
   print("MinAbs = ", MinAbs);
   print("FracZero = ", FracZero);
   print("FracNaN = ", FracNaN);

   return
