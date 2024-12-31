
import sys
import os
import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from multiprocessing import Pool

# ../bin/AnalysisVSPAERO.py をモジュールとしてインポート
sys.path.append(os.path.join('..')) # 親ディレクトリをモジュール探索パスに追加
from bin.AnalysisVSPAERO import *
from bin.ISAspecification import *

import openvsp as vsp

def make_CDo_correction(vsp, trimed_polar, Weight, CDp0, CDpCL, altitude=0, dT=0):

    reynolds = trimed_polar['Re_1e6'].values * 1e6
    CDo = (0.455 / (np.log10(reynolds) ** 2.58) - 1700/ reynolds) * 2
    Sref = vsp.GetDoubleAnalysisInput('VSPAEROSweep', 'Sref')[0] # [m^2]
    density = get_density(altitude=altitude, dT=dT) # [kg/m^3]

    trimed_polar['CDo_corr'] = CDo + CDp0 + CDpCL * (trimed_polar['CL'].values) ** 2
    trimed_polar['CDtot_corr'] = trimed_polar['CDo_corr'].values + trimed_polar['CDi'].values
    trimed_polar['L_D_corr'] = trimed_polar['CL'].values / trimed_polar['CDtot_corr'].values

    trimed_polar['gamma'] = np.arctan(1/trimed_polar['L_D_corr'].values)
    trimed_polar['Velocity'] = np.sqrt((2*Weight*g)/(density*Sref*trimed_polar['CL'].values*np.cos(trimed_polar['gamma'].values)))
    trimed_polar['Vx'] = trimed_polar['Velocity'].values*np.cos(trimed_polar['gamma'])*(3.6)
    trimed_polar['Vz'] = trimed_polar['Velocity'].values*np.sin(trimed_polar['gamma'].values)
    
    return trimed_polar

if __name__=='__main__':
    
    # Clear the current VSP model and read the new VSP file
    vsp.ClearVSPModel()
    vsp.Update()
    vsp.ReadVSPFile('G103A.vsp3')
    vsp.Update()

    # Define the list of alpha angles and the Mach number for the sweep
    alpha_list = [-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 10, 12]
    mach = 0.1
    
    Weight = 580 # [kg]
    trimed_polar = vsp_trimed_sweep(vsp=vsp, alpha_list=alpha_list, Weight=Weight)
    trimed_polar = make_CDo_correction(vsp, trimed_polar, Weight=Weight, CDp0=0.0015, CDpCL=0.0018, altitude=0, dT=0)
    trimed_polar.to_csv('G103A_DegenGeom_trimed.polar', index=False)
