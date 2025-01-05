
import sys
import os

import numpy as np

# ../bin/AnalysisVSPAERO.py をモジュールとしてインポート
sys.path.append(os.path.join('..')) # 親ディレクトリをモジュール探索パスに追加
from bin.AnalysisVSPAERO import *
from bin.ISAspecification import *

import openvsp as vsp

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
    trimed_polar = make_CDo_correction(vsp, trimed_polar, Weight=Weight, xTr=(0.5, 0.7), CDpCL=0.0016, thickness=0.19, interference_factor=1.14, altitude=0, dT=0)
    trimed_polar.to_csv('G103A_DegenGeom_trimed.polar', index=False)
