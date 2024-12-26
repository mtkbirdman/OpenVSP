import sys
import os

import numpy as np

# ../bin/AnalysisVSPAERO.py をモジュールとしてインポート
sys.path.append(os.path.join('..')) # 親ディレクトリをモジュール探索パスに追加
from bin.AnalysisVSPAERO import *

import openvsp as vsp

if __name__=='__main__':

    # Close and open the file
    vsp.ClearVSPModel()
    vsp.Update()
    vsp.ReadVSPFile('G103A.vsp3')
    vsp.Update()

    # Set control surface
    set_control_surface(geom_name='WingGeom',  deflection=10,  cs_group_name='AILERON_GROUP', gains=(1,1))
    set_control_surface(geom_name='HTailGeom', deflection=0, cs_group_name='Elevator_Group', gains=(1,-1))
    set_control_surface(geom_name='VTailGeom', deflection=0,  cs_group_name='Rudder_Group')
    vsp.Update()

    # Execute VSPAEROSweep
    alpha = np.linspace(-4, 12, 9)
    mach = [0.1]
    vsp_sweep(vsp=vsp, alpha=alpha, mach=mach)
