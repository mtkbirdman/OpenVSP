# https://openvsp.org/pyapi_docs/latest/openvsp.html
# https://github.com/OpenVSP/OpenVSP/blob/fc4a97b75c92399fecf74b21f2a57087f89e12b7/examples/scripts/TestAnalysisVSPAERO.vspscript

import os
import sys
from contextlib import contextmanager

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

import openvsp as vsp

from .ISAspecification import *

@contextmanager
def suppress_stdout():

    # 標準出力を一時的に/dev/nullにリダイレクト
    with open(os.devnull, 'w') as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout

def vsp_sweep(vsp, alpha, mach, reynolds=[1e6], verbose=1):

    if verbose:
        print('\n-> Calculate alpha & mach sweep analysis\n')

    # //==== Analysis: VSPAero Compute Geometry to Create Vortex Lattice DegenGeom File ====//
    
    # Set defaults
    compgeom_name = 'VSPAEROComputeGeometry'
    vsp.SetAnalysisInputDefaults(compgeom_name)

    # List inputs, type, and current values
    if verbose:
        print(compgeom_name)
        vsp.PrintAnalysisInputs(compgeom_name)
        print('')

    # Execute
    if verbose:
        print('\tExecuting...')
    compgeom_resid = vsp.ExecAnalysis(compgeom_name)
    if verbose:
        print('\tCOMPLETE')

    # Get & Display Results
    if verbose:
        vsp.PrintResults(compgeom_resid)
        print('')

    # //==== Analysis: VSPAero Sweep ====//

    # Set defaults
    analysis_name = 'VSPAEROSweep'
    vsp.SetAnalysisInputDefaults(analysis_name)

    # Reference geometry set
    geom_set = [0]
    vsp.SetIntAnalysisInput(analysis_name, 'GeomSet', geom_set, 0)
    ref_flag = [1]
    vsp.SetIntAnalysisInput(analysis_name, 'RefFlag', ref_flag, 0)
    wid = vsp.FindGeomsWithName('WingGeom')
    vsp.SetStringAnalysisInput(analysis_name, 'WingID', wid, 0)

    # Freestream Parameters
    alpha_npts = [len(alpha)]
    vsp.SetDoubleAnalysisInput(analysis_name, "AlphaStart", [alpha[0]], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "AlphaEnd", [alpha[-1]], 0)
    vsp.SetIntAnalysisInput(analysis_name, "AlphaNpts", alpha_npts, 0)
    mach_npts = [len(mach)]
    vsp.SetDoubleAnalysisInput(analysis_name, 'MachStart', [mach[0]], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, 'MachEnd', [mach[-1]], 0)
    vsp.SetIntAnalysisInput(analysis_name, 'MachNpts', mach_npts, 0)
    reynolds_npts = [len(reynolds)]
    vsp.SetDoubleAnalysisInput(analysis_name, 'ReCref', [reynolds[0]], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, 'ReCrefEnd', [reynolds[-1]], 0)
    vsp.SetIntAnalysisInput(analysis_name, 'ReCrefNpts', reynolds_npts, 0)
    vsp.Update()

    # List inputs, type, and current values
    if verbose:
        print(analysis_name)
        vsp.PrintAnalysisInputs(analysis_name)
        print('')

    # Execute
    if verbose:
        print('\tExecuting...')
        rid = vsp.ExecAnalysis(analysis_name)
        print('\tCOMPLETE')
    else:
        with suppress_stdout():
            rid = vsp.ExecAnalysis('VSPAEROSweep')

    # Get & Display Results
    # vsp.PrintResults(rid)

    return rid

def MyFindParm(cs_group_container_id, parm_name):

    # Find param_id from cs_group_container_id and parm_name
    parm_ids = vsp.FindContainerParmIDs( cs_group_container_id )
    for parm_id in parm_ids:
        if parm_name == vsp.GetParmName(parm_id):
            return parm_id
    return ''

def set_control_surface(geom_name, deflection, cs_group_name, sub_id=0, gains=(1,1), verbose=1):

    # Set control surface group
    current_group_names = [vsp.GetVSPAEROControlGroupName(group_index) for group_index in range(vsp.GetNumControlSurfaceGroups())]
    if cs_group_name in current_group_names:
        group_index = current_group_names.index(cs_group_name)
    else:
        group_index = vsp.CreateVSPAEROControlSurfaceGroup()
        vsp.SetVSPAEROControlGroupName(cs_group_name, group_index)

    # Add subsurface to control surface group
    cs_name_vec = vsp.GetAvailableCSNameVec(group_index)
    selcted = [1 + i for i, cs_name in enumerate(cs_name_vec) if geom_name in cs_name]
    vsp.AddSelectedToCSGroup(selcted, group_index)
    
    # Set gain (deferential)
    if verbose:
        print('\n', cs_group_name, ':', vsp.GetActiveCSNameVec(group_index))
        print('\t', f"{'container_id':14s}", f"{'parm_name':24s}", f"{'group_name':22s}", f"{'parm_id':14s}", f"{'value':8s}")
    geom_id = vsp.FindGeomsWithName(geom_name)[0]
    group_name = 'ControlSurfaceGroup_' + str(group_index)
    cs_group_container_id = vsp.FindContainer('VSPAEROSettings', 0)    
    for i, gain in enumerate(gains):
        parm_name = 'Surf_' + vsp.GetSubSurf(geom_id, sub_id) + '_' + str(i) + '_Gain'
        parm_id = MyFindParm(cs_group_container_id, parm_name)
        # parm_id = vsp.FindParm(cs_group_container_id, parm_name, group_name) # なぜか反応しない...
        if parm_id != '':
            vsp.SetParmVal(parm_id, gain)
            if verbose:
                print('\t', f'{cs_group_container_id:14s}', f'{parm_name:24s}', f'{group_name:22s}', f'{parm_id:14s}', f'{gain:8.3f}')
        else:
            if verbose:
                print('\t', 'Failed to find ' + parm_name)

    # Set defrection angle
    parm_name = 'DeflectionAngle'
    parm_id = vsp.FindParm(cs_group_container_id, parm_name, group_name)
    vsp.SetParmVal(parm_id, deflection)
    vsp.Update()

    # Check if deflection angle successfully applied
    if deflection == vsp.GetParmVal(parm_id):
        if verbose:
            print('\t', f'{cs_group_container_id:14s}', f'{parm_name:24s}', f'{group_name:22s}', f'{parm_id:14s}', f'{deflection:8.3f}')
    else:
        if verbose:
            print('Failed to set deflection angle')
        exit()

    return parm_id

def get_polar_result(result_ids):

    # Loop through the results associated with the given result_ids
    for result_id in vsp.GetStringResults(result_ids, 'ResultsVec'):
        # Check if the result name is 'VSPAERO_Polar'
        if vsp.GetResultsName(result_id) == 'VSPAERO_Polar':
            data, columns = [], []
            # Loop through all available data names for this result_id
            for data_name in vsp.GetAllDataNames(result_id):
                # Get double results for the data_name and append if available
                double_results = vsp.GetDoubleResults(result_id, data_name, 0)
                if double_results:
                    data.append(double_results)
                    columns.append(data_name)
            # If data was found, return it as a pandas DataFrame
            if data:
                return pd.DataFrame(np.array(data).T, columns=columns)
    
    # Return an empty DataFrame if no data is found
    return pd.DataFrame()

def _get_trim(x, vsp, alpha, mach, reynolds, CMy_tol=5e-4):

    # Set the elevator deflection and control surface parameters
    _ = set_control_surface(geom_name='HTailGeom', deflection=x, cs_group_name='Elevator_Group', gains=(1,-1), verbose=0)
    
    # Perform a sweep analysis with the given alpha and mach
    result_ids = vsp_sweep(vsp=vsp, alpha=[alpha], mach=[mach], reynolds=[reynolds], verbose=0)
    
    # Retrieve polar results from the sweep analysis
    df = get_polar_result(result_ids)
    
    # Get the CMy (moment coefficient about the y-axis) value
    if list(df.columns):

        CMy = df['CMy'].values[0]
        f = np.abs(CMy)
    
        # Print the current deflection (x) and CMy value
        print(f'\rx = {x:8.4f} ', f'CMy = {CMy:10.6f}', end='')
        
        # If the absolute value of CMy is smaller than the tolerance, raise StopIteration to exit
        if f < CMy_tol:  
            raise StopIteration(x)
        
        # Return the absolute value of CMy to be minimized
        return f
    
    else:
        raise StopIteration(999)

def vsp_trimed_sweep(vsp, alpha_list, Weight, altitude=0, dT=0):

    g = get_gravity() # [m/s^2]
    density = get_density(altitude=altitude, dT=dT) # [kg/m^3]
    Sref = vsp.GetDoubleAnalysisInput('VSPAEROSweep', 'Sref')[0] # [m^2]
    cref = vsp.GetDoubleAnalysisInput('VSPAEROSweep', 'cref')[0] # [m]

    # Create an empty DataFrame to store trimmed polar results
    trimed_polar = pd.DataFrame()
    
    # Iterate over the alpha angles
    for alpha in alpha_list:

        result_ids = vsp_sweep(vsp=vsp, alpha=[alpha], mach=[0], verbose=0)
        df = get_polar_result(result_ids)
        velocity = np.sqrt((2*Weight*g)/(density*Sref*df['CL'].values[0]))
        mach = velocity_to_mach(velocity=velocity, dT=0, altitude=0)
        reynolds = velocity_to_reynolds(velocity=mach_to_velocity(mach, altitude=altitude, dT=dT), length=cref, altitude=altitude, dT=dT)

        try:
            # Minimize the absolute CMy to find the optimal elevator deflection (x) for the given alpha
            res = minimize_scalar(fun=_get_trim, args=(vsp, alpha, mach, reynolds), method='bounded', bounds=(-20, 20))
            de = res.x
        except StopIteration as err:
            # If StopIteration is raised, use the optimal deflection (err.value) and print 'Done'
            if err.value == 999:
                print('\tError\t', 'Alpha 'f'{alpha:8.3f}, velocity {velocity:8.3f}, mach {mach:8.3f}, reynolds {reynolds:8.3e}')
                de = None
            else:
                print('\tDone\t', 'Alpha 'f'{alpha:8.3f}, velocity {velocity:8.3f}, mach {mach:8.3f}, reynolds {reynolds:8.3e}')
                de = err.value
        
        if de:
            # Set the elevator deflection with the optimized value
            _ = set_control_surface(geom_name='HTailGeom', deflection=de, cs_group_name='Elevator_Group', gains=(1,-1), verbose=0)
            
            # Perform another sweep with the optimized elevator deflection
            result_ids = vsp_sweep(vsp=vsp, alpha=[alpha], mach=[mach], reynolds=[reynolds], verbose=0)
            
            # Retrieve polar results and add the deflection value to the results DataFrame
            df = get_polar_result(result_ids)
            df['de'] = de
            
            # Append the results to the trimmed polar DataFrame
            trimed_polar = pd.concat([trimed_polar, df])
    
    trimed_polar['gamma'] = np.arctan(1/trimed_polar['L_D'])
    trimed_polar['Velocity'] = np.sqrt((2*Weight*g)/(density*Sref*trimed_polar['CL']*np.cos(trimed_polar['gamma'])))
    trimed_polar['Vx'] = trimed_polar['Velocity']*np.cos(trimed_polar['gamma'])*(3.6)
    trimed_polar['Vz'] = trimed_polar['Velocity']*np.sin(trimed_polar['gamma'])

    return trimed_polar

def make_CDo_correction(vsp, trimed_polar, Weight, CDpCL=0.0065, thickness=0.12, interference_factor=1, xTr=(0,0), altitude=0, dT=0):

    # Function to calculate the frictional drag coefficient in turbulent regions
    def Cf_turbulance(reynolds):
        return 0.455/(np.log10(reynolds)**2.58)

    # Function to calculate the frictional drag coefficient in the laminar flow regime
    def Cf_laminar(reynolds):
        return 1.32824 / np.sqrt(reynolds)
    
    g = get_gravity()
    Sref = vsp.GetDoubleAnalysisInput('VSPAEROSweep', 'Sref')[0] # Get wing area [m^2]
    density = get_density(altitude=altitude, dT=dT) # Get density based on altitude [kg/m^3]
    reynolds = trimed_polar['Re_1e6'].values * 1e6 # Get Reynolds number

    # Correct drag based on percentage of laminar flow area
    CDo = 0
    for laminar_persent in xTr:
        reynolds_laminar = np.maximum(reynolds * laminar_persent, 1e3)  # Reynolds number in laminar flow region
        CDo += Cf_turbulance(reynolds) - Cf_turbulance(reynolds_laminar) * laminar_persent + Cf_laminar(reynolds_laminar) * laminar_persent

    # form factor
    form_factor = 1 + 2 * thickness + 60 * (thickness ** 4)

    # Corrected CD0, CDtotal, and lift-drag ratio are calculated and added to the data frame
    trimed_polar['CDo_corr'] = CDo * form_factor * interference_factor + CDpCL * (trimed_polar['CL'].values) ** 2
    trimed_polar['CDtot_corr'] = trimed_polar['CDo_corr'].values + trimed_polar['CDi'].values
    trimed_polar['L_D_corr'] = trimed_polar['CL'].values / trimed_polar['CDtot_corr'].values

    # Calculate angle of attack from modified lift-drag ratio
    trimed_polar['gamma'] = np.arctan(1/trimed_polar['L_D_corr'].values)

    # Calculate velocity
    trimed_polar['Velocity'] = np.sqrt((2*Weight*g)/(density*Sref*trimed_polar['CL'].values*np.cos(trimed_polar['gamma'].values)))

    # Calculates horizontal velocity (Vx) and vertical velocity (Vz) and adds them to the data frame
    trimed_polar['Vx'] = trimed_polar['Velocity'].values*np.cos(trimed_polar['gamma'])*(3.6)  # horizontal velocity [km/h]
    trimed_polar['Vz'] = trimed_polar['Velocity'].values*np.sin(trimed_polar['gamma'].values)  # vertical velocity [m/s]

    return trimed_polar  # Returns corrected data

def vsp_sweep_wig(vsp, alpha, mach, reynolds, height, AnalysisMethod=0, verbose=1):

    if verbose:
        print('\n-> Calculate alpha & mach sweep analysis\n')

    # //==== Analysis: VSPAero Compute Geometry to Create Vortex Lattice DegenGeom File ====//
    
    # Set defaults
    compgeom_name = 'VSPAEROComputeGeometry'
    vsp.SetAnalysisInputDefaults(compgeom_name)
    if AnalysisMethod:
        vsp.SetIntAnalysisInput(compgeom_name, 'AnalysisMethod', [AnalysisMethod], 0)

    # List inputs, type, and current values
    if verbose:
        print(compgeom_name)
        vsp.PrintAnalysisInputs(compgeom_name)
        print('')

    # Execute
    if verbose:
        print('\tExecuting...')
    compgeom_resid = vsp.ExecAnalysis(compgeom_name)
    if verbose:
        print('\tCOMPLETE')

    # Get & Display Results
    if verbose:
        vsp.PrintResults(compgeom_resid)
        print('')

    # //==== Analysis: VSPAero Sweep ====//

    # Set defaults
    analysis_name = 'VSPAEROSweep'
    vsp.SetAnalysisInputDefaults(analysis_name)

    # Reference geometry set
    geom_set = [0]
    vsp.SetIntAnalysisInput(analysis_name, 'GeomSet', geom_set, 0)
    ref_flag = [1]
    vsp.SetIntAnalysisInput(analysis_name, 'RefFlag', ref_flag, 0)
    wid = vsp.FindGeomsWithName('WingGeom')
    vsp.SetStringAnalysisInput(analysis_name, 'WingID', wid, 0)

    df = pd.DataFrame()
    for re in reynolds:
        for ma in mach:
            for he in height:
                for al in alpha:
                    # Freestream Parameters
                    vsp.SetDoubleAnalysisInput(analysis_name, "AlphaStart", [al], 0)
                    vsp.SetIntAnalysisInput(analysis_name, "AlphaNpts", [1], 0)
                    vsp.SetDoubleAnalysisInput(analysis_name, 'MachStart', [ma], 0)
                    vsp.SetIntAnalysisInput(analysis_name, 'MachNpts', [1], 0)
                    vsp.SetDoubleAnalysisInput(analysis_name, 'ReCref', [re], 0)
                    vsp.SetIntAnalysisInput(analysis_name, 'ReCrefNpts', [1], 0)
                    vsp.Update()

                    # Ground Effect
                    vsp.SetIntAnalysisInput(analysis_name, 'GroundEffectToggle', [1], 0)
                    vsp.SetDoubleAnalysisInput(analysis_name, 'GroundEffect', [he], 0)

                    # List inputs, type, and current values
                    if verbose:
                        print(analysis_name)
                        vsp.PrintAnalysisInputs(analysis_name)
                        print('')

                    # Execute
                    if verbose:
                        print('\tExecuting...')
                        rid = vsp.ExecAnalysis(analysis_name)
                        print('\tCOMPLETE')
                    else:
                        with suppress_stdout():
                            rid = vsp.ExecAnalysis('VSPAEROSweep')
                    
                    tmp = get_polar_result(rid)
                    bref = vsp.GetDoubleAnalysisInput('VSPAEROSweep', 'bref')[0] # [m]
                    tmp['Alpha'], tmp['Height'], tmp['H_bref'] = al, he, he/bref
                    df = pd.concat([df, tmp])

    return df