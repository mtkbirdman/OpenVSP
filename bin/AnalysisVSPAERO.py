# https://openvsp.org/pyapi_docs/latest/openvsp.html
# https://github.com/OpenVSP/OpenVSP/blob/fc4a97b75c92399fecf74b21f2a57087f89e12b7/examples/scripts/TestAnalysisVSPAERO.vspscript

import openvsp as vsp

def vsp_sweep(vsp, alpha, mach):

    print('\n-> Calculate alpha & mach sweep analysis\n')

    # //==== Analysis: VSPAero Compute Geometry to Create Vortex Lattice DegenGeom File ====//
    
    # Set defaults
    compgeom_name = 'VSPAEROComputeGeometry'
    vsp.SetAnalysisInputDefaults(compgeom_name)

    # List inputs, type, and current values
    print(compgeom_name)
    vsp.PrintAnalysisInputs(compgeom_name)
    print('')

    # Execute
    print('\tExecuting...')
    compgeom_resid = vsp.ExecAnalysis(compgeom_name)
    print('\tCOMPLETE')

    # Get & Display Results
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
    alpha_npts = len(alpha)
    vsp.SetDoubleAnalysisInput(analysis_name, "AlphaStart", [alpha[0]], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "AlphaEnd", [alpha[-1]], 0)
    vsp.SetIntAnalysisInput(analysis_name, "AlphaNpts", [alpha_npts], 0)
    mach_npts = [len(mach)]
    vsp.SetDoubleAnalysisInput(analysis_name, 'MachStart', mach, 0)
    vsp.SetIntAnalysisInput(analysis_name, 'MachNpts', mach_npts, 0)

    vsp.Update()

    # List inputs, type, and current values
    print(analysis_name)
    vsp.PrintAnalysisInputs(analysis_name)
    print('')

    # Execute
    print('\tExecuting...')
    rid = vsp.ExecAnalysis(analysis_name)
    print('COMPLETE')

    # Get & Display Results
    # vsp.PrintResults(rid)

def MyFindParm(cs_group_container_id, parm_name):

    # Find param_id from cs_group_container_id and parm_name
    parm_ids = vsp.FindContainerParmIDs( cs_group_container_id )
    for parm_id in parm_ids:
        if parm_name == vsp.GetParmName(parm_id):
            return parm_id
    return ''

def set_control_surface(geom_name, deflection, cs_group_name, sub_id=0, gains=(1,1)):

    # Set control surface group
    group_index  = vsp.CreateVSPAEROControlSurfaceGroup()
    group_name = 'ControlSurfaceGroup_' + str(group_index)
    vsp.SetVSPAEROControlGroupName(cs_group_name, group_index)

    # Add subsurface to control surface group
    cs_name_vec = vsp.GetAvailableCSNameVec(group_index)
    selcted = [1 + i for i, cs_name in enumerate(cs_name_vec) if geom_name in cs_name]
    vsp.AddSelectedToCSGroup(selcted, group_index)
    
    # Set gain (deferential)
    print('\n', cs_group_name, ':', vsp.GetActiveCSNameVec(group_index))
    print('\t', f"{'container_id':14s}", f"{'parm_name':24s}", f"{'group_name':22s}", f"{'parm_id':14s}", f"{'value':8s}")
    geom_id = vsp.FindGeomsWithName(geom_name)[0]
    cs_group_container_id = vsp.FindContainer('VSPAEROSettings', 0)    
    for i, gain in enumerate(gains):
        parm_name = 'Surf_' + vsp.GetSubSurf(geom_id, sub_id) + '_' + str(i) + '_Gain'
        parm_id = MyFindParm(cs_group_container_id, parm_name)
        # parm_id = vsp.FindParm(cs_group_container_id, parm_name, group_name) # なぜか反応しない...
        if parm_id != '':
            vsp.SetParmVal(parm_id, gain)
            print('\t', f'{cs_group_container_id:14s}', f'{parm_name:24s}', f'{group_name:22s}', f'{parm_id:14s}', f'{gain:8.3f}')
        else:
            print('\t', 'Failed to find ' + parm_name)

    # Set defrection angle
    parm_name = 'DeflectionAngle'
    parm_id = vsp.FindParm(cs_group_container_id, parm_name, group_name)
    vsp.SetParmVal(parm_id, deflection)

    # Check if deflection angle successfully applied
    if deflection == vsp.GetParmVal(parm_id):
        print('\t', f'{cs_group_container_id:14s}', f'{parm_name:24s}', f'{group_name:22s}', f'{parm_id:14s}', f'{deflection:8.3f}')
    else:
        print('Failed to set deflection angle')
        exit()

    return geom_id, group_index