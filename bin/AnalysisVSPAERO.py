# https://openvsp.org/pyapi_docs/latest/openvsp.html
# https://github.com/OpenVSP/OpenVSP/blob/fc4a97b75c92399fecf74b21f2a57087f89e12b7/examples/scripts/TestAnalysisVSPAERO.vspscript

import openvsp as vsp

def vsp_sweep(alpha, mach):
    print('-> Calculate alpha & mach sweep analysis\n')

    # //==== Analysis: VSPAero Single Point ====//

    # Close and open the file
    vsp.ClearVSPModel()
    vsp.Update()
    vsp.ReadVSPFile('G103A.vsp3')  # Sets VSP3 file name
    vsp.Update()

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

# Usage sample
# if __name__ == '__main__':
#     alpha = list(range(-4, 13, 2))  # -4 to 12 with step of 2
#     mach = [0.1]
#     vsp_sweep(alpha=alpha, mach=mach)
