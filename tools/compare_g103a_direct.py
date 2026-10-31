"""Compare one native API Run with the supplied 2026-10-31 direct-run evidence.

This reports differences, never automatically approves V&V or numerical convergence.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path

# Resolve src layout from this repository without editable installation.
import sys
_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT / 'src') not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT / 'src'))
from openvsppy.vv.VVPostprocess import read_history, read_polar


def hashfile(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    if len(sys.argv)!=2:
        raise SystemExit('Usage: python -m tools.compare_g103a_direct <Run/raw path>')
    run=Path(sys.argv[1]).resolve()
    evidence=Path(__file__).resolve().parents[1]/'vv_study/02_reference/20261031_api_vs_direct'
    api_ref=evidence/'api_failed_inputs'
    direct_ref=evidence/'direct_success_outputs'
    result={'run_directory':str(run),'binary_identity':'NOT_CHECKED_BY_THIS_COMPARISON',
            'hash_comparison':{},'coefficients':{},'history':{}}
    for ext in ('.vsp3','.vspgeom','.vspaero','.vkey','.csf'):
        file=run/('model_effective'+ext)
        ref=api_ref/('model_effective'+ext)
        result['hash_comparison'][ext]={
            'new_sha256':hashfile(file) if file.is_file() else None,
            'previous_api_sha256':hashfile(ref) if ref.is_file() else None}
        result['hash_comparison'][ext]['exactly_equal']=(file.is_file() and ref.is_file() and hashfile(file)==hashfile(ref))
    refpolar=read_polar(direct_ref/'model_effective.polar')[0]
    cur=run/'model_effective.polar'
    if cur.is_file():
        newpolar=read_polar(cur)[0]
        for k in ('CLtot','CDtot','CMytot','CLiw','CDiw','L/D'):
            a,b=refpolar.get(k),newpolar.get(k)
            result['coefficients'][k]={'direct':a,'api':b,
                                       'delta':float(b)-float(a) if a is not None and b is not None else None}
    refhist=read_history(direct_ref/'model_effective.history')
    newhist=run/'model_effective.history'
    result['history']['direct_last_iter']=int(refhist[-1]['Iter'])
    if newhist.is_file():
        h=read_history(newhist)
        result['history']['api_last_iter']=int(h[-1]['Iter']) if h else None
    print(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False))


if __name__=='__main__':main()
