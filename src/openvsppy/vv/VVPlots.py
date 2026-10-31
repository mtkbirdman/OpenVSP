"""Offline reviewer plots from frozen V&V result snapshots; no PASS/FAIL inference."""
from __future__ import annotations

import re
from pathlib import Path


def create_study_plots(snapshot, destination):
    """Plot raw coefficients, histories, and reference comparisons, with source IDs.

    Keep missing/mismatched units separate. Charts never reinterpret a native
    stability derivative as a canonical aerodynamic derivative.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    destination=Path(destination)
    destination.mkdir(parents=True,exist_ok=True)
    charts=[]
    families=sorted({case['test_id'] for case in snapshot['cases']})
    for test in families:
        cases=[case for case in snapshot['cases'] if case['test_id']==test and case['run_id']]
        by_quantity={}
        for case in cases:
            for q in case['quantities']:
                # A reprocessed polar should not duplicate its native counterpart.
                key=(q['qoi'],q['variant'],q['component'],q['unit'],q['status'])
                by_quantity.setdefault(key,{})[case['case_id']]=(case,q)
        for (qoi,variant,component,unit,status),items in sorted(by_quantity.items()):
            if len(items)<2:continue
            ordered=[items[c['case_id']] for c in cases if c['case_id'] in items]
            if len(ordered)<2:continue
            fig,ax=plt.subplots(figsize=(8.0,4.0))
            ax.plot(range(len(ordered)),[pair[1]['value'] for pair in ordered],'-o')
            ax.set_xticks(range(len(ordered)),[pair[0]['role'] for pair in ordered],rotation=35,ha='right')
            ax.set_ylabel(f'{qoi} [{unit}]')
            ax.set_title(f'{test}: {qoi} / {variant} / {component}\n{status}')
            ax.grid(True,alpha=0.3)
            fig.tight_layout()
            label=re.sub(r'[^a-zA-Z0-9_-]','_',f'{qoi}_{variant}_{component}_{unit}')[:100]
            image=destination/f'{test}_{label}.png'
            fig.savefig(image,dpi=140);plt.close(fig)
            charts.append({'test_id':test,'kind':'CASE_SWEEP','qoi':qoi,'path':str(image),
                'sources':[{'case_id':c['case_id'],'run_id':c['run_id'],
                            'quantity_id':q['quantity_id'],'dataset_id':q['dataset_id']} for c,q in ordered]})
        # A one-Case test still needs a reviewer figure. Plot native values
        # only when they share a unit; never imply a derivative conversion.
        if len(cases)==1:
            case=cases[0]
            by_unit={}
            for q in case['quantities']:
                if q['variant']=='DEFAULT' and q['status']=='AVAILABLE':
                    by_unit.setdefault(q['unit'],[]).append(q)
            for unit,qs in by_unit.items():
                if len(qs)<2:continue
                fig,ax=plt.subplots(figsize=(max(8,len(qs)*0.45),4.5))
                ax.bar(range(len(qs)),[q['value'] for q in qs])
                ax.set_xticks(range(len(qs)),[q['qoi'] for q in qs],rotation=70,ha='right')
                ax.set_ylabel(f'Native values [{unit}]')
                ax.set_title(f'{test} / {case["case_id"]} — native coefficients (not compared)')
                ax.grid(True,axis='y',alpha=.3);fig.tight_layout()
                image=destination/f'{test}_native_values_{re.sub("[^A-Za-z0-9_-]","_",unit)}.png'
                fig.savefig(image,dpi=140);plt.close(fig)
                charts.append({'test_id':test,'kind':'SINGLE_CASE_VALUES','qoi':'native coefficients',
                    'path':str(image),'sources':[{'case_id':case['case_id'],'run_id':case['run_id'],
                    'dataset_id':q['dataset_id'],'quantity_id':q['quantity_id']} for q in qs]})
        if test.startswith('WAKE'):
            for case in cases:
                rows=case['history']
                if not rows:continue
                fields=[name for name in ('CLi','CDiw','CMiy','CLiw','CDi','l2_residual','max_residual')
                        if name in rows[0]]
                if not fields:continue
                fig,axes=plt.subplots(len(fields),1,figsize=(8.0,2.3*len(fields)),sharex=True)
                if len(fields)==1:axes=[axes]
                for ax,name in zip(axes,fields):
                    ax.plot([p['iteration'] for p in rows],[p[name] for p in rows],'-o',ms=2)
                    ax.set_ylabel(name);ax.grid(True,alpha=.3)
                axes[-1].set_xlabel('Iteration')
                fig.suptitle(f'{case["case_id"]}: native history')
                fig.tight_layout()
                image=destination/f'{case["case_id"]}_history.png'
                fig.savefig(image,dpi=140);plt.close(fig)
                charts.append({'test_id':test,'kind':'HISTORY','qoi':', '.join(fields),
                               'path':str(image),'sources':[{'case_id':case['case_id'],'run_id':case['run_id']}]})
    for test in ('DRAG-01','DRAG-02'):
        items=[r for r in snapshot.get('drag_polars',[]) if r['test_id']==test and r['cl'] is not None]
        if not items:continue
        if len({r['cl_qoi'] for r in items})!=1:continue
        items.sort(key=lambda x:x['cl'])
        fig,ax=plt.subplots(figsize=(8,4.5))
        if all(r['cd_native_wake'] is not None for r in items):
            ax.plot([r['cl'] for r in items],[r['cd_native_wake'] for r in items],'-o',label='Native CDiw')
        ax.plot([r['cl'] for r in items],[r['cd'] for r in items],'-s',label='Corrected '+items[0]['qoi'])
        ax.set_xlabel(items[0]['cl_qoi']);ax.set_ylabel('CD [1]')
        ax.set_title(f'{test}: derived drag polar (no engineering approval)')
        ax.grid(True,alpha=.3);ax.legend();fig.tight_layout()
        image=destination/f'{test}_corrected_drag_polar.png'
        fig.savefig(image,dpi=140);plt.close(fig)
        charts.append({'test_id':test,'kind':'DRAG_POLAR','qoi':items[0]['qoi'],
            'path':str(image),'sources':[{'dataset_id':r['dataset_id'],
                'quantity_id':r['quantity_id'],'polar_dataset_id':r['polar_dataset_id'],
                'cd0_dataset_id':r['cd0_dataset_id']} for r in items]})
    for test in sorted({x['test_id'] for x in snapshot['comparisons']}):
        for qoi in sorted({x['qoi'] for x in snapshot['comparisons'] if x['test_id']==test}):
            pairs=[x for x in snapshot['comparisons'] if x['test_id']==test and x['qoi']==qoi]
            if not pairs:continue
            if len({x['unit'] for x in pairs})!=1:continue
            fig,ax=plt.subplots(figsize=(8,4))
            indexes=list(range(len(pairs)))
            ax.plot(indexes,[x['predicted'] for x in pairs],'-o',label='VSPAERO')
            ax.plot(indexes,[x['reference'] for x in pairs],'-s',label='Reference')
            ax.set_xticks(indexes,[x['comparison_id'][-8:] for x in pairs],rotation=35)
            ax.set_ylabel(f'{qoi} [{pairs[0]["unit"]}]');ax.set_title(f'{test}: reference comparison')
            ax.legend();ax.grid(True,alpha=.3);fig.tight_layout()
            image=destination/f'{test}_comparison_{re.sub("[^A-Za-z0-9_-]","_",qoi)}.png'
            fig.savefig(image,dpi=140);plt.close(fig)
            charts.append({'test_id':test,'kind':'REFERENCE_COMPARISON','qoi':qoi,
                           'path':str(image),'sources':[{'comparison_id':x['comparison_id'],
                            'predicted_dataset':x['predicted_dataset'],
                            'reference_dataset':x['reference_dataset']} for x in pairs]})
    return charts
