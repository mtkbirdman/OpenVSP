import sys
import os

# ../bin/AnalysisVSPAERO.py をモジュールとしてインポート
sys.path.append(os.path.join('../..')) # 親ディレクトリをモジュール探索パスに追加

from src.AnalysisVSPAERO import validate_vsp3_for_stability_derivatives, vsp_stability_derivatives

if __name__ == '__main__':

    validation_report = validate_vsp3_for_stability_derivatives(
        r'../models/G103A/G103A.vsp3',
        verbose=2,
    )

    print('\nValidation result')
    print('-----------------')
    print('passed:', validation_report['passed'])

    if validation_report['errors']:
        print('\nErrors')
        for error in validation_report['errors']:
            print('-', error['code'], error['message'])

    if validation_report['warnings']:
        print('\nWarnings')
        for warning in validation_report['warnings']:
            print('-', warning['code'], warning['message'])

    if not validation_report['passed']:
        sys.exit(1)

    stability_report = vsp_stability_derivatives(
        r'../models/G103A/G103A.vsp3',
        alpha=2.0,
        mach=0.1,
        reynolds=4.4e6,
        verbose=2,
    )

    print('\nStability derivative result')
    print('---------------------------')
    print('passed:', stability_report['passed'])

    if stability_report['errors']:
        print('\nErrors')
        for error in stability_report['errors']:
            print('-', error['code'], error['message'])

    if stability_report['warnings']:
        print('\nWarnings')
        for warning in stability_report['warnings']:
            print('-', warning['code'], warning['message'])

    if not stability_report['passed']:
        sys.exit(1)

    print('\nStability derivatives')
    print(stability_report['derivatives'].T)
