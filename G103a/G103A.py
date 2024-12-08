import sys
import os

# ../bin/AnalysisVSPAERO.py をモジュールとしてインポート
sys.path.append(os.path.join(os.path.dirname(__file__), '..')) # 親ディレクトリをモジュール探索パスに追加
from bin.AnalysisVSPAERO import vsp_sweep

if __name__ == '__main__':
    alpha = list(range(-4, 13, 2))  # -4 to 12 with step of 2
    mach = [0.1]
    vsp_sweep(alpha=alpha, mach=mach)
