# OpenVSP / VSPAERO Python解析ツール

OpenVSP / VSPAEROをPythonから利用し、航空機の空力解析・飛行力学解析・可視化を行うリポジトリです。解析コード、Notebook、サンプルモデル、テストを収録しています。

> **対象と制約**：このリポジトリはOpenVSP本体の配布物ではありません。OpenVSP本体、互換性のあるPython API、VSPAERO実行ファイルは必要な機能に応じて別途用意してください。解析結果の妥当性は個別に検証してください。

## 1. 主な機能

| 分野 | 機能 | 主な実装 |
|---|---|---|
| VSPAERO解析 | 迎角・Mach数・Reynolds数のSweep、縦トリム、地面効果、安定微係数 | `src/openvsppy/vspaero/` |
| 空力・飛行力学 | 寄生抗力補正、定常旋回トリム、ラダー応答、横風突風の6自由度応答 | `src/openvsppy/aerodynamics/` |
| メッシュ | 3次元エッジ長に基づくVSPAEROメッシュ分割の調整 | `src/openvsppy/vspaero/VSPAEROMesh.py` |
| 可視化 | `.vsp3`の複数方向撮影、VSPAERO ADBの圧力分布4面図、GIF・MP4アニメーション | `src/openvsppy/visualization/` |

解析モデル・Notebook・設定例は`examples/`、回帰テストに必要な固定データは`tests/fixtures/`で管理します。

## 2. 動作環境

- **Python**：`pyproject.toml`の指定は3.10以上。ただし、OpenVSP Python APIを使用する場合は、**インストールされたOpenVSP 3.50.5 APIのビルド対象Pythonと互換性のある環境**が必要です。Pythonのバージョンだけでは共有ライブラリの動作まで保証されません。
- **OpenVSP / VSPAERO**：ネイティブテストは3.50.5を対象としています。オフラインテストの一部はOpenVSPなしで実行できます。
- **主なPythonライブラリ**：NumPy、pandas、SciPy、Matplotlib。撮影・動画関連機能にはPillow、imageio-ffmpeg、NotebookにはJupyterが必要です。
- **実行環境**：以下はWindows PowerShellの例です。パスとConda環境名は利用環境に合わせてください。

公式情報：[OpenVSP](https://openvsp.org/)／[Python APIリファレンス](https://openvsp.org/pyapi_docs/latest/)／[ドキュメント](https://openvsp.org/docs.shtml)

### 既存リポジトリからの起動

```powershell
# 自分のリポジトリルートへ移動（実際のパスへ変更）
$repo = 'C:\path\to\OpenVSP'
Set-Location $repo
conda activate vsppytools.3.50.5
$env:PYTHONPATH = (Resolve-Path '.\src').Path

python -m pytest -q -rs
```

環境に合わせて依存ライブラリを追加する場合の例：

```powershell
python -m pip install numpy pandas scipy matplotlib pillow imageio-ffmpeg jupyter pytest
```

`pip install -e .`は必須ではありません。直接`src/`を`PYTHONPATH`へ追加できます。Notebookは先頭セルでチェックアウトの`src/`を`sys.path`に追加します。

OpenVSP APIの利用前には、対応する環境で次を確認してください。

```powershell
python -c "import openvsp as vsp; print(vsp.GetVSPVersion())"
```

## 3. ディレクトリ構成

```text
OpenVSP/
├── README.md
├── pyproject.toml              # パッケージ・pytest設定
├── src/openvsppy/
│   ├── core/                   # 共通ユーティリティ・大気モデル
│   ├── geometry/               # 形状関数
│   ├── vspaero/                # 解析API・メッシュ・安定微係数
│   ├── aerodynamics/           # 飛行力学・抗力補正・設計計算
│   └── visualization/          # OpenVSPの撮影・画像処理
├── examples/
│   ├── models/                 # 解析モデルと付属データ
│   └── notebooks/              # 解析・可視化の実行例
├── tests/                      # オフライン・ネイティブ・GUIテスト
│   └── fixtures/               # テストに必要な固定データ
├── tools/
│   └── run_native_pytest.py    # ネイティブテストの隔離実行
└── tree.ipynb                  # ディレクトリ構成の確認用Notebook
```

`.gitignore`の対象となる解析結果、Notebook生成物、旧VV実行データなどは通常のGitコミットに含めません。必要なテストfixtureのみを収録しています。

## 4. Python APIとNotebookの使用例

リポジトリルートで`src/`がインポート可能な状態から利用します。

```python
from openvsppy.vspaero.VSPAEROMesh import equalize_vspaero_tessellation

report = equalize_vspaero_tessellation(
    input_vsp3_path="examples/models/G103A/G103A.vsp3",
    output_vsp3_path="results/G103A.equalized.vsp3",
    reference_wing_name="WingGeom",
    reference_tess_w=65,
)
print(report[["geom_name", "u_target_ratio", "w_target_ratio", "cap_target_ratio"]])
```

この例はOpenVSP Python APIと解析モデルを必要とします。実行前に`results/`などの保存先を作成してください。

主なNotebookは以下のとおりです。各Notebookの入力・前提条件を確認してください。

| ディレクトリ | 主題 |
|---|---|
| `trimmed_polar/` | 縦トリムを考慮したポーラー |
| `trimmed_ground_effect/` | 地面効果 |
| `turn-maneuver/` | 定常旋回・ラダー応答・横風突風 |
| `vv_gamma_chart/` | 垂直尾翼容積比・有効上反角の設計検討 |
| `wake_convergence/` | Wake反復・後流節点の感度 |
| `mesh_equalization/` | メッシュ分割の調整 |
| `view_capture/` | 多方向撮影・アニメーション |

これらは`examples/notebooks/`内にあります。一部のNotebookは実行結果や`.stab`など、ZIPに含まれない追加ファイルを必要とします。

### VSPAERO圧力分布の4面図（ADBを直接読み込み）

`src/openvsppy/visualization/VSPAEROCpViewCapture.py` は、**OpenVSP 3.50.5で生成した定常Version 3の`.adb`**をPythonで読み込み、NumPyの正射影・Z-bufferとMatplotlibのAggバックエンドを使って圧力分布をPNGにします。OpenVSP GUI、`vspaero.exe`、Viewer、C++コンパイル、PyVista/VTKは不要です。

同じ場所に`Aircraft.adb.cases`が必要です。解析結果を再実行せず、ケース番号（**1始まり**）を選択できます。

```python
from openvsppy.visualization.VSPAEROCpViewCapture import (
    capture_adb_cp_views,
    list_adb_cases,
)

for item in list_adb_cases('Aircraft.adb'):
    print(item.index, item.mach, item.alpha_deg, item.beta_deg, item.comment)

capture_adb_cp_views(
    'Aircraft.adb',
    'Aircraft_Cp.png',
    case=1,
    size=(1600, 1200),
    cp_limits=(-1.5, 1.0),  # None: 表示対象の最小値・最大値から自動設定
    align_orthographic_views=True,
    reflect_symmetry='auto',
)
```

既定配置は以前の形状4面図と同じく、上面（機首下向き）／左側面（機首下向き）／正面／左鳥瞰です。直交投影図は共通の画素縮尺とモデル中心を使用し、鳥瞰図はセルに独立にフィットします。全ビューに共通のCpカラーバーを表示します。撮影方向・表示範囲・背景透過などは関数引数で変更できます。Notebook例は`examples/notebooks/view_capture/example_adb_cp_capture.ipynb`です。

**制約**：初期対応は、3.50.5のWriter/Viewerに整合するVersion 3、定常VLM/Panel、**ロータ・ノズルレコードのない**ADBです。非定常、古いVersion 3レイアウト、混在Thin/Thickの面単位の圧力意味の識別、舵面のViewer上の追加変形、Wake描画には対応しません。VLMは`Delta Cp`、Panelは`Cp`と表示しますが、混在モデルの意味は別途確認が必要です。対称モデルの反転はY=0の半機体を検出した場合に限定します。ビューの色はViewerのPer-Triangleカラーマップに基づきますが、Viewerの照明やアンチエイリアスを完全に再現するものではありません。

ソース参照：OpenVSP 3.50.5 `src/vsp_aero/Solver/VSP_Solver.C` の`WriteOutAerothermalDatabaseHeader/Geometry/Solution`と、`src/vsp_aero/Viewer/glviewer.C`の`LoadMeshData`、`UpdateMeshData`、`LoadExistingSolutionData`、`DrawShadedSolutionPerTri`、`percent_to_rgb`。

## 5. テスト

通常の回帰テストでは、OpenVSP API、VSPAERO、GUIを使用するテストは**選択対象外（DESELECTED）**です。解析環境なしで実行できる部分と、ネイティブ実行を伴う試験結果を区別してください。

```powershell
python -m pytest -q -rs
```

OpenVSP API、VSPAERO、GUIの実行が可能な環境では、必要な区分だけを明示的に選択します。

```powershell
$env:VSPAERO_BINARY = 'C:\OpenVSP-3.50.5\vspaero.exe'  # 実際のパスに変更
python tools/run_native_pytest.py --marker all-native --list-only
python tools/run_native_pytest.py --marker native_openvsp --timeout 600
python tools/run_native_pytest.py --marker native_vspaero --timeout 3600
python tools/run_native_pytest.py --marker native_solver --timeout 240
# GUI環境でのみ実行
python tools/run_native_pytest.py --marker gui --timeout 300
```

`tools/run_native_pytest.py`はネイティブ試験を別プロセスで隔離実行し、タイムアウト時に子プロセスも終了させます。結果はGit管理外の`native_test_results/`に保存します。自動試験で検査できた範囲と、実機による解析性能の検証範囲を混同しないでください。

## 6. 注意事項・ライセンス

- 解析値の妥当性は、形状・解析モード・座標系・基準量・重心・Control Surface Group、メッシュ・Wake条件、参照データおよび数値収束に依存します。
- `.stab`による線形モデルの外挿や、VSPAERO内蔵の寄生抗力と外部CD0補正の組合せには適用範囲上の制約があります。
- OpenVSP、Python API、VSPAEROの仕様はバージョンで異なる場合があります。
- このリポジトリには現時点でルートの`LICENSE`ファイルがありません。権利・再配布条件は別途確認してください。第三者由来のモデル、フライトマニュアル、エアフォイル等には個別の利用条件が適用される場合があります。
