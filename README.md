# OpenVSP / VSPAERO analysis utilities

OpenVSP / VSPAERO を Python から実行し、航空機の空力解析・安定微係数解析・定常旋回トリム・6自由度応答・垂直尾翼／上反角の設計検討までをつなぐための解析コードと実行例です。

> [!NOTE]
> このリポジトリは NASA/OpenVSP 本体ではありません。OpenVSP 本体と Python API は別途インストールしてください。

## できること

主な解析の流れは次のとおりです。

```text
OpenVSP model (.vsp3)
        │
        ├─ VSPAERO sweep / pitch-trimmed polar / ground effect
        │      └─ src/AnalysisVSPAERO.py
        │
        └─ VSPAERO stability derivatives
               └─ src/AnalysisVSPAERO.py
                      │
                      ▼
                   .stab
                      │
                      ├─ .stab parsing / linear aero model
                      │      └─ src/VSPAEROStab.py
                      │
                      ├─ steady turn trim
                      │      └─ src/TrimTurnSolver.py
                      │
                      ├─ rudder-step / crosswind-gust response
                      │      └─ src/RollRudderGain.py
                      │
                      └─ Vv–Gamma design study
                             └─ src/VvGammaChart.py
```

具体的には、次の処理を含みます。

- VSPAERO の迎角・Mach・Reynolds 数 sweep
- Control Surface Group の設定
- エレベータを用いた pitch-trimmed polar
- ground effect sweep
- VSPAERO stability derivatives の実行と `.stab` 取得
- `.stab` の読み取りと線形空力モデルの評価
- 定常高度維持旋回・定常滑空旋回・ラダー限界旋回のトリム計算
- ラダーステップ入力に対する横・方向応答と 6DoF 応答
- 1-cosine 横風突風に対する 6DoF 応答
- 垂直尾翼容積比 `Vv` と翼端たわみ量を変化させた VSPAERO stability sweep
- `Vv–Gamma_eff` 設計チャートの後処理・描画
- wake iteration / wake node 数の収束確認
- reference Wing の実3次元 edge size を基準にした VSPAERO surface mesh equalization
- `.vsp3` モデルの複数方向スクリーンショット撮影と PNG グリッド合成
- 複数 `.vsp3` の連続撮影と GIF / MP4 アニメーション作成

## 必要な環境

### OpenVSP

OpenVSP 本体と、その OpenVSP build に対応する Python API が必要です。

OpenVSP Python API は、**OpenVSP がビルドされた Python バージョンと同じ Python バージョンを使用する必要があります**。

- OpenVSP: <https://openvsp.org/>
- OpenVSP Python API documentation: <https://openvsp.org/pyapi_docs/latest/>
- OpenVSP documentation: <https://openvsp.org/docs.shtml>

まず、使用する Python 環境で次が成功することを確認してください。

```python
import openvsp as vsp
print(vsp.GetVSPVersion())
```

### Python packages

このリポジトリは Python package として配布していません。リポジトリルートを作業ディレクトリとして、そのまま `src` を import して使用します。

OpenVSP Python API に加えて、解析内容に応じて次の package を使用します。

```bash
python -m pip install numpy pandas scipy matplotlib pillow imageio-ffmpeg jupyter
```

主な依存関係は次のとおりです。

- NumPy
- pandas
- SciPy
- Matplotlib
- Pillow（OpenVSP スクリーンショットの回転・合成）
- imageio-ffmpeg（GIF / MP4 のエンコード）
- Jupyter（Notebook を使用する場合）

## Quick start

リポジトリを clone し、**リポジトリルートから**実行します。

```bash
git clone https://github.com/mtkbirdman/OpenVSP.git
cd OpenVSP
```

### Tests を実行する

自動検証は `tests/` に統一しています。リポジトリルートから次を実行します。

```bash
python -m pytest
```

現在、次の test があります。

```text
tests/test_vspaero_mesh.py
tests/test_turn_trim.py
tests/test_sweep.py
tests/test_control_surface.py
tests/test_ground_effect.py
tests/test_trimmed_polar.py
tests/test_stability_derivatives.py
tests/test_openvsp_view_capture.py
```

`test_vspaero_mesh.py` の入力 validation、`test_turn_trim.py`、`test_openvsp_view_capture.py` の通常テストは OpenVSP を実行せずに検証できます。`test_turn_trim.py` は同梱されている `SampleGlider.stab` を使い、定常滑空旋回、定常高度維持・ラダーのみ旋回、ラダー限界旋回について、力・モーメント残差、無推力条件、高度維持条件、ラダー限界の選択処理を assertion で確認します。

`test_vspaero_mesh.py` の実 mesh equalization を含む残りのケースは OpenVSP / VSPAERO integration test です。単に解析が終了することだけでなく、sweep の入力と出力の対応、aileron による rolling moment の反転、ground effect による induced drag の低下、pitch trim 後の `CMytot ≈ 0`、stability-derivative preflight と出力データの有限性などを確認します。OpenVSP Python API が import できない環境では、これらの test は pytest 上で skip されます。

`STABILITY_ADJOINT` の実計算は非常に時間がかかる場合があるため `slow` marker を付け、通常の `python -m pytest` から除外しています。必要な場合だけ次のように実行します。

```bash
python -m pytest -m slow tests/test_stability_derivatives.py
```

スクリーンショットの実 GUI 結合テストは、画面を開ける環境で明示的に有効化します。

```bash
RUN_OPENVSP_GUI_TESTS=1 python -m pytest -m slow tests/test_openvsp_view_capture.py
```

`examples/scripts/` に置かれていた手動確認 script はすべて pytest 化したため、同ディレクトリは廃止しました。

## Notebook

より詳細な解析例は `examples/notebooks/` にあります。

| Notebook | 内容 |
|---|---|
| `trimmed_polar/example_trimmed_polar.ipynb` | G103A の pitch-trimmed polar |
| `trimmed_ground_effect/ground_effect.ipynb` | SampleGlider の ground effect |
| `turn-maneuver/example_turn_trim.ipynb` | G103A の定常旋回・横滑り・ラダー限界旋回 |
| `turn-maneuver/example_rudder_step_response.ipynb` | `.stab` を用いたラダーステップ 6DoF 応答 |
| `turn-maneuver/example_crosswind_gust_response.ipynb` | 1-cosine 横風突風に対する 6DoF 応答 |
| `vv_gamma_chart/test_vv_gamma_chart.ipynb` | `Vv × wing tip deflection` の stability sweep |
| `vv_gamma_chart/plot_vv_gamma_chart.ipynb` | `Vv–Gamma_eff` chart の後処理・描画 |
| `wake_convergence/run_wake_convergence.ipynb` | WakeNumIter / NumWakeNodes の収束計算 |
| `wake_convergence/plot_wake_convergence.ipynb` | wake convergence 結果の可視化 |
| `mesh_equalization/run_mesh_equalization.ipynb` | reference Wing の実3次元 edge size を基準にした surface mesh equalization |
| `view_capture/example_view_capture.ipynb` | G103A の上面・左側面・正面・鳥瞰図を1枚の PNG に合成 |
| `view_capture/example_view_animation.ipynb` | 複数 `.vsp3` を同じ視点で連続撮影し、MP4 / GIF を作成 |

Notebook は、基本的にリポジトリ内の `src/` と `examples/models/` を利用する構成です。

### OpenVSP view capture

`src/OpenVSPViewCapture.py` の `capture_vsp3_views()` は、GUI 対応版 OpenVSP
を専用の子プロセスで起動し、指定した標準ビューを個別に PNG 撮影してから
1枚に合成します。子プロセス内で `openvsp_config.LOAD_GRAPHICS` と
`LOAD_FACADE` を `openvsp` の import より前に設定するため、既存の Jupyter
カーネルにある OpenVSP のモデル状態を変更しません。撮影中は GUI が表示され、
撮影終了時に自動的に閉じます。

```python
from src.OpenVSPViewCapture import capture_vsp3_views

capture_vsp3_views(
    "examples/models/G103A/G103A.vsp3",
    "G103A_four_views.png",
    size=(1600, 1200),
    render_mode="hidden",
)
```

既定配置は、左上が機首下向きの上面図、右上が機首下向きの左側面図、
左下が正面図、右下が左鳥瞰図です。OpenVSP の `ScreenGrab` が対応する
保存形式に合わせ、出力は PNG のみに制限しています。`render_mode` は
`preserve`、`wire`、`hidden`、`shade`、`texture` から選べます。`preserve`
以外を指定しても、モデル内で非表示の Geom は非表示のままです。
レンダーモードの変更には OpenVSP 3.50 系にも存在する `SET_SHOWN` を使用します。
各パネルの撮影前に OpenVSP のビューポート寸法を出力寸法へ合わせるため、
3.50 系の `ScreenGrab` でも縦横比を保ちます。90度回転するビューは撮影時点で
幅と高さを入れ替え、回転後に再拡大・再縮小しません。

複数モデルでは `create_vsp3_animation()` に順序付きのパスを渡します。
同じパスの重複もそのまま1フレームとして扱います。出力拡張子が `.gif`
なら GIF、`.mp4` なら H.264 MP4 を作成します。

```python
from src.OpenVSPViewCapture import create_vsp3_animation

result = create_vsp3_animation(
    ["case_001.vsp3", "case_002.vsp3", "case_003.vsp3"],
    "mesh_change.mp4",
    size=(1600, 1200),
    render_mode="hidden",
    fps=3,
)
```

既定では `mesh_change_frames/` に連番 PNG と撮影条件の manifest を残します。
同じモデル・表示条件で再実行すると正常な既存フレームを再利用し、欠損・破損
フレームだけを撮り直します。`fps` や GIF / MP4 の違いは撮影画像を変えないため、
同じ `frames_dir` から速度や形式を変えて再エンコードできます。モデルの一覧を
処理する間は OpenVSP GUI を1回だけ起動し、完了時に閉じます。
H.264 MP4 はアルファチャンネルを保持できないため、MP4 と
`transparent_background=True` の組み合わせは受け付けません。

## ディレクトリ構成

```text
OpenVSP/
├─ README.md
├─ src/
│  ├─ AnalysisVSPAERO.py
│  ├─ VSPAEROMesh.py
│  ├─ VSPAEROStab.py
│  ├─ OpenVSPViewCapture.py
│  ├─ TrimTurnSolver.py
│  ├─ RollRudderGain.py
│  ├─ VvGammaChart.py
│  ├─ util.py
│  ├─ ISAspecification.py
│  └─ CST.py
└─ examples/
   ├─ models/
   │  ├─ G103A/
   │  ├─ SampleGlider/
   │  └─ Boeing_777-9x_mod/
   └─ notebooks/
```

## 主要モジュール

### `src/AnalysisVSPAERO.py`

OpenVSP Python API を通して VSPAERO 解析を実行します。

主な公開関数は次のとおりです。

- `vsp_sweep()` — 通常の VSPAERO sweep
- `vsp_trimmed_sweep()` — 迎角を固定し、`CMytot = 0` となるエレベータ舵角を求める trimmed polar
- `vsp_sweep_wig()` — ground effect を含む sweep / trimmed sweep
- `vsp_stability_derivatives()` — stability derivatives の計算
- `validate_vsp3_for_stability_derivatives()` — stability analysis 前のモデル設定確認
- `make_CDo_correction()` — skin-friction / profile drag の補正

### `src/VSPAEROMesh.py`

VSPAERO surface tessellation の **physical mesh equalization** を担当します。

主な公開関数は `equalize_vspaero_tessellation()` です。reference Wing に指定した `Tess_W` を適用し、その実3次元 W-edge median を target mesh size として、各 surface Geom の `Tess_W`、`SectTess_U` / `Tess_U`、active Wing end cap の `CapUMinTess` を調整します。clustering parameter は変更しません。

```python
from src.VSPAEROMesh import equalize_vspaero_tessellation

report = equalize_vspaero_tessellation(
    input_vsp3_path="examples/models/G103A/G103A.vsp3",
    output_vsp3_path="results/G103A.equalized.vsp3",
    reference_wing_name="WingGeom",
    reference_tess_w=65,
)

print(report[[
    "geom_name",
    "u_target_ratio",
    "w_target_ratio",
    "cap_target_ratio",
]])
```

処理は、geometry から初期 tessellation count を決めた後、OpenVSP が実際に生成した 3D tessellation edge を再測定し、指定 tolerance を外れる方向だけを少数回補正します。VSPAERO convergence study や post-intersection mesh diagnostics はこの module の責務に含めません。

### `src/VSPAEROStab.py`

VSPAERO `.stab` ファイルの読み取りを担当します。

- reference quantities
- base flight condition
- base aerodynamic coefficients
- stability derivatives
- Control Surface Group と `ConGrp_*` の対応

を共通形式へ変換し、定常旋回 solver、6DoF simulation、設計チャートから共通利用します。

### `src/TrimTurnSolver.py`

`.stab` の線形空力モデルを使って定常旋回を解きます。

- `solve_steady_level_turn()`
- `solve_steady_gliding_turn()`
- `solve_rudder_limit_turn()`

solver の body axes は `+x` forward、`+y` right、`+z` down です。OpenVSP の body axes とは向きが異なるため、`.stab` の空力係数を solver 内で変換しています。

### `src/RollRudderGain.py`

`.stab` を飛行力学モデルへ接続し、ラダー入力・横風突風に対する応答を計算します。

主な処理は次のとおりです。

- 線形横・方向応答
- ラダーステップ 6DoF simulation
- 1-cosine crosswind gust 6DoF simulation
- 時刻歴 CSV 出力
- 時刻歴 plot
- roll response / gust response の評価指標

### `src/VvGammaChart.py`

垂直尾翼容積比 `Vv` と主翼の等価上反角 `Gamma_eff` を設計変数として、横・方向特性を比較するための処理をまとめています。

大きな流れは次の2段階です。

```text
base .vsp3
  ↓
Vv / wing tip deflection を変更
  ↓
VSPAERO stability derivatives
  ↓
case .vsp3 / .stab
```

```text
case .vsp3 / .stab
  ↓
Gamma_eff・安定性・旋回・6DoF 指標を後処理
  ↓
Vv–Gamma_eff chart
```

## OpenVSP model の前提

解析関数は任意の `.vsp3` を無条件に処理できるものではありません。特に stability derivatives や trim では、OpenVSP GUI 側で次の設定を確認してください。

- VSPAERO reference area / chord / span
- reference point / CG
- VSPAERO geometry set
- Thin / Thick geometry の扱い
- Control Surface Group
- `Control Group Angles` の有効化状態
- control surface gain と舵角符号
- VSPAERO analysis method
- wake settings

G103A の trimmed-polar workflow では、既存モデルの `WingGeom` と `ELEVATOR_GROUP` を使用します。`ELEVATOR_GROUP` は VSPAERO の `Control Group Angles` で有効化されている必要があります。`vsp_trimmed_sweep()` はこの設定を自動変更せず、無効な場合は解析前にエラーにします。

Control Surface Group の名前は `.stab` の `ConGrp_*` と実際の舵を対応させるためにも重要です。G103A と SampleGlider の例では、次の名前を使用しています。

```text
AILERON_GROUP
ELEVATOR_GROUP
RUDDER_GROUP
```

## `.vsp3` と VSPAERO 生成物

`examples/models/` には、解析の再現や後処理例に必要な `.vsp3` と一部の VSPAERO 出力を含めています。

特に `.stab` は、VSPAERO を再実行せずに `VSPAEROStab.py`、`TrimTurnSolver.py`、`RollRudderGain.py` の処理を確認するための reference result として使用します。

`.vspaero`、`.stab`、`.flt`、`.polar` などを同一解析 run の一式として扱う場合は、**Mach、AoA、Reynolds 数、reference quantities、CG、wake settings、control groups などの解析条件をそろえてください**。

`.history`、`.lod`、`.vspgeom` などの一時的な VSPAERO 出力の多くは `.gitignore` の対象です。

## 数値精度と wake convergence

VSPAERO の安定微係数は wake discretization / wake iteration の影響を受けます。特に `CL_p`、`Cl_beta`、`Cn_beta` などを設計指標として使用する場合、1つの wake 設定を無条件に採用せず、対象機体・解析条件ごとに収束を確認してください。

このため、`examples/notebooks/wake_convergence/` に `WakeNumIter` と `NumWakeNodes` を変化させる Notebook を置いています。

計算時間が大きい場合も、単に mesh / wake を細かくすればよいとは限りません。必要な微係数が十分収束する最小限の設定を選ぶことを推奨します。

## 関連記事

実装の背景、OpenVSP GUI の設定、空力・飛行力学の理論は以下の記事で詳しく説明しています。

### OpenVSP / VSPAERO

- [【まとめ】OpenVSP入門](https://mtkbirdman.com/openvsp-index)
- [【OpenVSP入門】インストール](https://mtkbirdman.com/openvsp-install)
- [【Windows】OpenVSP Python APIのインストール](https://mtkbirdman.com/openvsp-python-api-installation)
- [OpenVSPのPythonAPIでポーラーカーブを計算する](https://mtkbirdman.com/openvsp-python-api-sweep-analysis)
- [OpenVSPのPythonAPIで舵角を設定して計算する](https://mtkbirdman.com/openvsp-python-api-control-surface)
- [OpenVSPのPythonAPIでトリムドポーラーを計算する](https://mtkbirdman.com/openvsp-python-api-trimed-polar)
- [OpenVSPのPythonAPIで地面効果を計算する](https://mtkbirdman.com/openvsp-python-api-ground-effect)
- [OpenVSP の Python API で安定微係数を計算する](https://mtkbirdman.com/openvsp-python-api-stability-analysis)
- [GROB G 103 Twin II の OpenVSPモデル](https://mtkbirdman.com/openvsp-grob-g-103-twin-ii-example)
- [鳥コン滑空機を想定したSampleGliderのOpenVSPモデル](https://mtkbirdman.com/openvsp-br-sample-glider-example)

### 定常旋回・6DoF

- [航空機の定常旋回運動](https://mtkbirdman.com/turn-maneuver-basis)
- [航空機の定常旋回パラメータを計算する Python スクリプト](https://mtkbirdman.com/turn-maneuver-trim-parameter-python-script)
- [航空機のラダーのみ旋回の 6 自由度剛体運動をシミュレーションする Python スクリプト](https://mtkbirdman.com/turn-maneuver-only-rudder-roll-python-script)
- [航空機の横風突風応答の 6 自由度剛体運動をシミュレーションする Python スクリプト](https://mtkbirdman.com/turn-maneuver-crosswind-gust-python-script)

### 垂直尾翼・上反角設計

- [ラダーのみを用いて旋回する航空機の垂直尾翼・上反角設計](https://mtkbirdman.com/vv-gamma-chart)

## 注意事項

- VSPAERO の結果は、モデル形状だけでなく reference quantities、CG、control group、解析 method、wake settings に依存します。
- `.stab` を用いる後処理は、その基準状態近傍の線形空力モデルです。大迎角・大横滑り角・大舵角まで線形外挿した結果は、収束していても物理的妥当性を別途確認してください。
- `TrimTurnSolver.py`、`RollRudderGain.py` は `.stab` に記録された長さ・速度・密度などを暗黙に単位変換しません。質量・慣性モーメントなどの入力値を同じ単位系にそろえてください。
- OpenVSP / VSPAERO の API や出力形式はバージョンによって変わる可能性があります。

## License

現時点では、このリポジトリ全体に対する `LICENSE` ファイルは置いていません。

また、`examples/models/` に含まれる flight manual、airfoil data、その他第三者由来の資料・データについては、それぞれの権利者・配布元の条件に従ってください。
## 技術・V&V 統合ベースライン

過去の OpenVSP / VSPAERO mesh 診断、optimizer、rule-calibration dataset を含む技術履歴は、`docs/archive/OpenVSP_VSPAERO_integrated_technical_baseline_2026-09-23_v2.md` に保存しています。現在の実装は `src/VSPAEROMesh.py` の physical mesh equalization と `examples/notebooks/mesh_equalization/` を基準にしてください。
