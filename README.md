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
- VSPAERO `.adb` / `.vspgeom` を用いたパネル品質・`Cp` spike 診断
- `ThinWing / ThickWing / Hybrid / ThickAll` をユーザー指定した VSPAERO mesh-convergence optimization

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
python -m pip install numpy pandas scipy matplotlib jupyter
```

主な依存関係は次のとおりです。

- NumPy
- pandas
- SciPy
- Matplotlib
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
tests/test_vspaero_mesh_quality.py
tests/test_vspaero_mesh_optimizer.py
tests/test_vspaero_mesh.py
tests/test_turn_trim.py
tests/test_sweep.py
tests/test_control_surface.py
tests/test_ground_effect.py
tests/test_trimmed_polar.py
tests/test_stability_derivatives.py
```

`test_vspaero_mesh_quality.py`、`test_vspaero_mesh_optimizer.py`、`test_turn_trim.py` は OpenVSP を実行せずに検証できます。`test_turn_trim.py` は同梱されている `SampleGlider.stab` を使い、定常滑空旋回、定常高度維持・ラダーのみ旋回、ラダー限界旋回について、力・モーメント残差、無推力条件、高度維持条件、ラダー限界の選択処理を assertion で確認します。

残りは OpenVSP / VSPAERO integration test です。単に解析が終了することだけでなく、sweep の入力と出力の対応、aileron による rolling moment の反転、ground effect による induced drag の低下、pitch trim 後の `CMytot ≈ 0`、stability-derivative preflight と出力データの有限性などを確認します。OpenVSP Python API が import できない環境では、これらの test は pytest 上で skip されます。

`STABILITY_ADJOINT` の実計算は非常に時間がかかる場合があるため `slow` marker を付け、通常の `python -m pytest` から除外しています。必要な場合だけ次のように実行します。

```bash
python -m pytest -m slow tests/test_stability_derivatives.py
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
| `mesh_quality/example_mesh_quality.ipynb` | `.adb` / `.vspgeom` のパネル品質・`Cp` spike 診断 |
| `mesh_optimization/run_mesh_optimization.ipynb` | representation を指定した W→U mesh convergence search |
| `mesh_rule_dataset/run_mesh_rule_dataset.ipynb` | 全mesh自由度を対象にした append-only calibration dataset 生成 |
| `mesh_parameter_sweep/run_mesh_parameter_sweep.ipynb` | 4 representation を横断する manual mesh V&V / DOE |

Notebook は、基本的にリポジトリ内の `src/` と `examples/models/` を利用する構成です。

## ディレクトリ構成

```text
OpenVSP/
├─ README.md
├─ src/
│  ├─ AnalysisVSPAERO.py
│  ├─ VSPAEROMesh.py
│  ├─ VSPAEROMeshDataset.py
│  ├─ VSPAEROMeshQuality.py
│  ├─ VSPAEROStab.py
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

VSPAERO surface tessellation の初期化と mesh-convergence search を担当します。

- `equalize_vspaero_tessellation()` — 実3次元 edge length を基準に、`Tess_W` / `SectTess_U` / `Tess_U` / active cap tessellation を揃える初期化処理
- `resolve_vspaero_representation()` — `ThinWing / ThickWing / Hybrid / ThickAll` をモデル内の Geom Set に解決
- `optimize_vspaero_tessellation()` — geometry・representation・clustering・wake 条件を固定し、W方向→U方向の順に mesh convergence を探索

`optimize_vspaero_tessellation()` は Boeing や G103A の Geom 名を前提にしません。`representation` は必須で、既定では次の Set 名を使います。

| representation | ThinGeomSet | GeomSet |
|---|---|---|
| `ThinWing` | `ThinGeom` | none |
| `ThickWing` | none | `ThinGeom` |
| `Hybrid` | `ThinGeom` | `ThickGeom` |
| `ThickAll` | none | `ThickAll` |

Set 名が異なるモデルでは `lifting_set_name` / `body_set_name` / `thick_all_set_name` を変更できます。`mesh_targets=None` なら、選択した representation で active な Geom のうち `Tess_W` を持つ surface Geom を自動的に対象にします。特定 Geom に限定する場合は Geom 名または Geom ID を `mesh_targets` に指定します。各 case の `.vsp3` には選択した `GeomSet / ThinGeomSet` も保存するため、出力モデル自体にも optimization representation が残ります。

第一版の optimizer は、clustering や geometry を同時最適化しません。まず `Tess_W` を倍率的に変化させ、その selected case から `SectTess_U / Tess_U` を倍率的に変化させます。各 case は `VSPAEROMeshQuality.py` で診断し、次を分離して扱います。

- hard constraint: mapping / non-manifold / Kutta などの geometry-topology check
- primary convergence: ユーザー指定 QoI の mesh sensitivity
- regression guard: strong mesh advisory / local Cp spike / LOD outlier
- cost / provenance: NGon / triangle 数、wall time、OpenVSP version、各 raw artifact

使用例:

```python
from src.VSPAEROMesh import optimize_vspaero_tessellation

result = optimize_vspaero_tessellation(
    input_vsp3_path="examples/models/G103A/G103A.vsp3",
    output_dir="results/mesh_optimization",
    representation="Hybrid",
    qoi_tolerances={
        "CLiw": 0.005,
        "CDiw": 0.01,
        "CMytot": 0.01,
    },
    # near-zero quantity の正規化 scale が必要なら明示する
    qoi_scales={"CMytot": 0.05},
    fixed_wake_flag=True,
)

print(result["summary"]["fully_converged"])
print(result["selected_vsp3_path"])
```

`qoi_tolerances` は普遍的な OpenVSP 公式閾値ではないため、解析目的に応じて利用者が指定します。隣接 mesh 間の判定には、各 QoI について

```text
abs(q_fine - q_coarse) / max(abs(q_fine), abs(q_coarse), qoi_scale)
```

を使用します。ゼロ近傍の moment 等では `qoi_scales` を明示してください。optimizer は `CLi-CLiw` の一致自体を目的関数にはしません。

### `src/VSPAEROMeshDataset.py`

将来の形状へ転用できる mesh sizing rule を作るための **raw calibration dataset generator** です。optimizer ではありません。`ThinWing / ThickWing / Hybrid / ThickAll` の active Geom を読み、実際に存在する mesh Parm だけを whitelist から発見します。

対象は `Tess_W` / `Tess_U` / `SectTess_U` / `CapUMinTess` / `LECluster` / `TECluster` / `InCluster` / `OutCluster` / `FwdCluster` / `AftCluster` です。Geom type ごとに固定の section 数や名前を仮定しません。各 Parm は独立自由度として catalog 化され、既定では全自由度の one-factor sweep と、全自由度を同時に振る Latin-hypercube sample を生成します。必要なら pairwise extreme block や明示 case を後から同じ dataset に追加できます。

同じ `output_dir` への再実行は append-only です。source model hash、OpenVSP version、representation、飛行条件、solver/wake 条件、全 requested mesh Parm から deterministic `case_id` を作り、完了済み case は再計算しません。別モデルや別 representation も同じ dataset に追加できます。

永続データの source of truth は `cases/<case_id>/attempt_xx/` です。解析開始時に `request.json` を保存し、`.vsp3`、VSPAERO artifact、`mesh_quality/`、各 raw CSV を保存した最後にだけ `case.json` を atomic に作成します。`case.json` がない attempt は途中中断として扱われ、次回は新しい attempt で再実行されます。したがって集約 CSV は resume 判定には使用しません。

各 attempt に保存する raw data は `parameters.csv` / `geoms.csv` / `sections.csv` / `junctions.csv` / `polar.csv` です。dataset root の `cases.csv` などは派生集計で、正常終了時に一度だけ再構築されます。途中終了後でも `rebuild_mesh_dataset_tables(output_dir)` を呼べば、commit 済み `case.json` だけから再生成できます。

VSPAERO sweep の `verbose` は dataset campaign では常に 0 です。campaign 自身の `verbose` は進捗表示だけを制御し、現在時刻、case runtime、累積時間、直近20件の正常終了 case の runtime 中央値を使った estimated finish time を表示します。

生成される主な集約表は次のとおりです。

- `cases.csv` — 実行状態、開始/終了時刻、計算時間、mesh/Cp/LOD/topology の要約、artifact hash
- `parameter_catalog.csv` — source model に存在した全 mesh Parm と baseline / limits
- `parameters.csv` — 各 case の requested / effective Parm 値
- `geoms.csv` — geometry scale、曲率、実3D U/W edge distribution、adjacent growth、`small_panel_w` / `max_growth_w`
- `sections.csv` — section/cap の物理長、実U-edge、XSec寸法、section clustering
- `junctions.csv` — post-intersection junction の raw diagnostics
- `polar.csv` — VSPAERO polar result を long format で保存

`SmallPanelW` / `MaxGrowth` は入力自由度にはせず、生成後の実3D edge から dataset 側で再計算します。また dataset generator は `good_mesh` のような教師ラベルを作りません。topology family、convergence plateau、reference candidate、許容誤差などは raw data を保持したまま後段で再定義します。

使用例:

```python
from src.VSPAEROMeshDataset import (
    build_vspaero_mesh_rule_dataset,
    rebuild_mesh_dataset_tables,
)

result = build_vspaero_mesh_rule_dataset(
    input_vsp3_path="examples/models/Boeing_777-9x_mod/Boeing_777-9x_mod.vsp3",
    output_dir="results/mesh_rule_dataset",
    representation="ThickAll",
    lhs_samples=256,
    ncpu=8,
    wake_num_iter=12,
)

print(result["planned_case_count"])
print(result["parameter_catalog_path"])

# 途中停止後でも、commit 済み case だけから集約表を再構築できます。
rebuild_mesh_dataset_tables("results/mesh_rule_dataset")
```

長時間 campaign を拡張するときは、同じ `output_dir` を指定したまま `lhs_seed` / `lhs_samples` を変更する、`include_pairwise_extremes=True` を追加する、または `explicit_cases` を渡します。

### `src/VSPAEROStab.py`

VSPAERO `.stab` ファイルの読み取りを担当します。

- reference quantities
- base flight condition
- base aerodynamic coefficients
- stability derivatives
- Control Surface Group と `ConGrp_*` の対応

を共通形式へ変換し、定常旋回 solver、6DoF simulation、設計チャートから共通利用します。

### `src/VSPAEROMeshQuality.py`

VSPAERO の `.adb` v3 を読み、surface triangle の幾何品質と `Cp` を同じ ID 上で診断します。`.vspgeom v3` がある場合は alternate triangulation を original NGon に対応付け、NGon 単位の局所 `Cp`、Kutta / wake topology、実際に生成された 3D edge 長、junction の cut-edge 寸法を同時に確認できます。

`.history` / `.polar` / `.lod` / `.vspaero` / `.vsp3` がある場合は solution-level diagnostics と解析 provenance も同じ report に含めます。入力した `Tess_W` や `SectTess_U` の値だけでは mesh quality を判断せず、`.vspgeom` に実際に書かれた post-intersection mesh の物理寸法を出力します。

主な公開関数は `analyze_vspaero_mesh_quality()` です。

```python
from VSPAEROMeshQuality import MeshQualitySettings, analyze_vspaero_mesh_quality

result = analyze_vspaero_mesh_quality(
    "model.adb",
    "model.vspgeom",
    "mesh_quality",
    history_path="model.history",
    lod_path="model.lod",
    vspaero_path="model.vspaero",
    vsp3_path="model.vsp3",
    settings=MeshQualitySettings(),
)
```

v5 では API を整理しています。

- threshold 群は `MeshQualitySettings` に集約
- `component_ids` / `surface_ids` / `bbox` による main analyzer 内の部分 filter は削除
- `summary` は `provenance` / `condition` / `mesh` / `cp` / `topology` / `kutta` / `solution` / `lod` / `checks` に整理
- `structural_checks_passed` は廃止し、`summary["checks"]["geometry_topology_checks_passed"]` を使用
- v3 compatibility alias は削除

局所的な抽出は analyzer が返す `triangles` / `ngons` / `mesh_edges` / `junction_quality` を呼び出し側で filter してください。これにより、一つの summary 内で「一部領域の Cp」と「全機の Kutta topology」が混在することを避けています。

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
- control surface gain と舵角符号
- VSPAERO analysis method
- wake settings

G103A の trimmed-polar workflow では、既存モデルの `WingGeom` と `ELEVATOR_GROUP` を使用します。

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
