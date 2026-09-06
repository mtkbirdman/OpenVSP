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

### OpenVSP を使わずに定常旋回 solver を確認する

`G103A.stab` は定常旋回例の入力として同梱しています。OpenVSP / VSPAERO を再実行しなくても、NumPy・pandas・SciPy があれば後処理を確認できます。

```bash
python examples/scripts/test_turn_trim.py
```

この script では、既存の `TrimTurnSolver.py` を用いて、代表的な次のケースを計算します。

1. 定常滑空・協調旋回
2. 定常高度維持・ラダーのみ旋回
3. 定常滑空・ラダー限界旋回

### VSPAERO sweep を実行する

OpenVSP Python API が使用できる環境では、次のように G103A の sweep を実行できます。

```bash
python examples/scripts/test_sweep.py
```

同様に、次の script があります。

```text
examples/scripts/test_control_surface.py
examples/scripts/test_ground_effect.py
examples/scripts/test_stability_derivatives.py
examples/scripts/test_trimmed_polar.py
examples/scripts/test_turn_trim.py
```

各 script は `__file__` からリポジトリルートとモデル位置を求めるため、リポジトリルートから直接実行できます。

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

Notebook は、基本的にリポジトリ内の `src/` と `examples/models/` を利用する構成です。

## ディレクトリ構成

```text
OpenVSP/
├─ README.md
├─ src/
│  ├─ AnalysisVSPAERO.py
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
   ├─ scripts/
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
