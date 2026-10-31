# OpenVSP / VSPAERO Python解析ツール

OpenVSP / VSPAEROをPythonから利用し、航空機の空力解析・飛行力学解析・可視化と、VSPAEROのVerification & Validation（V&V）を支援するリポジトリです。解析コード、Notebook、再現可能なV&V定義、試験用データ、技術文書を収録しています。

> **対象と制約**：このリポジトリはOpenVSP本体の配布物ではありません。OpenVSP本体、互換性のあるPython API、VSPAERO実行ファイルは別途用意してください。V&V基盤を実装していることは、VSPAEROの予測精度や工学的使用可能性が認定されたことを意味しません。

## 1. 主な機能

| 分野 | 機能 | 主な実装 |
|---|---|---|
| VSPAERO解析 | 迎角・Mach数・Reynolds数のSweep、縦トリム、地面効果、安定微係数 | `src/openvsppy/vspaero/` |
| 空力・飛行力学 | 寄生抗力補正、定常旋回トリム、ラダー応答、横風突風の6自由度応答 | `src/openvsppy/aerodynamics/` |
| メッシュ | 3次元エッジ長に基づくVSPAEROメッシュ分割の調整 | `src/openvsppy/vspaero/VSPAEROMesh.py` |
| 可視化 | `.vsp3`の複数方向撮影、PNG合成、GIF・MP4アニメーション | `src/openvsppy/visualization/` |
| V&V | 試験定義、実行管理、SQLite記録、結果比較・評価、図表・報告書生成 | `src/openvsppy/vv/` |

解析モデル・Notebook・設定例は`examples/`、V&Vの仕様とテンプレートは`docs/vv/`、実際のStudy設定と実行時の参照データは`vv_study/`で管理します。回帰テストに必要な固定データは`tests/fixtures/`に配置します。

## 2. 動作環境

- **Python**：`pyproject.toml`の指定は3.10以上。ただし、OpenVSP Python APIを使用する場合は、**インストールされたOpenVSP 3.50.5 APIのビルド対象Pythonと互換性のある環境**が必要です。Pythonのバージョンだけでは共有ライブラリの動作まで保証されません。
- **OpenVSP / VSPAERO**：V&Vのネイティブ試験は3.50.5を対象としています。オフラインテストの一部はOpenVSPなしで実行できます。
- **主なPythonライブラリ**：NumPy、pandas、SciPy、Matplotlib。撮影・動画関連機能にはPillow、imageio-ffmpeg、NotebookにはJupyterが必要です。V&V報告書の追加依存は`requirements-vv-report.txt`を参照してください。
- **実行環境**：以下はWindows PowerShellの例です。パスとConda環境名は利用環境に合わせてください。

公式情報：[OpenVSP](https://openvsp.org/)／[Python APIリファレンス](https://openvsp.org/pyapi_docs/latest/)／[ドキュメント](https://openvsp.org/docs.shtml)

### 既存リポジトリからの起動

```powershell
# 自分のリポジトリルートへ移動（実際のパスへ変更）
$repo = 'C:\path\to\OpenVSP'
Set-Location $repo
conda activate vsppytools.3.50.5
$env:PYTHONPATH = (Resolve-Path '.\src').Path

python -m openvsppy.vv --help
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
├── requirements-vv-report.txt
├── src/openvsppy/
│   ├── core/                   # 共通ユーティリティ・大気モデル
│   ├── geometry/               # 形状関数
│   ├── vspaero/                # 解析API・メッシュ・安定微係数
│   ├── aerodynamics/           # 飛行力学・抗力補正・設計計算
│   ├── visualization/          # OpenVSPの撮影・画像処理
│   └── vv/                     # V&V CLI・DB・ワーカー・報告
├── examples/
│   ├── models/                 # 解析モデルと付属データ
│   ├── notebooks/              # 解析・可視化の実行例
│   └── vv/                     # V&V設定・参照データの使用例
├── tests/                      # オフライン・ネイティブ・GUIテスト
│   └── fixtures/               # テストに必要な固定データ
├── tools/
│   └── run_native_pytest.py    # ネイティブテストの隔離実行
├── docs/vv/                    # V&V計画・仕様・運用・履歴
│   └── templates/              # 報告書・試験記録・設定・調査表のひな形
└── vv_study/                   # 現行Study設定・試験カタログ・実行時データ
```

`.gitignore`の対象となるSQLite DB、解析Run、生成レポート、Notebook結果などは、通常のGitコミットには含めません。回帰テストに必要な選別済みfixtureのみを収録し、一時的な比較調査の原本は別アーカイブに保管します。

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

## 5. VSPAERO V&Vの実行

V&Vには、計画書、試験カタログ、レビュー済み入力設定、SQLite、解析結果の証拠ファイルを利用します。解析条件が未設定のCaseや、比較に必要な実験データがないCaseを推測で補完しません。

### CLIの確認と新規データベース

```powershell
# リポジトリルートにおいて、上記の環境設定後に実行
$db = 'vv_study\01_database\vv.sqlite'
python -m openvsppy.vv --help
python -m openvsppy.vv --db $db init
python -m openvsppy.vv --db $db plan
```

`init`は新規DBの作成に使います。旧版SQLite DBの移行機能はありません。旧`vv`モジュールのインポートパスにも互換性はありません。**既存DBを直接上書きしたり、試験Runの証拠ファイルを削除したりしないでください。**

### ネイティブ解析と一括実行

OpenVSP 3.50.5のAPIとVSPAERO実行ファイルが使用可能な環境では、個別Caseの確認後に一括処理を実施できます。以下はコマンド例です。

```powershell
$solver = 'C:\OpenVSP-3.50.5\vspaero.exe'  # 実際のパスに変更
$pythonExe = (Get-Command python).Source
$config = 'vv_study\execution.reviewed.json'
python -m openvsppy.vv --db $db check --case GEO-03-G3 --binary $solver
python -m openvsppy.vv --db $db execute-all --config $config --binary $solver --python $pythonExe --timeout 3600 --ncpu 4 --verbose
```

個別のCaseは事前にレビュー済み条件で`prepare`・登録する必要があります。`--ncpu`はVSPAEROソルバーのOpenMPスレッド数であり、同時Case数ではありません。詳細な手順、解析モード別タイムアウト、Runの再利用条件、結果確認は[CLI運用ガイド](docs/vv/cli_guide.md)を参照してください。

**`execute-all`による報告書生成やRunの`SUCCEEDED`登録は、V&Vの合格や工学的使用認定を意味しません。** 実機での検証状況・未確認事項は[実装状況](docs/vv/implementation_status.md)と[変更履歴](docs/vv/CHANGELOG.md)に記載しています。

## 6. テスト

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

## 7. ドキュメント

| 文書 | 用途 |
|---|---|
| [V&V計画書](docs/vv/plan.md) | 試験目的・範囲・技術的判定方針。文書版数は本文とStudy設定で管理 |
| [DB詳細仕様書](docs/vv/database_spec.md) | データモデル・状態遷移・証拠管理 |
| [評価報告書テンプレート](docs/vv/templates/evaluation_report.md) | 13章構成のレビュー草稿 |
| [試験記録テンプレート](docs/vv/templates/test_record.md) | Test ID単位の記録 |
| [DB・報告書フィールド対応表](docs/vv/db_report_field_mapping.md) | DBと報告書の追跡 |
| [方法論と参考文献](docs/vv/methods_and_references.md) | 採用した手法と一次資料の対応 |
| [CLI運用ガイド](docs/vv/cli_guide.md) | 最新の実行手順 |
| [実装状況](docs/vv/implementation_status.md) | 25試験の実装範囲と未検証項目 |
| [変更履歴](docs/vv/CHANGELOG.md) | 過去の開発記録・版別の変更点 |
| [設定・報告書テンプレート](docs/vv/templates/README.md) | 新規Study・Case・参照メタデータ、報告書のひな形 |
| [V&V設定の使用例](examples/vv/) | 具体的なG103A解析条件・参照値の例（正式承認とは別） |

計画書およびテンプレートの**ファイル名は固定**します。計画書の文書版数、Studyの`plan.version`・`plan.sha256`、特定Gitコミットの組合せで、実際に使用した計画の内容を識別できます。計画書本文を改訂した場合は対応するSHA-256を更新し、過去のRunに使用した計画の特定情報を保持してください。

## 8. 注意事項・ライセンス

- 解析値の妥当性は、形状・解析モード・座標系・基準量・重心・Control Surface Group、メッシュ・Wake条件、参照データおよび数値収束に依存します。
- `.stab`による線形モデルの外挿や、VSPAERO内蔵の寄生抗力と外部CD0補正の組合せには適用範囲上の制約があります。
- OpenVSP、Python API、VSPAEROの仕様はバージョンで異なる場合があります。
- このリポジトリには現時点でルートの`LICENSE`ファイルがありません。権利・再配布条件は別途確認してください。第三者由来のモデル、フライトマニュアル、エアフォイル等には個別の利用条件が適用される場合があります。
