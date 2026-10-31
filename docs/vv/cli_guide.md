# VSPAERO V&V — CLI運用ガイド

本書は**現在有効な操作手順**です。旧版の`python -m vv`、`src.VV*`、`vv_study/05_report_templates/`、SQLite `migrate`などの手順は廃止済みで、過去のリリースノートからは実行しないでください。実装の真の引数仕様は`src/openvsppy/vv/__main__.py`（`--help`）です。

## 1. 前提と配置

- Pythonパッケージ：`src/openvsppy/`、V&V CLI：`python -m openvsppy.vv`、APIワーカー：`python -m openvsppy.vv.worker`。
- Study設定・承認済み入力：`vv_study/`、DDLとCatalog：`vv_study/01_database/`。
- 正式文書・雛形：`docs/vv/`。**13章評価報告書テンプレートの章構成と見出しを維持する**。
- 生成されるSQLite、Run、図表、レポートはStudy側に保存し、Git管理対象外とする。
- OpenVSP Python APIとVSPAERO実行ファイルは別途必要。`python -m pip install -e .`は**任意**。Notebookは先頭セルの`sys.path`追加のみでソースをimportする。

## 2. Windows PowerShell — 新規DBからGEO-03-G3まで

新規展開したリポジトリ内で実行する例。ZIP展開先、Conda環境、solverパスは実機に合わせて変更してください。過去データとの互換性を保証するための移行は行いません。必要なら新しいDBで一連の試験を実行してください。

```powershell
Set-Location 'C:\path\to\OpenVSP'
conda activate vsppytools.3.50.5
$env:PYTHONPATH = (Resolve-Path '.\src').Path
$pythonExe = (Get-Command python).Source
$solver = 'C:\path\to\OpenVSP-3.50.5\vspaero.exe'
$db = 'vv_study\01_database\vv_20261031_11.sqlite'
$config = 'vv_study\execution.reviewed.json'

python -c "import openvsp as vsp; print(vsp.GetVSPVersion())"
Get-FileHash -LiteralPath $solver -Algorithm SHA256
python -m pytest -q -rs
if ($LASTEXITCODE -ne 0) { throw 'Offline tests failed' }
python -m openvsppy.vv --db $db init
if ($LASTEXITCODE -ne 0) { throw 'SQLite init failed (fresh database required)' }
python -m openvsppy.vv --db $db plan
python -m openvsppy.vv --db $db prepare --test GEO-03 --base vv_study\g103a_baseline.reviewed.json --register
python -m openvsppy.vv --db $db check --case GEO-03-G3 --binary $solver
python -m openvsppy.vv --db $db run --case GEO-03-G3 --binary $solver --python $pythonExe --timeout 3600 --ncpu 4
if ($LASTEXITCODE -ne 0) { throw 'Native GEO-03-G3 failed: inspect Run logs' }
```

**注意：** `init`は既存DBの初期化コマンドではありません。既存DBとのMigrationはありません。解析結果はワーカーがOpenVSP Python APIの`ExecAnalysis('VSPAEROComputeGeometry')`と`ExecAnalysis('VSPAEROSweep')`で生成します。VSPAEROへの直接起動フォールバックはありません。

## 3. 一括実行、集計、照合

GEO-03-G3の実機回帰結果を確認後に実施します。`--no-resume`で成功Runの再利用を抑制できます。目標はすべての実行可能Caseを解析・後処理・報告まで一気通貫で進めることであり、Runの蓄積自体は目的ではありません。`--no-run`は実解析を行わない報告書生成・動作確認用です。

```powershell
python -m openvsppy.vv --db $db execute-all --config $config --binary $solver --python $pythonExe --timeout 3600 --ncpu 4 --verbose
if ($LASTEXITCODE -ne 0) { throw 'execute-all workflow failed' }
python -m openvsppy.vv --db $db verify-artifacts
python -m openvsppy.vv --db $db assess-all --output vv_study\04_reports\evidence_inventory.json
python -m openvsppy.vv --db $db export --output vv_study\04_reports\manual_export
python -m openvsppy.vv --db $db report --output vv_study\04_reports\manual_review.md
```

`execute-all`が正常終了しても、Caseごとに`BLOCKED`・欠測・解析失敗があり得ます。**`REVIEW_REPORT_GENERATED`は全解析成功やV&V合格を意味しません**。レビュー済みマニフェストに指定されていない参照値・飛行条件は自動補完しません。

## 4. 個別調査・再処理

```powershell
python -m openvsppy.vv --db $db inspect-run --run <RUN_ID>
python -m openvsppy.vv --db $db reprocess-run --run <RUN_ID>
python -m openvsppy.vv --db $db run-test --test WAKE-01 --binary $solver --python $pythonExe --resume --jobs 1 --ncpu 4
python -m openvsppy.vv --db $db evaluate --test GRID-01
python -m openvsppy.vv --db $db gci --values 0.1 0.11 0.12 --widths 0.5 1.0 2.0
```

`<RUN_ID>`等は**実際の識別子へ置換**してから実行してください。GCIの数値はCLI構文例であり、G103Aの承認済み格子幅や収束判定ではありません。`reprocess-run`は元ネイティブファイルを変更せず派生データを追加します。

参照Datasetの登録・比較（`import-reference`、`compare`、`compare-test`）と抗力補正（`drag-run`、`drag`）は、レビュー済みメタデータ・座標系・正規化条件の一致を必要とします。具体的な引数は`python -m openvsppy.vv <command> --help`を参照してください。暗黙の軸変換・検証未済の物理量認定はしません。

## 5. VSPAERO NCPUの指定方法

`NCPU`は**各VSPAEROソルバーのOpenMPスレッド数**です。Studyの`runtime.processes=1`は同時に起動する解析ワーカープロセス数であり、NCPUとは異なります。`--jobs`も同時実行のためのパラメータです（現在は1のみ対応）。

**優先順位**：`run`／`run-test`／`execute-all`の`--ncpu`、Case JSONの`numerical_setting.requested.vspaero.NCPU`、省略時の従来値`1`の順に採用します。値は正の整数のみ受け付けます。通常は次のように指定します。

```powershell
python -m openvsppy.vv --db $db run --case GEO-03-G3 --binary $solver --python $pythonExe --ncpu 4
python -m openvsppy.vv --db $db run-test --test WAKE-01 --binary $solver --python $pythonExe --resume --ncpu 4
python -m openvsppy.vv --db $db execute-all --config $config --binary $solver --python $pythonExe --timeout 3600 --ncpu 4 --verbose
```

Caseごとに指定したい場合は、**Case登録前**のCase JSONで次の項目を追加します。すでに登録されたCaseのNumericalSetting/setting_idは不変です。既存の承認済み定義を書き換えず、必要なら新しい設定ID・新規DBを使用してください。

```json
{
  "numerical_setting": {
    "requested": {
      "vspaero": {"NCPU": 4}
    }
  }
}
```

これは**既存Case JSONに追加する断片例**であり、単独で登録可能な完全なCase定義ではありません。CLI `--ncpu`は登録済みCaseを書き換えず、実行時だけスレッド数を指定します。

- `worker.py`は`VSPAEROSweep`の整数入力`NCPU`に値を設定し、API読戻しを検証します。
- 実効値を`raw/api_diagnostics.json`の`ncpu`・期待する`-omp`引数、`worker_result.json`の`effective.sweep.NCPU`、`input/job.json`へ記録します。
- SQLiteの`run.environment_json.ncpu`と実行指紋に反映し、Runの再利用ではNCPU一致を必要条件とします。**既存NCPU=1の成功RunをNCPU=4として再利用しません**。
- 過去版のRunには`environment_json.ncpu`が記録されていないため、最新版では保守的に新しいRunを作成します。DBスキーマ変更や過去Runの書換えは行いません。
- 同じNCPUでも、解析ソース指紋・VSPAEROバイナリSHA-256・必須Artifactの照合など、従来の再利用条件を満たす必要があります。
- `--ncpu 4`は4スレッドの**要求**であり、実際のCPU占有率や処理時間の保証ではありません。Windows実機のv20261031.11 API経由解析・数値同一性は未検証です。

## 6. 生成ファイルと診断

- `vv_study/01_database/*.sqlite`：新しいStudy DB。
- `vv_study/03_run_raw/RUN-.../{input,logs,raw}`：入力・ワーカーログ・ネイティブ結果・API診断JSON。
- `vv_study/04_reports/execute_.../`：`progress.jsonl`、固定DBスナップショット、CSV、図表、13章の`evaluation_report.md`。
- `raw/vspaero_api_stdout.log`：C++側監視出力。`raw/api_diagnostics.json`：入力ハッシュ、APIモード、Run実行環境など。
- `inspect-run`と`verify-artifacts`はファイル・レコード整合性を検査するもので、ソルバー数値解や工学的妥当性を自動保証しません。

## 7. Notebook

全13件で先頭セルからリポジトリの`src`を`sys.path`へ追加します。Editable Installは不要です。ただしOpenVSP API、Notebook固有のライブラリや入力データは別途必要です。`example_turn_trim.ipynb`が参照する`G103A.stab`はZIPに含まれません。Notebookごとの結果保存先は`examples/notebooks/`以下の管理対象外ディレクトリを使用します。

## 8. 文書体系

- 正本：`plan.md`、`database_spec.md`
- 固定雛形：`templates/evaluation_report.md`、`templates/test_record.md`
- 付属：`db_report_field_mapping.md`、`methods_and_references.md`
- 運用：**本書**、`implementation_status.md`、`CHANGELOG.md`

## 9. 受入と未検証範囲

v20261031.3のWindows実機でGEO-03-G3が20/20 Wake反復、`.polar/.history/.lod/.adb`、SQLite `SUCCEEDED`に到達した記録があります。**v20261031.4以降のWindowsネイティブ実機試験、STA-01～03のネイティブ完走は未確認**です。履歴は`CHANGELOG.md`を参照してください。

### execute-all の実行進捗（v20261031.9）

標準ではCaseの確認・起動・終了（`SUCCEEDED`、`REUSED`、`BLOCKED`など）、処理時間、フェーズ開始・終了を表示します。60秒間ワーカーが終了しない場合は経過秒数の生存確認を出します。`--verbose`では、保存済みの`raw/vspaero_api_stdout.log`を読み取り、観測できたWake反復の開始番号・Run ID・必須出力ファイルサイズ・SQLite登録を追加表示します。反復表示は**反復開始の観測**であり、反復完了・数値収束の認定ではありません。ログ出力が遅延した場合は反復表示が省略または遅れる場合があります。

```powershell
python -m openvsppy.vv --db $db execute-all --config $config --binary $solver --python $pythonExe --timeout 3600 --ncpu 4 --verbose
```

画面進捗は**標準エラー出力（stderr）**、最終JSONは**標準出力（stdout）**へ分離します。全イベントを`vv_study/04_reports/execute_*/progress.jsonl`に逐次記録します。Windows PowerShellには両方が表示されます。既存のRun成功判定・ソルバー設定・`--timeout`は変更しません。`--no-run`は解析を実行せず、報告書生成のフェーズ進捗のみ表示します。結果が`REVIEW_REPORT_GENERATED`でも全Caseの成功を意味しません。

## 10. pytestの通常実行とNative実行

v20261031.9から、**通常の `python -m pytest -q -rs` はOpenVSP API・VSPAERO・GUIを起動しません**。`pyproject.toml` にマーカーを集約し、旧 `pytest.ini` および `RUN_VSPAERO_INTEGRATION` / `RUN_OPENVSP_GUI_TESTS` のオプトイン方式を廃止しました。通常の回帰でネイティブケースを含めないのは意図的な `DESELECTED` であり、API/Solverの実証結果ではありません。

| マーカー | 対象 | 現時点の個別ケース数 |
|---|---|---:|
| `native_openvsp` | OpenVSP Python APIのみを実使用 | 4 |
| `native_vspaero` | OpenVSP API経由でVSPAEROを起動 | 7 |
| `native_solver` | VSPAERO単体の診断スモーク | 1 |
| `gui` | GUI/OpenGLを使用する画像キャプチャ | 1 |
| `slow` | 長時間に及ぶ可能性を付記（選択区分とは独立） | 複数 |

通常の`pytest`でNativeを意図せず実行しない一方、Nativeを明示選択したときは、APIのimport・バージョン、実行ファイル、`SetVSPAEROPath`のエラーを**FAILとして報告**します（曖昧なSKIPにはしません）。Native APIのバージョンは3.50.5を必要とします。

```powershell
# PowerShell; リポジトリルート、正しいConda/OpenVSP環境で実行する
$env:VSPAERO_BINARY = 'C:\OpenVSP-3.50.5\vspaero.exe'  # 実際の配置へ変更
python -m pytest -q -rs
python tools/run_native_pytest.py --marker all-native --list-only
python tools/run_native_pytest.py --marker native_openvsp --timeout 600
python tools/run_native_pytest.py --marker native_vspaero --timeout 3600
python tools/run_native_pytest.py --marker native_solver --timeout 240
# GUI表示環境でのみ実行する
python tools/run_native_pytest.py --marker gui --timeout 300
```

`tools/run_native_pytest.py`は、収集した各テストを**別のpytest子プロセス**で実行します。1 Caseあたりの制限時間と60秒Heartbeatを設け、タイムアウト時にはWindowsの`taskkill /T /F`、Unixのプロセスグループ終了によってソルバー子プロセスも停止します。タイムアウト・例外・非ゼロ終了はPASSとして扱いません。`--nodeid`で単独Caseを選択でき、`--heartbeat`と`--report-dir`も指定できます。ログと各Caseの一時作業領域、`summary.json`は`native_test_results/`内に保存して再検証可能です。通常の直接`python -m pytest -m native_vspaero`ではこのプロセス単位の制限時間は働かないため、長時間試験にはランナーを推奨します。

Native試験は`isolated_model`・`native_workspace` fixtureを用いてモデル・補助入力をコピーし、作業ディレクトリと結果ファイルをリポジトリ内の承認済みファイル・Run証拠から隔離します。数値解析関数、G103A Gain、Wake、SQLite Schemaは変更しません。

**`native_solver`の診断モデルは正式なGEO-03-G3ではありません。** 添付Ubuntu版VSPAEROで成功したG103Aの診断コピーから、Control Group数0、Wake反復1の`.vspgeom/.vspaero/.vkey/.csf`を`tests/fixtures/native_solver_smoke/`に収録しています。ソルバーを独立起動して `.polar/.history/.lod/.adb`を検査し、`.history`を既存Parserで再読込します。これはAPI経由の20 Wake解析、25 Test IDの合格判定、数値解の収束性を証明しません。LinuxのOpenVSP 3.50.5 Python APIにはPython 3.12と`libcminpack.so.1`等の共有ライブラリが必要です。依存関係不足の環境ではAPI試験をPASSとしません。


## 11. 失敗診断・Adjoint・Case別タイムアウト

`execute-all`は各Caseに対して検証可能な状態を残し、正常終了済みのCaseだけでなく、失敗したCaseも報告書に記録します。V&Vの工学的合否は参照データ、数値収束、不確かさを検討した後に人間が判断します。

- **Geometry SetとControl Group：** ワーカーは実際のThin/Thick Setに含まれる形状からControl Surfaceの親形状を照合します。選択外の形状に属するControl Groupは偏向角0の場合に無効化します。Group内で一部の舵面だけが欠落、親形状が不明、または欠落Groupに非ゼロ偏向を要求した場合は、曖昧な入力をソルバーへ渡さず`EXECUTION_FAILED`とします。Gainと物理的符号は変更しません。
- **VSPGEOMとVKEY：** `VSPAEROComputeGeometry`で生成された実ファイルを検査し、極端に少ないLoop、空の格子、VKEYの表面部品数不正を`mesh_diagnostics.json`へ記録します。構造検査はGCIの合格判定ではありません。`NRef`は改変せず、異常値を使ったソルバー実行や自動補正は行いません。
- **ソルバー失敗：** `api_diagnostics.json`に最終観測ステージとログ由来の分類を保存します。OpenVSP Python APIが子ソルバーの終了コードを公開しないため、取得不能な値を推測しません。VSPAEROが出力を生成しない場合、`SUCCEEDED`として登録しません。
- **Adjoint進捗：** `--verbose`時にWakeと併せ、`Solving adjoint linear systems`、力・モーメント成分別の進行、およそ10反復単位のGMRES進捗を表示します（ネイティブログが利用できる場合）。進捗表示の欠落は単独では解析の成否を決定しません。
- **解析モード・Case別タイムアウト：** `vv_study/execution.reviewed.json`内に以下のキーを設定します。`timeouts_seconds_by_case`が最優先、次に`timeouts_seconds_by_mode`、最後にCLIの`--timeout`を適用します。値は秒単位の正の整数です。

```json
{
  "timeouts_seconds_by_mode": {"ADJOINT": 10800},
  "timeouts_seconds_by_case": {"STA-01-ADJOINT": 14400}
}
```

上記は設定キーの追加例であり、**既存の`execution.reviewed.json`の他の項目を消さない**でください。配布物の既定値は`ADJOINT=10800`、Case別設定は空です。Wake反復数・NCPU・GMRES閾値等の数値解析条件は変更しません。長時間のAdjoint実行が必ず収束する保証はありません。

**実機での受入：** v20261031.10のWindowsログではGEO-03-G3がNCPU=8で成功、NCPU=4の一括実行では複数Caseが成功する一方、GEO-01/GEO-02のControl Group不整合、`GRID-03-NREF_1`のメッシュ関連異常、Adjoint3600秒タイムアウト等が観測されました。v20261031.11の修正による成功は、Windows実機で再実行した後に初めて確認できます。
