# VSPAERO 3.50.5 V&Vデータベース作成プログラム 詳細仕様書

| 項目 | 値 |
|---|---|
| 文書番号 | VSP-VV-SW-SPEC-001 |
| 版数 | 1.0.1（上位計画書v1.1への参照追従。機能変更なし） |
| 作成日 | 2026-10-21 |
| 対象 | OpenVSP/VSPAERO 3.50.5、複雑な全機形状（胴体・主翼・水平尾翼・垂直尾翼） |
| 上位計画 | `plan.md`、SHA-256は`templates/study.json`に格納 |
| 既存コード | `OpenVSP-main_20260927_trimmed_sweep_activeflag_v1.zip` 内 `src/` と `tests/` |
| 仕様付属品 | `vv_study/01_database/schema.sql`、`templates/study.json`、`vv_study/01_database/test_catalog.csv` |
| 成果物 | SQLite正本、実行・参照の生データ、Excel/CSVレビュー出力、比較図および判定記録 |

> **本文の位置付け**：これはプログラムの実装仕様であってVSPAEROのValidation済み証明ではない。`schema.sql`はSQLiteで構文検査済み、実行ライフサイクル・OpenVSP API連携は未実装・未実機検証である。計画書で未決の飛行条件、基準面積、審査者、要求精度は勝手に補完しない。

## 1. 採用する設計判断

1. **SQLiteを唯一の編集可能な正本**とする。Excel/CSVは参照・配布用の生成物であり、Excelからの編集結果を暗黙にDBへ逆輸入しない。
2. **試験（Test ID）／予定ケース（Case ID）／解析実行（Run ID）／出力条件点（Point ID）を独立させる**。1回のRunの出力点を複数テストが参照できる。
3. **1 Run = 1明示的な飛行条件点**を初期版の運用原則とする。既存の等間隔VSPAERO Sweepは非等間隔指定をそのまま再現しないため、単点実行を基本にする。ただし既存結果の複数点取り込みに備え`case_point`を設ける。
4. **生出力を不可変にする**。再計算した場合は必ず新規Run ID。`.polar`、`.stab`等のネイティブデータ、解析設定、モデル、ログのハッシュを保持。
5. **実行成功、数値受容、物理Validation、用途認定を別々に管理する**。
6. **元のVSPAERO数値とCD0補正・微係数の差分算出などの派生値を別Datasetとして保持する**。内蔵`CDo`と外部`CD0`を二重計上しない。
7. 初期版は**Windows上の単一プロセス・逐次実行**。OpenVSP APIに保持されるモデル状態をRun間で引き継がない。必要になったらプロセス単位の並列化を追加する。
8. GUI、新規の汎用CFDファイルリーダー、全自動GCI/不確かさ解析の完成を初期必須条件としない。

## 2. 入出力とプログラム構成

### 2.1 再利用する既存コード

| 対象 | 観察した現状 | 採用する方針 |
|---|---|---|
| `src/AnalysisVSPAERO.py::vsp_sweep` | API解析のGeometry計算とSweep、`alpha[0]`・`alpha[-1]`・点数による等間隔指定、Beta引数なし | 単点解析を新たな小さなラッパで統制。既存関数は互換性を壊さず再利用または局所拡張 |
| `src/AnalysisVSPAERO.py::vsp_stability_derivatives` | Steady/Steady Adjoint、診断辞書、実効設定を返す | 安定微係数専用の実行アダプタとして再利用。成功フラグだけでV&V合格にしない |
| `src/AnalysisVSPAERO.py::vsp_trimmed_sweep` | 固有Geometry名・舵面グループ名への依存、参照面積モードの変更経路がある | 汎用Runの入口には直接用いない。トリム試験で明示的設定・読み戻し確認を条件として利用 |
| `src/VSPAEROMesh.py` | OpenVSP形状の分割寸法・細分化調整 | 前処理の分割状態の測定に利用。VSPAEROの交線処理後の最終メッシュ品質と同一視しない |
| `src/VSPAEROStab.py` | `.stab`を解析し派生微係数を扱う | パーサを再利用。`ConGrp_n`の舵面名への推定フォールバックはV&V認定用途で禁止 |
| `src/util.py` | パラメータ検索・入力設定・Results読出し | APIアダプタに再利用。要求値→読み戻し値の照合を追加 |
| `tests/` | 既存機能の回帰テスト | 変更に対する回帰テストに使用。CFDのV&Vとは明確に区別 |

### 2.2 新規機能の責務（小規模な4モジュール）

| 新規モジュール | 主な責務 |
|---|---|
| `VVRunner.py` | Test Caseの具体化、事前検査、Runの隔離起動、実効設定読戻し、成果物収集、状態イベント記録 |
| `VVDatabase.py` | SQLiteのスキーマ初期化、トランザクション、データの型・参照・ハッシュ整合性、書込み・読出し |
| `VVPostprocess.py` | `.polar`/`.stab`等の取り込み、QoI定義、数値収束、中央差分、ドラッグ補正、参照照合、比較と根拠生成 |
| `VVReport.py` | ケース一覧・Excel/CSV・図・状態/判定一覧の出力。正本は変更しない |

`schema.sql`はDB構造の正本、`test_catalog.csv`は上位計画書の25件のTest IDカタログ、`templates/study.json`は新規Study設定のひな形とする。

### 2.3 想定コマンド

```text
python -m vv init --study study.json
python -m vv plan --study study.json
python -m vv check --case GRID-01-MEDIUM
python -m vv run --case GRID-01-MEDIUM
python -m vv import-reference --file references_wt.csv --metadata wt.json
python -m vv evaluate --test GRID-01
python -m vv export --format xlsx
python -m vv verify-artifacts --all
```

実装時にCLIの正確なオプション名を固定する。ここに記したコマンドは**将来実装するインターフェース仕様**であり、現行添付コードで実行できると主張していない。

## 3. SQLiteスキーマ（DDLは同梱）

### 3.1 エンティティと責務

| テーブル | 主キー | 必須の考え方 |
|---|---|---|
| `study` | `study_id` | 対象バージョン3.50.5、計画版・計画ハッシュ・設定ハッシュ |
| `test_definition` | `test_id` | 計画書に存在する**25件**、必須／条件付き、実行可能区分 |
| `geometry` | `geometry_id` | G1/G2/G3/Thick化などの**物理形状・解析表現**。Thin/Thick SetをJSONで保持 |
| `flight_condition` | `condition_id` | Mach、Re、α、β、舵角、トリム、Sref/bref/cref、基準点、軸規約 |
| `numerical_setting` | `setting_id` | Tessellation、Wake、ソルバーモード、要求パラメータ一式とハッシュ |
| `planned_case` | `case_id` | Test ID＋ケース役割（Coarse/Medium/Fine、W0など）＋予定入力 |
| `run` | `run_id` | **1回の実行**、入力参照・状態・要求/実効/生成条件・実行バイナリハッシュ |
| `case_point` | `point_id` | `.polar`等の条件点行。Run内の順序、入力との条件一致結果 |
| `case_run_link` | `(case_id,point_id)` | 同一解析点を複数のテストへ使う多対多対応 |
| `dataset` | `dataset_id` | VSPAERO／パネル法／NS／風洞／理論／外部CD0／派生値、Calibration/Validation区分 |
| `artifact` | `artifact_id` | 生ファイル、設定、ログ、図の相対パス・SHA-256・サイズ・帰属 |
| `quantity` | `quantity_id` | 定常係数・微係数・CD0・分布要約値。生値と正規化後の値を分離 |
| `quantity_dependency` | `(derived_quantity_id,input_quantity_id,role)` | 方式A/B、中央差分、傾斜計算などの**値単位の入力依存性** |
| `acceptance_rule` | `rule_id` | QoI、数値／工学／適用の閾値・出典・暫定/承認状態 |
| `comparison` | `comparison_id` | 比較する2つの量、差、不確かさ、条件整合、評価結果 |
| `assessment` | `assessment_id` | 試験・QoI・適用領域ごとの数値/工学判定、使用認定、根拠、承認 |
| `run_event` | `event_id` | 状態変更・エラー・理由を時系列の追記専用イベントとして保持 |

### 3.2 テーブル間の主な関係

- `test_definition (1) - (N) planned_case`。
- `geometry`, `flight_condition`, `numerical_setting` はそれぞれ複数の `run` から参照可能。
- `run (1) - (N) case_point`。初期の新規実行はN=1とする。
- `planned_case (N) - (N) case_point` は `case_run_link`で接続。Runを作り直しても旧リンクを削除しない。
- `case_point (1) - (N) dataset (1) - (N) quantity`。例えば一つの解析点に `.polar` と `.stab` に由来するDatasetを別々に保持可能。
- 参照DatasetはRunを持たず、元資料・取得条件・証拠Classを持つ。理論値と風洞値を混ぜない。
- 補正後の値は`kind=DERIVED`のDatasetに格納し、`quantity_dependency`で入力数値を参照。
- `comparison`は2つの`quantity`を参照し、`assessment`はその判定・適用範囲の記録を担う。

### 3.3 生データの不変性と数値定義

- `.polar`の`CLi`, `CDi`, `CLiw`, `CDiw`, `CDo`, `CLtot`, `CDtot`等を**元の名前・元の値**で保存する。`CLiw`はwake側評価であり、`CLi`に加算しない。
- 係数の正規化と座標変換を施す場合、`native_name/value/unit`と`canonical_value/unit`, `convention_id`, `normalization_method_id`を区別する。
- 微係数には、Steady FD／Steady Adjoint／P/Q/R、rad当たりの定義、
  `p_hat=pb/(2V)`, `q_hat=qc/(2V)`, `r_hat=rb/(2V)`、定常・非定常組合せの違いを記録する。動安定の`q`と`q+alpha_dot`を同一QoIに正規化しない。
- `.stab`の内蔵寄生抗力由来成分の有無を確認できない場合は、`variant=INTERNAL_TOTAL_UNRESOLVED`等の表現を用い、非粘性単独微係数として扱わない。
- 欠測は`NOT_CALCULATED / NOT_AVAILABLE / NOT_APPLICABLE / INVALID`とし、ゼロを埋めない。原値があっても正規化定義が未確認なら`AVAILABLE`に昇格させない。
- `quantity`の数値がNaN・InfならSQLへ登録しない。元ファイルは保存し、該当値を`INVALID`にして理由をログへ残す。

### 3.4 データ整合性の検査

DBへ接続したら**毎回**`PRAGMA foreign_keys=ON`を実行する。DBはトランザクションで更新し、関連行の部分登録を避ける。SQLiteの制約だけでは足りないため、アプリケーション側で次を追加検査する。

1. `case_run_link`により結び付いた予定ケースとRun点が同じStudy、形状、飛行条件、解析設定を指していること（試験で許可した差分を除く）。
2. 同じ`run_id`の入力ファイルと出力ファイルの相対パスがそのRunディレクトリ配下にあり、パストラバーサルしないこと。
3. ハッシュとファイル内容の一致（SHA-256）。単にDB上にハッシュが記録されているだけでは「検証済み」としない。
4. `dataset.kind`と`point_id`・`usage_class`の整合性。計画調整に用いた風洞点を独立Validationとして流用しない。
5. `comparison`の両値のQoI、単位、軸、基準量、Mach/Re/α/β/舵角、補正条件が一致していること。
6. `quantity_dependency`の派生値の入力が存在し、循環参照がなく、計算式の版数が明記されること。
7. `assessment`で「適用可」を承認する場合、承認済み`acceptance_rule`、数値受容、実験Validation、根拠、レビュー者がそろうこと。

DBファイルの変更履歴はスキーマ版で管理し、変更時は事前バックアップとマイグレーションを実施する。通常の直接SQL編集は認めない。

## 4. 設定ファイル仕様（`study.json`）

### 4.1 必須ブロック

| ブロック | 内容 | バリデーション |
|---|---|---|
| `schema_version`,`study_id` | 設定のバージョン・Study識別子 | 既存Studyの版整合が必須 |
| `plan` | Markdown計画書の版・SHA-256 | 現物と照合、不一致なら実行拒否または明示的改訂承認 |
| `runtime` | 3.50.5のバイナリ、Python環境、プロセス数 | バージョン一致、実行ファイル実在、SHA-256記録 |
| `geometry` | G1/G2/G3等、Thin/Thickの構成、Set名 | Setの重複不可、名前とIDの対応を実行時読み直す |
| `reference` | 面積/長さ/重心・軸・単位 | **解析前に数値確定**。null不可 |
| `flight_condition` | Mach、Re、α、β、舵角、トリム状態 | 値が具体的なケース化までに確定すること |
| `numerics` | 面分割・Wake・解析モード等 | 3.50.5のAPI/Parmが実際に受け付けた値まで照合 |
| `policies` | 自動補間、追加出力行、欠測等の扱い | デフォルトは誤比較を防ぐ厳格側 |
| `screening_only` | 暫定的なMedium/Fine差の目安 | 工学的要求精度として承認済み扱いしない |
| `drag` | CD0の外部データID、A/Bの計算式・基準条件 | 未指定ならDRAG評価を保留、他の解析は可能 |

サンプルファイルでは機体固有の`mach`, `re_cref`, `sref`, `bref`, `cref`, `moment_xyz`, `axes_convention_id`等を`null`にしている。`init`と`plan`はそれを許容し、`check`と`run`は対象ケースの必要情報の不足を明確なエラーとして通知する。

### 4.2 パラメータの標準的な保持名とOpenVSPとの対応

| 論理項目 | OpenVSP/VSPAERO側 | 対応方針 |
|---|---|---|
| 主翼/尾翼弦方向 | `Tess_W` | Geometry Parmとして更新。設定した数字と最終パネル数を区別 |
| 各翼幅／胴体区間 | `SectTess_U` | XSec・セクション単位で設定し読み戻す |
| 前縁集中度 | `LECluster` | 0.25を基準、0.5/1.0を追加感度 |
| 胴体周方向 | `Tess_W` | 17/25/33など、モデルごとの制約を確認 |
| Wake反復 | `WakeNumIter` | 初期5、試験5/10/20、要求・実効・`.vspaero`実際の`WakeIters`を対応付け |
| Wake節点 | `RootWakeNodes`／解析入力`NumWakeNodes` | 命名差を1つの論理変数へ対応付け、生成された値も記録 |
| Set選択 | `GeomSet`, `ThinGeomSet` | Geometry生成とSweep/安定解析の双方で一致することを照合 |
| 参照係数 | `RefFlag`, `Sref`, `bref`, `cref`, `Xcg`等 | `MANUAL_REF`を基本。既存トリム関数が変更し得るため再検査 |
| 安定解析モード | `STABILITY_OFF`, `STABILITY_DEFAULT`, `STABILITY_ADJOINT`、P/Q/R | 実行環境の定数を取得し、保存は人可読のモード名と実効値 |
| 操舵 | `ActiveFlag`, `DeflectionAngle`, `Surf_<SubSurfID>_<Index>_Gain` | Control Groupと実際の左右偏向方向を対応付け、推定フォールバックを禁止 |

パラメータが現行API/Parmに見つからない場合は、**黙って無視せずPrecheck Failedにする**。事前設定で存在が確認できない補助設定は対象ケースから除外または実装保留とし、別名を推測して適用しない。

### 4.3 設定の優先関係

`study.json`の基準設定に、Test IDごとの明示的な変更を重ね、その結果を`requested_snapshot.json`として凍結する。OpenVSPのReadVSPFile後に**値を読み直して**`effective_snapshot.json`を作成する。さらに`.vspaero`・`.vspgeom`から分かる範囲で`generated_snapshot.json`を作る。

比較時の階層は以下を区別する。

- **Requested**：ユーザー／テストカタログが要求した条件。
- **Effective**：OpenVSP API・Parmから実際に読み取れた条件。
- **Generated**：VSPAEROへ実際に渡された設定・生成メッシュの情報。

要求値と実効値が一致しない場合は差分を残す。全てを読み出せないパラメータが存在する場合、当該項目は`UNVERIFIED`とし、設定不変性が要件のテストでは実行・認定を保留する。

## 5. Run IDと状態管理

### 5.1 ID／保存フォルダ

- Study ID：例`VSP3505-VV-001`（固定）。
- Test ID：計画書の`GRID-01`等、変更不可。
- Case ID：例`GRID-01-MEDIUM`、`WAKE-01-W2`。同一Test内で役割は一意。
- Run ID：`RUN-000001`からの単調増加番号（採番はDBトランザクションで直列化）。再試行は新ID。
- Point ID：`RUN-000001-P000`など。初期Runは1点だけ。
- Dataset ID、Quantity ID等も一意IDとし、行番号やファイル名から機械的に同一性を推測しない。

### 5.2 推奨状態遷移

```text
PLANNED
  | precheck OK
  v
PRECHECKED
  | isolated process starts
  v
RUNNING
  | solver finished and native files collected
  v
CAPTURED
  | parsing, point-to-condition and file-integrity checks OK
  v
PARSED
  | all mandatory result families present and accepted
  v
SUCCEEDED
```

例外遷移：

- `PLANNED → PRECHECK_FAILED`（設定・形状・ソフトウェア不一致）。
- `PRECHECKED/RUNNING → EXECUTION_FAILED`（プロセス異常、ソルバー実行失敗）。
- `RUNNING/CAPTURED/PARSED → OUTPUT_INVALID`（必須ファイル欠落、解析点不一致、パーサ失敗、非有限値）。
- `PRECHECKED/RUNNING → INTERRUPTED`（タイムアウト、クラッシュ、強制中断などで実行完了を確認できない）。
- 非完了状態から明示的な操作で`CANCELLED`。
- **完了・失敗状態のRunを再びRUNNINGへ戻さない**。再実行は新しいRun ID＋`retry_of_run_id`で記録。

すべての状態変更を`run_event`へ追記する。状態変更とイベント追記は一つのSQLiteトランザクションで実施する。`SUCCEEDED`の意味は「結果取得・条件一致・解析出力の検査を通過」であり、Numerical AcceptanceやValidationのPASSとは無関係。

### 5.3 事前検査の判定順序

1. ファイル実在とVSPAERO 3.50.5版・配布バイナリのハッシュ確認。
2. 必要な形状・条件・参照値がすべて具体値かを確認。
3. Thin/Thick Setが重複せず、期待する胴体・主翼・尾翼が選択されていること。
4. 舵面テストの場合、Control Group、Gain、ActiveFlag、実偏向方向が期待どおりであること。
5. `.vsp3`をRun専用作業フォルダへ複製し、当該コピーのみに変更を適用。
6. OpenVSPのAnalysis API設定とParm設定を読み戻し、単位・有効範囲を照合。
7. 書き出しモデル・Geometry/Wakeの生成条件を記録し、意図しない変更がないこと。
8. `FixedWakeFlag`有効時にWake反復試験が指定されていないことなど、試験別の論理矛盾を検出。
9. 全件合格で`PRECHECKED`。実行前に物理的な妥当性を完全自動判定するものではないため、メッシュ交線品質などは目視レビューも保持。

### 5.4 成果物の確定

Runフォルダは`03_run_raw/RUN-xxxxx/{input,output,logs}`とし、出力収集後にSHA-256を計算し`manifest.json`とDBへ記録する。事後の検査は`verify-artifacts`で行う。元ファイルの差替えを検出した場合、関連するComparison/Assessmentを再評価待ちにする。ファイル変更を自動で「同じRunの新結果」として受け入れない。

## 6. Test ID別の実行手順と判定ロジック

**全25件の詳細カタログは`test_catalog.csv`。** 本節は判定フローとして位置付け、各Test IDの必要条件を省略せずに定義する。

### 6.1 VER／GEO

| Test ID | 実行 | 判定（必要な情報が足りなければINCONCLUSIVE） |
|---|---|---|
| VER-01 | Hershey Barの理論参照と対応する3.50.5ケース | `CLi`・誘導抗力・分布の誤差と分割感度を記録。理論解の適用条件・承認された閾値が必要 |
| VER-02 | Ellipsoid、Thick、理論`Cp` | 位置対応した`Cp`のRMS/最大差。面圧規約と座標対応が不明なら合否不可 |
| GEO-01 | G1：主翼単体 | メッシュの健全性と基礎係数を確認。単独で実験Validationとはしない |
| GEO-02 | G2：主翼＋胴体 | G1と同じ基準面積・基準点で`ΔCL`,`ΔCm`等を計算。モデル差の符号と整合を確認 |
| GEO-03 | G3：全機 | G2からの尾翼追加差分、6分力と翼尾干渉。物理形状差を数値誤差へ混入させない |
| GEO-04 | Thin/Thick変更（条件付き） | 新Geometry IDとして形態差を測定。数値分割収束とは別レポート |

### 6.2 GRID／WAKE

| Test ID | 予定Run・条件 | 数値評価 |
|---|---|---|
| GRID-01 | 全機G3同一条件、Coarse/Medium/Fine | `ΔCM=QM−QC`、`ΔMF=QF−QM`、メッシュ要素数と品質を算出。非単調・格子非相似ならGCIは保留 |
| GRID-02 | 干渉が強い第2飛行条件、同じ3水準 | GRID-01と同じ処理。代表1条件のみで全条件へ一般化しない |
| GRID-03 | `LECluster`、局所分割、`NRef`等の1因子変更（条件付き） | メッシュ品質と局所応答。形状トポロジーの変化は独立に注記 |
| WAKE-01 | Wake反復5、10、20、その他同一 | `Δ5,10`,`Δ10,20`、残差とWake品質、各係数の安定性。反復が増えるだけでPASSとはしない |
| WAKE-02 | Wake節点8、12、16。Wake回数はWAKE-01で選定 | 節点数のみ変更した係数差。`CLi/CLiw`差に一律合格閾値を設けない |
| WAKE-03 | Explicit/Implicit/Fixedを代表条件で比較（条件付き） | 方式の差として保存。方式間一致は必須ではなく、各方式の適用性評価とする |

GRIDでは`CLi`,`CMiy`,`CDiw`のMedium/Fine絶対差の**暫定スクリーニング値**（0.005/0.002/0.0005）を表示する。ただし閾値が暫定の間、正式なNumerical PASSを自動付与しない。`CLiw`、`CDi`も保存し、同一系列の収束を別々に評価する。

GCIの自動算出に必要な条件：格子幅の定義・細分化比・計算量の対応、漸近収束の妥当性、QoI差の数値的有意性が確かめられたときに限る。判断できなければ`GCI_NOT_APPLICABLE`を記録して単なる差分表を出す。

### 6.3 XCD／VAL

| Test ID | データ準備・比較 | 判定 |
|---|---|---|
| XCD-01 | パネル法結果を共通CSVで取込み、形状/後流/軸/基準量を対応付け | 差はDiagnostic。両コード一致を実験ValidationのPASSとして扱わない |
| XCD-02 | NS結果を取込み、乱流/遷移/格子収束情報を登録 | 数値誤差と物理モデル差を分ける。NS計算結果を真値と仮定しない |
| VAL-01 | 風洞の定常縦特性、迎角3～5点等 | `CL,Cm,CD`の差・要求精度・実験不確かさに基づく。未対応の誤差補正は明記 |
| VAL-02 | 風洞の横滑り特性 | `CY,Cl,Cn`および傾斜を比較。横滑りデータがなければ証拠不足 |

参照の利用区分は`CALIBRATION`（モデル・補正の調整）、`VALIDATION`（設定固定後の独立検証）、`REFERENCE`（その他）を保持する。`Class C`のデータだけで工学的なPASSを確定しない。

参照との数値比較は同一QoI・同一定義・同条件を原則とし、単純に最近傍の迎角データを選択しない。補間する場合は補間元点、アルゴリズムの版数、補間範囲、不確かさを別途記録し、`comparison_kind=INTERPOLATED`とする。

### 6.4 STA／CTRL

| Test ID | 実行 | 判定 |
|---|---|---|
| STA-01 | 定常Steady AdjointとSteady有限差分を同条件で実行 | `Cm_q`,`Cl_p`,`Cn_r`等の符号・単位・揚力/抗力成分・収束を照合。相互一致のみでValidation PASSにはしない |
| STA-02 | Medium/Fine、Wake設定、必要なら有限差分幅変更 | 微係数の絶対差、符号、飛行力学的感度を評価。設計閾値未承認なら数値判定保留 |
| STA-03 | 条件付きP/Q/R。非定常NS/風洞があれば比較 | 角振動と`q+αdot`等の組合せ微係数を別QoIに登録。動的参照なしなら動的Validationは証拠不足 |
| CTRL-01 | エレベータ、0/±2/±5/±10°（舵角制限内） | `Cm,CL`の有限角変化と±2/±5中央差分（rad^-1）を独立保存 |
| CTRL-02 | エルロン同様 | `Cl,Cn,CY`、逆ヨー、左右Gainの符号、有限角の非線形性 |
| CTRL-03 | ラダー同様 | `Cn,CY,Cl`、偏向正負、有限角/微小舵角の比較 |

舵面の`ConGrp_n`を`AILERON/ELEVATOR/RUDDER`へ自動で割り当てることはしない。OpenVSPのControl Groupの名前、所属SubSurface ID、Gain、実際の偏向方向が確認された場合のみ正規化した舵効きとして登録する。

中央差分：

\[
Q_{\delta}\approx \frac{Q(+\Delta\delta)-Q(-\Delta\delta)}{2\Delta\delta},\qquad \Delta\delta\ \text{はrad}
\]

差分元のRun/Point/Quantity ID、±の舵角、基準飛行条件を`quantity_dependency`へ記録する。角速度微係数も無次元化定義を確認しない限り単位変換しない。

### 6.5 DRAG／ACC

| Test ID | 実行 | 判定 |
|---|---|---|
| DRAG-01 | 外部Parasite Drag Tool等から`CD0_build`を取り込み、`CDiw + CD0_build`を計算 | Mach/Re/基準翼面積/構成要素の積み上げ一致を確認。出力は別のDERIVED Dataset |
| DRAG-02 | `CDiw+(CDo-CDo_baseline)+CD0_build`、Bは比較候補 | 同じMach/Reで定義した`CDo_baseline`を明示。未定なら未計算、Aとの比較を風洞で検証 |
| ACC-01 | 主要Test/QoI/飛行領域を集約しレビュー | Numerical・Engineeringの独立判定と制約を承認。自動認定は禁止 |

`CD0`の**寄生抗力成分**と`CDmin`の**最小全抗力**を同一視しない。`CDtot`,`CDwtot`へ外部CD0を直接足さない。外部CD0は揚力・方向微係数への局所寄与を表せない場合があるため、`.stab`中の内蔵粘性寄与を外部CD0体系へ無条件に流用しない。

## 7. 比較量、誤差、受容・認定の判定

### 7.1 比較誤差

比較可能な同一QoIに対して

\[
E_Q = Q_{\mathrm{VSP}}-Q_{\mathrm{ref}},\qquad
|E_Q|\ \text{と}\ E_Q\ \text{の両方を記録}
\]

相対誤差は基準値がゼロ近傍のとき未定義とし、`NULL`を保存する。適当な微小定数を勝手に分母へ加えない。条件・符号規約が未一致なら`DIRECT`比較不可とする。

### 7.2 不確かさ

独立で互換的な**標準不確かさ**がある場合のみ、例えば

\[
u_{\rm val}^2=u_{\rm num}^2+u_{\rm input}^2+u_{\rm exp}^2
\]

を計算する。相関がある場合は共分散を含める。GCIの値を標準不確かさへ無条件に読み替えない。`comparison.uncertainty_metadata_json`に、各項の意味、標準／拡張、包含係数、信頼水準、共分散の考え方を記録する。情報がなければ不確かさ列はNULLとする。

### 7.3 Numerical Acceptance

1. Run/Result状態がSUCCEEDED、形状・Set・メッシュ・解析点が正しい。
2. 主要係数について少なくとも指定水準の分割・Wake感度試験を実施。
3. 評価量の数値変動または適切に推定した数値不確かさが、**承認済み**の数値受容基準以内。
4. 不規則収束・非有限値・重大なメッシュ不良がない。
5. 閾値未承認、収束判断不能、特定QoI欠測なら`INCONCLUSIVE`。暫定数値は「スクリーニング結果」でありPASSの根拠にしない。

### 7.4 Engineering Acceptance（例示的な保守的ルール）

`T_Q`を用途別に承認された許容誤差、`U_Q`を同一定義で扱える比較不確かさの**幅**とする。内部の保守的な判定ルールとして、次を採用することを提案する。

- \(|E_Q|+U_Q\le T_Q\) なら、その他の品質要件が満たされる場合に限り工学的PASSの候補。
- \(|E_Q|-U_Q>T_Q\) ならFAILの候補。
- その中間は`INCONCLUSIVE`として技術レビューへ回す。

これはASME V&V 20がすべての用途に要求する公式の合否判定式ではなく、本プロジェクトの**安全側の意思決定案**である。`T_Q`、`U_Q`、証拠Class、対象飛行範囲、モデル変更履歴が承認されていない場合は自動PASSにしない。Class Bのデータしかない場合は根拠の制約を明示し、`CONDITIONAL`か`INSUFFICIENT_EVIDENCE`を技術判断する。Class Cでは正式な「適用可」を認めない。

### 7.5 総合判断ロジック

```text
if run_output_invalid:
    numerical = INCONCLUSIVE or FAIL (原因による)
    engineering = NOT_EVALUATED
    disposition = INSUFFICIENT_EVIDENCE
elif numerical_approved_rule_missing or numerical_evidence_missing:
    disposition = INSUFFICIENT_EVIDENCE
elif numerical == FAIL:
    disposition = NOT_APPLICABLE
elif engineering_reference_missing or engineering_approved_rule_missing:
    disposition = INSUFFICIENT_EVIDENCE
elif numerical == PASS and engineering == PASS:
    disposition = APPLICABLE  # さらに人手で対象領域・条件を承認
elif any_required_result_is_FAIL:
    disposition = NOT_APPLICABLE
else:
    disposition = INSUFFICIENT_EVIDENCE

# 限定用途なら技術レビューでCONDITIONALとし、制限/補正/承認者を必須記録
```

検証対象が定常データのみで、動的参照がない場合は、動安定微係数の実験ValidationをPASSにしない。`ACC-01`は**評価量×形状×飛行条件×用途**単位で認定範囲・除外範囲を出す。

## 8. 外部参照データの共通インポート

初期は共通CSVとメタデータJSONを採用する。ソルバーごとの独自ファイル自動読込は将来機能とする。

```text
source_type,source_case_id,mach,re_cref,alpha_deg,beta_deg,controls_deg_json,
qoi,value,unit,axes_convention_id,sref,bref,cref,moment_x,moment_y,moment_z,
uncertainty_value,uncertainty_type,coverage_factor,usage_class,evidence_class
```

単位系、基準量、モーメント基準点、風洞支持干渉・壁面補正、NS乱流モデル・収束、参照ファイル出典をメタデータJSONへ含める。Reference Datasetを登録するときに元CSV/JSONのハッシュを取る。入力不足であれば記録は可能でも比較可否は`ILLUSTRATIVE`、使用認定は証拠不足とする。

## 9. 試験結果のエクスポート

- `TestMatrix`：25件すべてのTest ID、必須/条件付き、必要Run、欠測、実行状況。
- `RunRegister`：Run ID、モデル、飛行条件、数値設定、ソフトウェア・ハッシュ・ログ・状態。
- `Coefficients`：生の`CLi/CDi/CLiw/CDiw/CDo`等と補正後の別値。
- `Derivatives`：静/動微係数・舵効き、無次元化/符号/角速度・舵角幅。
- `GridWakeStudy`：格子数・節点数・Wake反復・収束図の参照先。
- `ReferenceComparison`：比較可否、誤差、不確かさ、補間条件。
- `Assessment`：QoI×飛行条件×用途のNumerical/Engineering判定、制約と承認者。
- `ArtifactManifest`：原本/生成物ファイルのパスとSHA-256。

Excelは**閲覧・レビュー用途のビュー**であり、正本の代用にしない。CSVはUTF-8（Excelでの扱いが必要な場合はUTF-8 with BOM）とする。大量の`Cp`や時系列は別ファイルに置き、DBには相対パスとハッシュだけを登録する。

## 10. ソフトウェア受入試験（最低限）

| ID | 目的 | 合格条件 |
|---|---|---|
| SW-01 | 同一条件の再実行 | 新Runを発行し旧結果は変更されない |
| SW-02 | Runの状態遷移 | 禁止遷移を拒否、イベント欠落なし、失敗原因保存 |
| SW-03 | 1 Run 1点条件確認 | 2点以上の意図しない出力でOUTPUT_INVALID |
| SW-04 | Thin/Thickと参照値の事前検査 | Set重複・欠損・無効な舵面等を検出 |
| SW-05 | `.polar`読込 | サンプル値・出力順・生列名・有限値が一致 |
| SW-06 | `.stab`読込 | 実名と`ConGrp_n`の対応不明では認定値を出さない |
| SW-07 | 単位・符号 | 軸変換、deg/rad、角速度無次元化の既知例が一致 |
| SW-08 | GRID/WAKE | 許可された設定以外の変更を検出、差分は正しく計算 |
| SW-09 | CD0 A/B | 既知の手計算値と一致、CDoの二重計上なし |
| SW-10 | CSV参照比較 | 条件不一致時にDIRECT比較を拒否、未承認補間を行わない |
| SW-11 | 不確かさ／判定 | 未定閾値・欠測・Class Cは自動PASSにならない |
| SW-12 | ファイル完全性 | 改変/欠落を検出、関係する判定は要再評価 |
| SW-13 | Excel/CSV出力 | DBの値、行数、試験・Runリンクと一致 |
| SW-14 | 25試験登録 | 計画書とカタログのID集合が一致し、欠落・重複なし |

SW-01～14のうちAPI不要の検査は単体テストで実施。WindowsのOpenVSP 3.50.5実機でGeometry生成、Wake設定、解析、出力、`.stab`取得、舵面偏向を含む結合試験を別途実施する。

## 11. 実装の段階と完了ゲート

| 段階 | 実装対象 | 完了ゲート |
|---|---|---|
| P1 | DDL・Study/25テスト登録・Run状態管理・ファイル保全・`.polar`/`.stab`パーサ | サンプル取込/書込み、SW-01～07・SW-12/14が通る |
| P2 | Geometry/GRID/WAKEのケース生成、単点実行、数値感度レポート | Windows 3.50.5で同じ条件の実行・比較を再現、SW-08が通る |
| P3 | パネル/NS/風洞CSVと比較、数値/工学判定 | 参照条件の不一致検出、証拠不足判定、SW-10/11が通る |
| P4 | Steady/Adjoint・CTRL・外部CD0方式A/B・Excel | 動安定/舵効き/抗力補正が原値と分離され、SW-09/13が通る |
| P5（条件付き） | P/Q/R解析自動化、適格GCI、将来の並列化 | 対応する専用テスト・実験データで受入できた場合のみ追加 |

初期段階で未承認の飛行条件・誤差閾値があっても、Study初期化、試験計画、パーサ、Database機能は開発できる。ただし不足する具体値をダミーで埋めて実行成功・Validation PASSを演出しない。

## 12. 未確定事項と仕様として決めておくこと

| ID | 未確定事項 | プログラム上の取り扱い |
|---|---|---|
| O-01 | 実機Mach/Re/α/β/舵角範囲 | `study.json`のnull／未指定。ケースを具体化するときに必須 |
| O-02 | `Sref`,`bref`,`cref`、基準点、軸定義 | 実行Precheckの必須条件。値は推測しない |
| O-03 | 数値/工学許容誤差 | `acceptance_rule`の`PROVISIONAL`扱い。承認前は最終PASS不可 |
| O-04 | 風洞の不確かさと使用区分 | データClassとCalibration/Validation区分を明示。不十分なら証拠不足 |
| O-05 | P/Q/Rの非定常データと必要度 | `STA-03`はCONDITIONAL。対応データがなければ未認定 |
| O-06 | 3.50.5の実行バイナリとAPI可用性 | 実機で解析入力名・値の型・生成ファイルを照合して実装を最終確定 |
| O-07 | Runごとの実際のメッシュ品質判定方法 | 初期はログ・要素数・画像によるレビュー。定量メッシュ診断は将来拡張 |

本仕様が規定したDBの内部構造は、機体固有値を確定しなくても実装可能である。一方、**飛行条件、比較条件、承認済み要求精度なしに正式なV&V合否を出さない**という運用上のゲートは維持する。

## 13. 参考・トレーサビリティ

- 上位計画書：`plan.md` 第3～12章および付録。
- 利用者提供のコード：`OpenVSP-main_20260927_trimmed_sweep_activeflag_v1.zip`、特に`AnalysisVSPAERO.py`、`VSPAEROMesh.py`、`VSPAEROStab.py`、`util.py`。
- OpenVSP 3.50.5ソース：https://github.com/OpenVSP/OpenVSP/tree/OpenVSP_3.50.5 （ソースと配布バイナリの対応は別途実機確認）。
- 公式API：https://openvsp.org/api_docs/latest/group___analysis.html 、https://openvsp.org/api_docs/latest/group___results.html （「latest」は3.50.5固定ではないため入力の最終確認は対象実行環境で行う）。
- V&V方法論：ASME V&V 20、AIAA G-077、Celik et al. (2008)、Eça & Hoekstra (2014)。詳細は上位計画書の参考文献による。

## 改訂履歴

| 版 | 日付 | 内容 |
|---|---|---|
| 1.0 | 2026-10-21 | 上位計画書に基づくDB構造、設定ファイル、Run管理、25試験、判定、実装段階、受入試験の詳細設計初版 |
| 1.0.1 | 2026-10-21 | 上位計画書のv1.1化に伴う参照ファイル名の更新。schema.sql・test_catalog.csv・判定処理・DB設計は変更なし |
