# VSPAERO 3.50.5 V&V評価報告書 第1～4章 自動入力機能 詳細実装仕様書

**仕様書番号**：VSP-VV-SPEC-REPORT-AUTOFILL-001  
**版**：1.0（実装前仕様・レビュー草稿）  
**作成日**：2026-11-04  
**対象リポジトリ**：`OpenVSP-main_v20261104.0`  
**対象実生成物**：`offline_regenerated_20261104T035320Z_8a052a77.zip`  
**対象機能**：`src/openvsppy/vv/VVReport.py` と報告書テンプレート、必要な読み取り専用データ抽出処理  
**目的**：第1～4章および表紙の「未設定」「未確認」を、証拠により確定できる事実に置き換える。認定・承認を自動で作成しない。

> **位置づけ**：本書は修正実装に着手するための詳細仕様であり、現在のコードに実装済みであること、VSPAEROのV&Vが完了したこと、または解析結果が工学的に承認されたことを示す文書ではない。

---

## 1. 監査で確認した現行状態

### 1.1 監査対象と確認方法

- `offline_regenerated_20261104T035320Z_8a052a77.zip` 内の `evaluation_report.md`、`database_snapshot.sqlite`、`results_snapshot.json`、`evidence_manifest.json` を照合した。
- `OpenVSP-main_v20261104.0.zip` の `VVReport.py`、`VVPostprocess.py`、`VVAutomation.py`、公式レポートテンプレートおよび `db_report_field_mapping.md` を確認した。
- SQLiteは読み取り専用で照会し、過去のRun、解析モデル、参照資料、評価状態は変更していない。
- 表紙の未入力候補は7か所、第1章は17、第2章は28、第3章は19、第4章は19か所（`《...》`出現数）。**この数は空欄の個数であって、自動確定可能項目数ではない**。加えて第4章には「未確認」が9か所ある。

### 1.2 現在の固定DBで確認された登録件数

| テーブル・情報 | 件数 | 報告書での使用上の注意 |
|---|---:|---|
| `study` | 1 | 計画・DB識別子。コードの実行時検証済み版とは区別 |
| `test_definition` | 25 | 計画したV&V試験の一覧 |
| `planned_case` | 69 | 成功した計算Case数とは異なる |
| `geometry` | 24 | G103AとVER系ベンチマークが混在 |
| `flight_condition` | 25 | 異なるCase IDで同一数値条件も存在 |
| `numerical_setting` | 59 | Requested設定を表す。実行実績ではない |
| `run` | 60 | SUCCEEDED 43、EXECUTION_FAILED 9、INTERRUPTED 8 |
| `dataset` | 46 | VSPAERO 43、DERIVED 2、EXTERNAL_CD0 1 |
| `quantity` | 1,984 | 採用Runに紐づく報告書側のQuantityは1,937。母集団を混同しない |
| `artifact` | 1,935 | 記録内容を原本の完全性／解析の妥当性と混同しない |
| `acceptance_rule`／`comparison`／`assessment` | 各0 | 工学的PASS／使用認定を作り出してはならない |

採用成功Caseは42件。単一Runの複数Case紐付けや過去Runの重複があり得るため、**成功Run総数43と採用成功Case数42を等価に扱わない**。

### 1.3 認められた既存実装の特徴

- `write_evaluation_report()` は公式テンプレートの `<!-- VV:META:... -->`／`<!-- VV:CHAPTER:... -->` マーカーを置換する。
- 第1章の冒頭、5章の実施状況、第6～13章の一部は自動生成されるが、1.1～4.4のテンプレート本体の大半が固定の`《未設定》`のままである。
- `collect_study_results()` はCase単位で採用成功Runを選び、`results_snapshot.json`に固定する。これを第1～4章にも利用する。
- `run.requested_snapshot_json`、`effective_snapshot_json`、`generated_snapshot_json` にそれぞれ要求値、API取得実効値、生成ソルバー設定が保存されている。
- 旧Studyの計画書ハッシュと現行コードの計画書は一致しないため、**過去のStudyを再ロードして計画書の版を強制更新してはならない**。

### 1.4 追加で確認した特有の不整合

1. `GEO-01-G1-G` と `GEO-02-G2-G` は `representation_json.case_role` がそれぞれ `G1`／`G2`である一方、`geometry.stage` はどちらも `G3`である。**Geometry IDからStageを無条件に確定してはならない**。
2. `GEO-04-THICK-G` の `physical_model_json.model.wing` は `Thin` のまま。一方、成功Case `GEO-04-THICK` の `numerical_setting.requested_json.geometry_membership` は主翼をThick側に置き、`effective_snapshot_json.geometry_sets` でも主翼に対応するGeometry IDがThick側で記録されている。**宣言メタデータと実行設定の差異**として明示する。物理的な同等性承認は自動で行わない。
3. 解析条件の`Mach=0.001`、`0.1`等はVERベンチマーク由来。G103A全機に対して計算・妥当性確認された包絡に含めない。
4. `WakeNumIter`（Requested）と`WakeIters`（Generated）、`RootWakeNodes`と`NumWakeNodes`など、同一概念の保存キー名が異なる。単純なキー一致判定は不適切。
5. 一部の`FixedWakeFlag`はRequestedにはあるが、Generatedのパーサー結果には同じフィールドが存在しない。**実際に生成されたことを確認できない項目は「要求値のみ」とする**。
6. `run.environment_json.git_commit`には`null`が存在する。現在のGitのHEADを、2026-10-31の実行時Commitとして補記してはならない。
7. `geometry_equivalence_status`は現状`UNVERIFIED`。Geometry Setに含まれた事実だけで、実機／風洞模型との形状等価性は認定できない。

---

## 2. 修正の境界・不変条件

**F-001（対象）**：表紙、1.1～1.3、2.1～2.3、3.1～3.5、4.1～4.4に自動入力ブロックを追加する。既存の13章構成、25 Test ID、既存CLIの基本的な呼出し形式は維持する。

**F-002（読み取り専用）**：報告書の生成は固定SQLiteを参照し、原本SQLite、元Run、Artifact、Study JSON、計画書、Acceptance Rule、Assessmentを変更しない。必要なSQLite接続には`PRAGMA query_only=ON`を使用できるようにする。

**F-003（母集団の固定）**：採用Run・採用Quantity・比較対象は`results_snapshot.json`の選択結果に従う。`planned_case`は予定Caseの母集団、`run`全件は履歴の母集団である。双方を同じ列や見出しで混在させない。

**F-004（自動判定の禁止）**：`SUCCEEDED`、`reviewed=true`、`geometry_membership_reviewed=true`、一致するSHA-256、数値設定の一致のみを理由に、`PASS`、`APPLICABLE`、Geometry同等性確認済み、承認済み等を出力しない。

**F-005（独立性）**：外部Validation、Cross-code、要求精度、用途承認、出典文献との物理的同等性は、記録された人間の審査結果がない限り未評価とする。

**F-006（後方互換）**：既存報告書の章番号、表の意味、図・試験別Markdownへの相対リンク、`evidence_manifest.json`の基本形式を維持する。追加フィールドは拡張として実装する。

**F-007（過去の環境）**：実行時環境と報告書再生成時環境を別の表にする。過去Runのバイナリを現在のPCに存在するバイナリから推定しない。

---

## 3. 自動入力フィールドの状態モデルと出典

### 3.1 状態と表示規則

以下は**データフィールドの埋め方を示す状態**であり、工学的な合否とは無関係とする。

| 状態コード | 定義 | 表示例 | 自動採用 |
|---|---|---|---|
| `RECORDED` | 固定DBに1個以上の原記録が存在 | `20（Requested）` | 出典付き記載のみ可 |
| `CROSSCHECKED` | 選択Runの対応する異なる記録が一致 | `20（Requested/Generated整合）` | 設定値の事実として可 |
| `CONFLICT` | 同一意味の記録が不一致、またはメタデータ矛盾 | `要確認：G1 / stage=G3` | 単一値の確定は禁止 |
| `NOT_RECORDED` | 記録がない、記録値がnull | `DB未記録` | 補完・類推禁止 |
| `NOT_EXECUTED` | 対象Caseに選択された成功Runがない | `実効値なし（成功Runなし）` | Requestedのみ別表示可 |
| `REVIEW_REQUIRED` | 記録はあるが、物理的妥当性を要審査 | `符号規約：要確認` | 原記録としてのみ可 |
| `MANUAL_DECISION` | 使用用途・許容誤差・承認等の判断 | `未承認（担当者入力）` | 自動確定禁止 |
| `NOT_APPLICABLE` | 対象外と確認できる項目 | `該当なし（理由）` | 理由と適用範囲が必要 |

同じ行で複数状態を持つ場合、個々の値に状態を付す。例：`WakeIters: 20[CROSSCHECKED] / FixedWakeFlag: false[RECORDED・要求のみ]`。Markdownの色や絵文字だけに依存せず、文字列で判読可能にする。

### 3.2 原記録の役割と選択順位

| レイヤー | 正本または一次的に参照すべき情報 | 役割 |
|---|---|---|
| Study | `study` | Study ID、計画書版／SHA、ターゲット版、DB版 |
| 計画 | `test_definition`、`planned_case`、`flight_condition`、`numerical_setting` | 何を試験する予定だったか |
| 採用解析 | `results_snapshot.cases[].run_id` と対応する `run`、`case_run_link`、`case_point` | どの実行を報告書に採用するか |
| 実行設定 | `run.requested_snapshot_json` | 当該Runに要求した設定。後日の変更の影響を受けない |
| API実効値 | `run.effective_snapshot_json` | APIで読出した実効値 |
| 生成値 | `run.generated_snapshot_json`（特に`solver_settings`） | 生成されたソルバー入力をパースした値 |
| 物理量 | 選択Runに紐づく`dataset`、`quantity` | 取得結果。採用Run以外の量を混入させない |
| 解析原本 | `artifact`、元RunのRaw（必要な場合のみ） | 出典照会・ハッシュ検証。原本改変禁止 |
| 規則と認定 | `acceptance_rule`、`comparison`、`assessment` | 存在し承認要件を満たすものだけを判定情報として参照 |
| 文書判断 | 明示的なレビュー設定ファイル／発行管理記録 | 機密区分、担当者、正式なIntended Use、承認者等 |

**出典の優先順位は不一致を黙って上書きする順位ではない。** 例えばRequested=5、Generated=10なら`Generated=10`を単独の正式値とせず、不一致を残して`CONFLICT`とし、両値とRun IDを提示する。

### 3.3 出典証明レコード（内部データ型）

各自動入力値は少なくとも次の構造を持つ。実装は`dataclass`でも辞書でもよいが、明示的な型と検証を持たせる。

```json
{
  "field_id": "4.3.wake_iterations",
  "display_value": "20",
  "unit": "iterations",
  "state": "CROSSCHECKED",
  "basis": "SELECTED_SUCCESSFUL_RUN",
  "source_values": [
    {"source": "run.requested_snapshot_json.setting.requested_json.vspaero.WakeNumIter", "value": 20, "run_id": "RUN-..."},
    {"source": "run.generated_snapshot_json.solver_settings.WakeIters", "value": 20, "run_id": "RUN-..."}
  ],
  "case_ids": ["GEO-03-G3"],
  "notes": "工学判定ではなく入力設定の一致"
}
```

`RUN-...`はデータ形式の説明であり、実在するRun IDではない。

### 3.4 数値・単位・比較規則

- ブール値は`0/1`、`false/true`の表現を正規化して比較する。`0`は「未記録」と見なさない。
- 同一の概念名と確認できる場合にのみキーを対応付ける。別物の指標は一致比較しない。
- 数値の厳密一致が保証されない浮動小数点項目は、意味・単位に基づく事前定義された絶対・相対許容差で比較する。提案既定値は `abs_tol=1e-9, rel_tol=1e-8`（単位付きの量や離散整数では個別規則を優先）。**この許容差は記録値同士の照合用であり、空力V&Vの合格基準ではない**。
- Mach、Re、角度、速度、面積、長さ、重量を無条件に相互変換しない。Reの基準長は必ず`Cref`等の定義を添える。
- レポート本文の数値表示桁と内部の生値は別にする。本文は通常6～7有効桁程度を上限候補にし、出典の全精度はCSV／SQLiteに保持する。
- 原データがnull／欠測／非有限値の場合は `NOT_RECORDED`、数値0の場合は`0`と表示する。

---

## 4. 表紙・文書管理欄の実装仕様

| Field ID | 項目 | 情報源 | 自動処理 | 残すべき判断 |
|---|---|---|---|---|
| C-01 | 対象機体・形態 | 選択されたG103A側Geometryの`physical_model_json`とCase紐付け | 機体名とG1/G2/G3/THICKなどを区別して要約。BENCHMARKとの混合を禁止 | 形状等価性 |
| C-02 | 対象ソフトウェア | `study.target_openvsp_version`、Runのバイナリ情報 | `計画対象=3.50.5`と`実行ファイルのハッシュ`を分ける | API/バイナリ実証版の確認 |
| C-03 | Git Commit／ソース | `run.environment_json.git_commit`、実行時fingerprint、再生成コードfingerprint | NULLなら`実行時Commit未記録`。再生成コードは別欄 | 推測による過去Commit補完禁止 |
| C-04 | Study/Plan/DB/Manifest | `study`、固定DBのSHA-256 | 現行の自動入力を維持 | 計画実ファイルとの同一性（未検証なら分ける） |
| C-05 | 作成者・検証者・承認者 | オプションのレビュー設定ファイル | 値が明示的に与えられたときのみ転載 | 権限・署名の妥当性 |
| C-06 | 機密区分・配布範囲 | レビュー設定ファイル | 入力がなければ`未設定（文書管理者判断）` | 自動推定禁止 |
| C-07 | 発行版・発行日 | 発行管理（レビュー設定） | 草稿は引き続き`未発行` | 自動発行・承認禁止 |

---

## 5. 第1章のフィールド仕様

### 5.1 第1.1節：審査対象と設計用途

| Field ID | 現在の項目 | 自動入力する情報 | 根拠／ガード | 結果状態 |
|---|---|---|---|---|
| 1.1-01 | Intended Use | 計画のTest FamilyとQoIを**参考候補**として併記 | 正式な使用目的はDBに確定登録なし | `MANUAL_DECISION` |
| 1.1-02 | 機体・構成 | G103A、各GEOの登録役割、部品別Thin/ThickのRequested/Effective | G1/G2 stage矛盾・GEO-04 THICK矛盾は警告 | `RECORDED/CONFLICT` |
| 1.1-03 | 飛行条件の評価範囲 | 選択成功CaseのConditionを機体・試験Family別に集計 | VERベンチマークをG103Aに含めない | `RECORDED` |
| 1.1-04 | 評価対象QoI | 選択されたQuantityのQoI／単位／Method別一覧。取得済み/計画のみを区別 | `quantity.value_status`、Datasetを考慮 | `RECORDED` |
| 1.1-05 | 主要根拠 | 成功RunがあるTest Family、Theory/WIND_TUNNEL/PANEL/NS Datasetの有無 | テスト定義の存在だけで参照取得済みとはしない | `RECORDED` |
| 1.1-06 | 結果の要約 | 69 Case、採用成功42、未成功17、未解決10、評価0など | 計画・実行・工学判定を分離 | `RECORDED` |
| 1.1-07 | 制限事項 | 未成功Case、外部Validation未登録、未承認Rule、座標系未整合、Geometry不一致を診断リスト化 | 自動列挙は事実とリスクの記載に限定 | `RECORDED/REVIEW_REQUIRED` |
| 1.1-08 | 使用総合判断 | 承認済みApplicability Assessmentがなければ`INSUFFICIENT_EVIDENCE` | ユーザー承認なしにPASSへ移行しない | `MANUAL_DECISION` |

### 5.2 第1.2節：決定事項

`DEC-01`～`DEC-05`を削除・自動承認しない。各Decisionに`候補・実績`列（または注記）を追加して、関係するCase／数値結果／Acceptance Rule／Assessmentの有無を表示する。適用例：

- `DEC-01`：GRID-01/02/03の採用成功Runとメッシュ実測量は表示可能。しかし「基準メッシュ採用」は未評価。
- `DEC-02`：WAKE-01/02/03のRequested/Generated値を列挙可能。しかし「標準Wake設定の確定」は未評価。
- `DEC-03`：成功RunのCondition群は表示可能。「全区間で使用可能」は未認定。
- `DEC-04`：STA/CTRL取得済みQuantityを表示可能。「安定微係数・舵効きの採用」は未評価。
- `DEC-05`：CD0方法A/Bの派生Dataset有無と前提は表示可能。「方法選定」は未評価。

### 5.3 第1.3節：未解決事項

`OPEN-*`という**報告書内の自動診断ID**を用意し、理由・影響・対応先を行で表示する。手動の正式Action IDに勝手に昇格させない。最低限、(a) 成功RunがないCase、(b) `geometry`内の記録矛盾、(c) 承認済みRule欠如、(d) 外部参照不足、(e) 実行版／座標系不確定、(f) 未認定領域を扱う。詳細は`tables/report_field_audit.csv`に出力する。

---

## 6. 第2章のフィールド仕様

### 6.1 第2.1節：使用シナリオと評価要求

- `test_definition.description`等から「必要候補QoI／試験群」を提示することは可能。
- `UC-01`の「使用目的・設計判断・精度・未達時影響・要求書頁」は**人間が定義する設計要求**であり、DBに直接対応する承認レコードがない。いずれも自動確定禁止。
- 外部設定ファイルに使用シナリオを定義する場合は、出典Document ID、Revision、Reviewer、承認状態と適用Studyを要求する。`PROPOSED`なら承認済みとして表示しない。

### 6.2 第2.2節：対象形状と評価段階

| 登録上の形態 | DB候補 | 事前検査 | 表示要件 |
|---|---|---|---|
| G1 | `GEO-01-G1-G` | `stage=G3`／`case_role=G1`の不一致 | 両値・Geometry ID・選択Runなしを記載 |
| G2 | `GEO-02-G2-G` | `stage=G3`／`case_role=G2`の不一致 | 両値・Geometry ID・選択Runなしを記載 |
| G3 | `GEO-03-G3-G` | 成功Runの形状Setと登録設定を照合 | Geometry IDと対応Run、構成を表示 |
| THIN/THICK | `GEO-04-THIN-G`／`GEO-04-THICK-G` | 登録されたwing=Thinと、THICK側Runの主翼Thickが不一致 | `CONFLICT`を表示。自動修復しない |

G1/G2などで成功Runがない場合は、実行値や生成メッシュの実績を捏造しない。`planned_case`および`numerical_setting`由来の**予定構成**のみを表示する。G3とTHIN/THICKの実効Geometry IDは名前・Setインデックスの対応検証ができる場合のみ部品名へ翻訳する。

### 6.3 第2.3節：包絡・除外範囲

計画条件、成功実行条件、外部参照条件、承認済み適用条件を**独立した列**にする。

| 列 | 入力元 | 規則 |
|---|---|---|
| 計画上の条件 | `planned_case`と`flight_condition`（必要に応じ明示的な計画書） | Plannedと記載。予定範囲が未確定ならその旨を記載 |
| 実行済み条件 | 選択成功Runの`flight_condition`および生成入力 | ベンチマーク／G103A／感度試験を別グループで集計 |
| 参照データ範囲 | 適格な外部`dataset`と条件 | 参照資料がゼロなら`未登録` |
| 認定済み適用範囲 | 承認済み`assessment.applicability_domain_json` | なければ`未認定` |
| 除外事項 | 計画書の除外事項＋実績の不足 | 自動算出事実と手動判断を分ける |

- **離散条件を連続した検証済み包絡と表現しない。** 例：α=0°/2°/6°という成功実績は、αの全区間が認定済みであることを意味しない。
- βについて成功Runが0°のみなら`β=0°の記録のみ`とする。横方向領域がValidationされたと主張しない。
- 舵角はエレベータ・エルロン・ラダー別に表示し、選択成功Runにない点と未成功点を混同しない。
- `reference_length_unit`と`cref`が異なる場合にはReを統合しない。
- VER系のMach=0.001、0.1、Re=2.00×10^6をG103Aの検証包絡として混ぜない。

---

## 7. 第3章のフィールド仕様

### 7.1 第3.1節：V&Vの分解

既存の分類表は保持する。各区分に計画Test数、Case数、採用成功Case数、参照・Comparison・Assessment件数を追加可能。**分類表の方法論そのものは変更しない**。Code Verification等の呼称を機械的に格上げしない。

### 7.2 第3.2節：QoI定義・符号・正規化

自動入力対象は、実際のQuantity集合、`qoi`、`component`、`variant`、`native_unit`、`canonical_unit`、`convention_id`、`normalization_method_id`、`flight_condition.sref/bref/cref/moment_x/y/z`、`axes_convention_id`。以下を厳守。

1. ネイティブ名称と正規化後の名称を区別し、単位`1`を**座標軸・符号規約が検証されたことの代替**にしない。
2. 係数の基準寸法は同一形状・同一Conditionの値と結合する。G103AのSref等をHershey Bar／Ellipsoidに転記しない。
3. `CLi/CDi`、`CLiw/CDiw`、`CMix/CMiy/CMiz`の物理的意味は、根拠のあるコード規約／マニュアル確認がない限り`要技術確認`と記載。
4. 動的・静的微係数、操舵微係数について、`rad^-1`等の変換や符号変更を根拠なしに行わない。
5. ベンチマーク、解析条件、微係数算出方法別に、異なるNorm／Methodを統合しない。

### 7.3 第3.3節：誤差・不確かさ

評価式の説明を保持し、`comparison`テーブルに実際の結果がある場合のみ差分・不確かさを追加する。今回の固定DBはComparisonが0なので、**比較実績なし**と明記する。理論式が書かれていることから、不確かさ推算が実施済みと見なさない。

### 7.4 第3.4節：Ruleの適用と承認状態

`acceptance_rule`が0件なら「Rule未登録」と明示し、現在の空白行は削除または状態行へ置換する。Ruleが登録されている場合でも`approval`、適用QoI、Test、条件領域、Thresholdの意味が確認できない場合は承認済みとして表示しない。`assessment`のNumerical/Engineering/Applicabilityを分離し、相互に代用しない。

### 7.5 第3.5節：Validation資料

`dataset.kind`、`usage_class`、`evidence_class`、`metadata_json`の登録有無を自動一覧化する。今回の固定DBには`WIND_TUNNEL`、`NS`、`PANEL`、`THEORY`のDatasetが存在せず、風洞Validation等の正式参照比較は**未登録**である。`EXTERNAL_CD0`は風洞Validationと同義ではない。

---

## 8. 第4章のフィールド仕様（最優先）

### 8.1 第4.1節：実行ソフトウェアと再現性

| Field ID | 報告項目 | データ出典 | 合成・判定規則 |
|---|---|---|---|
| 4.1-01 | Study対象OpenVSP版 | `study.target_openvsp_version` | `計画対象:3.50.5`と記載。実行時実測版ではない |
| 4.1-02 | VSPAERO実行バイナリ | 採用Runの`environment_json.binary_path`、`binary_sha256` | 全Run同一か検査、差異があれば複数行。ハッシュの一意性を版番号と同義にしない |
| 4.1-03 | OpenVSP API版／ライブラリ | 実行時に保存したArtifact・API診断の明示値 | 現在のDBに該当情報がなければ`未記録` |
| 4.1-04 | 実行環境 | `run.environment_json.ncpu`、Run Snapshot | 記録のある実行条件のみ出力。Python/OS/CPUはなければ未記録 |
| 4.1-05 | 実行時プログラムの版 | `run.environment_json.git_commit`、`execution_source_sha256`、`runner_sha256` | Commit=nullを隠さない。SHA-256は別欄で記載 |
| 4.1-06 | 再生成プログラムの版 | `results_snapshot`の可視化指紋、Manifest | 実行時コードとは別枠で明示 |
| 4.1-07 | DBスキーマ版 | `study.schema_version` | `1.0`を**スキーマ版**と表示。「VVDBプログラム版」と誤記しない |
| 4.1-08 | 計画・Studyハッシュ | `study.plan_sha256/config_sha256` | 実ファイルとの照合結果とDB内記録を区別 |
| 4.1-09 | 固定SQLite SHA-256 | `database_snapshot.sqlite`、Manifest | 現行機能を維持 |

既存60 Runで`binary_sha256`が一致していても、それだけでOpenVSP API・GUIとソルバーがすべて3.50.5と保証されたとは記載しない。

### 8.2 第4.2節：Geometryと最終メッシュ

**情報の表示粒度**：Geometry表（役割・計画）、選択成功Run別のSet実績表、メッシュ情報表の3種類に分離する。

- Geometry表：`geometry_id`、`stage`、`representation_json.case_role`、`physical_model_json`、`source_vsp3_sha256`、`geometry_equivalence_status`。
- Set実績表：採用Runの`effective_snapshot_json.geometry_sets`、`GeomSet`、`ThinGeomSet`と、`numerical_setting.requested_json.geometry_membership`。ID→部品名の変換は追跡可能な対応表で確認済みの場合に限る。
- メッシュ表：Run原本のメッシュ証拠から、存在確認できた`triangles_total`等とRun IDを列記。可能なら`quantity`／診断JSONと照合する。入力Tess値と完成メッシュ要素数は別の列で表示する。
- `GEO-01/02`のStage矛盾、`GEO-04-THICK`のwing表記矛盾は警告出力。元DBを修正しない。
- 成功RunがないGeometryは`予定構成のみ・実効Setは未確認`と表示する。

### 8.3 第4.3節：代表数値設定

列を`項目／計画候補／要求値／API実効値／生成ソルバー値／整合状態／Run・Case証拠`へ変更する。**基準Runは明示的に選択されたG103A代表Caseを使い、WAKE等の感度試験は別表に分離**する。

| 同一概念 | Requested | Effective | Generated | 注意 |
|---|---|---|---|---|
| Wake反復数 | `requested_json.vspaero.WakeNumIter` | 存在時のみ | `generated_snapshot_json.solver_settings.WakeIters` | 代表と感度試験を区別 |
| Wake節点数 | `requested_json.vspaero.RootWakeNodes` | 存在時のみ | `solver_settings.NumWakeNodes` | 代表と感度試験を区別 |
| Implicit Wake | `requested_json.vspaero.ImplicitWake` | `effective_snapshot_json.ImplicitWake` | `solver_settings.ImplicitWake` | boolean化して一致検査 |
| Fixed Wake | `requested_json.vspaero.FixedWakeFlag` | 記録がある場合のみ | 現行パーサーで対応フィールドなし | `要求値のみ`と表示 |
| Symmetry | `requested_json.vspaero.Symmetry` | `effective_snapshot_json.Symmetry` | `solver_settings.Symmetry` | 0/1⇔bool正規化 |
| Reference Flag | 設定に記録ある場合 | `effective_snapshot_json.RefFlag` | 確実な同等フィールドがある場合のみ | `RefFlag=0`の意味は仕様に基づき注記 |
| Geom Set | `thin_set_index`／`thick_set_index` | `ThinGeomSet`／`GeomSet` | — | -1は未設定ではなく選択なし |
| 基準寸法 | `flight_condition.sref/bref/cref` | `Sref/bref/cref` | `Sref/Bref/Cref` | 単位とベンチマーク区分必須 |
| モーメント基準点 | `flight_condition.moment_x/y/z` | `Xcg/Ycg/Zcg` | `X_cg/Y_cg/Z_cg` | 幾何基準・座標系を表示 |
| 主翼`Tess_W` | `requested_json.geom_parms`に含む場合 | `effective_snapshot_json.geom_parms`に該当する場合 | — | 値がなければ`記録なし`。計画書31で自動補完不可 |
| `SectTess_U` | 同上 | 同上 | — | 各セクション単位で識別、単一値への集約禁止 |
| メッシュ要素数 | — | — | 直接の生成実績がある場合 | 生成後の診断値でありTess値とは別 |

**基準Caseの選定規則**：`report_profile.json`の`baseline_case_id`または明示的な選択設定を最優先する。指定がなければ「基準Case未確定」として、`GEO-03-G3`等の候補を提案表示するにとどめる。`最新成功Run`や`役割名BASE`だけで、採用すべき標準設定の認定をしてはならない。

### 8.4 第4.4節：飛行条件・操舵・トリム

- 選択成功Caseを`condition_id`でグループ化し、G103A基準・感度解析・VERベンチマークの3群を分ける。
- 表示列：`Test ID/Case ID/Run ID`、`Condition ID`、Mach、Re(Cref)、α、β、エレベータ、エルロン、ラダー、`trim_mode`、Sref/Bref/Cref、基準点、軸規約ID、`source_state`。
- 同一数値飛行条件でも`condition_id`が違う場合、表の重複表示を抑えるときは異なるIDの一覧を脚注に残す。
- G103A基準Condition `FC-G103A-ISA-SL-105KMH-A2-20261031-REVIEWED`には、Mach 0.0857102、Re(Cref) 2.12×10^6、α=2°、β=0°、Sref 17.88185312625759 m²、Bref 17.536552870371967 m、Cref 1.0624862547003904 m、基準点(2.87,0,0.164) mの原記録がある。この値を**該当Condition以外へ無断転記しない**。
- 舵角の符号は`controls_deg_json`から出すが、物理偏向方向の確認記録（Gainなど）と結び付く場合も、別途承認状態を表示する。
- 有効なConditionに設定があっても、採用成功Runがなければ`実行実績`としては計上しない。
- `effective_snapshot_json.sweep`に複数計算点があれば、Report上のCase/Pointと適切に結合し、1 Run=1 Conditionと決めつけない。

---

## 9. 実装構成・API詳細

### 9.1 推奨ファイル変更

| ファイル | 変更案 | 備考 |
|---|---|---|
| `src/openvsppy/vv/VVReport.py` | 既存`write_evaluation_report()`に1～4章のMarker置換を追加 | エントリポイントを維持。400行級の巨大な単一関数をさらに肥大化させない |
| `src/openvsppy/vv/VVReportFacts.py`（新規提案） | 固定DBと選択RunからField/Provenanceを抽出し、矛盾診断する純粋・読取専用処理 | 責務分割を行う場合のみ新設。過剰な細分化は不要 |
| `docs/vv/templates/evaluation_report.md` | 1.1～4.4の固定プレースホルダーを明示的な`VV:AUTO`マーカーに差替え | 13章は維持 |
| `docs/vv/templates/report_profile.json`（新規提案） | 基準Case・用途等を明示する任意のレビュー用入力雛形 | 原DBを上書きしない |
| `tests/test_vv_report_revision.py` | 既存の主要テストに新仕様の回帰テストを追加 | Windowsの`read_text(encoding='utf-8')`徹底 |
| `tests/test_vv_report_facts.py`（新規提案） | Case結合・設定不一致・欠測・工学承認禁止の詳細テスト | 必要なら分割 |
| `Regenerate-OpenVSP-VV-Report.py`（配布スクリプト） | `--report-profile`などを任意引数として追加可能 | 既存`--repo --study-root --db`は変更しない |

### 9.2 処理フロー

1. **固定データ読込**：過去DBをSQLite backupで複製、整合性検査、選択Run固定を完了する（既存の再生成処理を利用）。
2. **結合検証**：選択CaseのRun IDがDBの`case_run_link`・`case_point`と整合することを確認。
3. **Fact抽出**：Study、計画条件、採用実行、Quantity、判定規則、Geometryを別々に抽出する。
4. **整合性診断**：同一概念のRequested／Effective／Generated、Geometry Stage／Role、基準寸法／座標規約を照合する。
5. **出力モデル作成**：章ごとに短い表示用値と詳細出典ID、状態を生成する。
6. **テンプレート差込み**：厳格なMarker個数チェック後、表紙・1～4章へ置換。その他章の既存レンダリングを保持。
7. **監査出力**：`tables/report_field_audit.csv`と`report_field_provenance.json`を生成。
8. **品質ゲート**：空欄・矛盾・参照・状態の検査。`CONFLICT`があっても草稿は出力可。ただし`CROSSCHECKED`へ自動変換してはならない。
9. **Manifest**：新たなファイルを既存のEvidence Manifestに追加。全成果物は同じ固定DBに紐づける。

### 9.3 テンプレートのMarker仕様

`<!-- VV:CHAPTER:1 -->`など既存Markerは維持し、その下位に`<!-- VV:AUTO:1.1 -->`、`<!-- VV:AUTO:1.2 -->`、`<!-- VV:AUTO:1.3 -->`、同様に`2.1`～`4.4`を追加する。各Markerは**ちょうど1件**。例えば：

```markdown
### 4.3 代表数値設定

<!-- VV:AUTO:4.3 -->
```

ただし、既存の静的表と自動表を重複掲載しない。段階的移行の場合、まず自動表の生成結果を確かめた後、旧固定表を置換する。

### 9.4 擬似コード（実装意図）

```python
selected = validate_selected_runs(frozen_db, results_snapshot)
report_facts = collect_report_facts(frozen_db, selected, optional_profile)
report_facts = audit_cross_source_consistency(report_facts)
chapter_blocks = render_autofill_sections(report_facts)
report_document = render_existing_13_chapter_template(chapter_blocks, existing_chapters)
write_utf8(report_document)
write_csv(report_facts.audit, 'tables/report_field_audit.csv')
write_json(report_facts.provenance, 'report_field_provenance.json')
```

`collect_report_facts()`と`audit_cross_source_consistency()`は、ソルバー起動・SQLite更新・ファイル削除・Git状態変更をしてはならない。

### 9.5 任意レビュー設定ファイル

```json
{
  "schema_version": "REPORT_PROFILE_V1",
  "study_id": "VSP3505-VV-001",
  "baseline_case_id": "GEO-03-G3",
  "intended_uses": [],
  "document_control": {
    "classification": null,
    "author": null,
    "reviewer": null,
    "approver": null,
    "approval_status": "DRAFT"
  }
}
```

`baseline_case_id`の指定は代表情報の選択を意味し、採用基準メッシュのEngineering Acceptanceを意味しない。レビュー設定はReport ManifestにSHA-256と適用日時を保存する。レビュー設定がない場合は原則として機械的な候補表示にとどめる。

---

## 10. 視覚設計・Markdown要件

- 第1章は**Executive Summaryとして1画面程度で重要事実と未解決事項を認識できること**を目指す。列幅の大きいハッシュは短縮表記＋証拠一覧参照とし、完全SHA-256は別表に残す。
- 第2章は計画・実行実績・Validation参照・認定範囲を4列に分ける。未認定の値をグラフの連続範囲として描画しない。
- 第3章は「既存の方法論」「登録済みの客観情報」「技術審査が必要な項目」を視覚的に分ける。
- 第4章は表が肥大化しやすいため、本文に**基準Case概要＋差異・矛盾一覧**、`test_reports/`やCSVに**全Run・全設定**を置く。25条件をすべて長い表で本文に積まない。
- GFM表の各行は列数を保証し、`|`、バックスラッシュ、改行、識別子をエスケープする。
- 表中の長いWindows絶対パスやJSONを直接表示せず、Artifact ID・短いパス・SHA256の先頭などを用いる。完全値はAudit JSON/CSVに置く。
- Markdownの色・HTML独自タグに依存しない。GitHub／VS Codeの通常プレビューで表、内部リンク、目次を読めること。
- 現行の「発行済み報告書ではない」という注意書きを維持し、`CONFLICT`がある場合は第1章にも未解決課題を明示する。

---

## 11. テスト仕様

各IDは新規に追加または既存テストを拡張する。期待結果に成功条件と禁止事項を含める。

| Test ID | 試験内容 | 合格基準 |
|---|---|---|
| AF-01 | 実固定DBで表紙と第1～4章を生成 | 13章構造、25試験、既存リンクを保持 |
| AF-02 | 69 Case/42採用成功と60 Run/43成功Runの区別 | 数値を混同しない |
| AF-03 | `quantity`全件1,984と採用1,937の区別 | 異なる母集団を混入しない |
| AF-04 | VERベンチマークとG103Aを分離 | G103AにMach=0.001/0.1を誤転記しない |
| AF-05 | G1/G2 `stage=G3`と`case_role`を照合 | `CONFLICT`を表示しStageを自動修正しない |
| AF-06 | GEO-04-THICKの記録不一致 | WingのThin/Thickを根拠なく単独確定しない |
| AF-07 | Requested=20、Generated=20 | `CROSSCHECKED`と出典Runを出力 |
| AF-08 | Requested=5、Generated=10の人工データ | `CONFLICT`として両値を表示 |
| AF-09 | `FixedWakeFlag`の生成値なし | `要求値のみ`、`CROSSCHECKED`にしない |
| AF-10 | 空文字・null・0・Falseを混在 | 欠測と数値0/Falseを区別 |
| AF-11 | `Tess_W`実測データ欠測 | 計画31を「実際値」へ転記しない |
| AF-12 | Caseに成功Runがない | Planned値のみ表示し、実行済みとしない |
| AF-13 | 一部CaseでRunを明示選択 | 選択Runを変更した場合だけ該当章の実績が変化 |
| AF-14 | 条件が同値、Condition IDが別 | 重複抑制しても全IDを追跡可能 |
| AF-15 | Rule/Comparison/Assessment各0件 | PASS・使用認定・不確かさ値を捏造しない |
| AF-16 | 人工の承認済みRuleと未承認Rule | 承認状態を正しく区別 |
| AF-17 | 元SQLiteハッシュ検査 | 再生成前後で同一。元Runも書込みなし |
| AF-18 | Markdown構文検査 | GFM表の列数・リンク・内部Markerの完全性を確認 |
| AF-19 | 出典ID照合 | 出力した全Run/Case/Geometry/Condition IDは固定DBに存在 |
| AF-20 | JSON／CSV／Manifestのハッシュ照合 | 生成物が一貫した固定DBと結び付く |
| AF-21 | Windows CP932環境でpytest | 全Markdown読み書きでUTF-8明示。UnicodeDecodeErrorがない |
| AF-22 | オフライン環境・APIなし | OpenVSPをインポート／起動せずReport生成成功 |
| AF-23 | `git_commit=null` | `未記録`と表示しHEADを誤用しない |
| AF-24 | 実行時と報告時のコード指紋が異なる | 2つのフィンガープリントを区別 |
| AF-25 | 任意レビュー設定未提供 | 手動欄は未承認・未設定を維持 |
| AF-26 | 任意レビュー設定のStudy ID不一致 | 例外で停止し、元データは保護 |
| AF-27 | 設定値整合の浮動小数点近似 | 記録間照合許容差とV&V受入基準を混同しない |
| AF-28 | `case_point`複数点 | 異なるCase／Pointの数値を誤結合しない |

**検証の順番**：新規単体試験 → 現行149件相当のオフライン回帰試験（環境差に注意）→ 実固定DBからのレポート再生成 → Markdownプレビューでの目視確認。OpenVSP 3.50.5の再計算はこの仕様の受入要件ではない。

---

## 12. 受入基準（Acceptance Criteria）

| ID | 完了条件 |
|---|---|
| AC-01 | 第1～4章の全自動対象フィールドで、記録がある場合は出典付きの値、欠測時は理由付きの状態が表示される |
| AC-02 | 既存の各章の静的な意味を保持し、同一表題の重複表を残さない |
| AC-03 | 自動入力値は`report_field_provenance.json`からDBカラム／JSONキー／Run IDに追跡できる |
| AC-04 | G1/G2・GEO-04-THICKの既知不整合が表示され、誤った一意の形状表現にならない |
| AC-05 | 基準Caseと感度試験、G103Aとベンチマーク、予定値と実効値が明確に分離される |
| AC-06 | `SUCCEEDED`だけでNumerical/Engineering/ApplicabilityのPASSを作成しない |
| AC-07 | 既存SQLite・Raw・幾何モデル・計画書のバイト列／SHA-256が変わらない |
| AC-08 | 13章構成、25 Test ID、相対図リンク、Manifest、既存CLIの後方互換を確認する |
| AC-09 | `AF-01`～`AF-28`の対象テストが成功する。外部依存テストは実行条件を記録する |
| AC-10 | `evaluation_report.md`をMarkdownプレビューで開いて、表が壊れず、長い識別子で本文の視線を妨げない |
| AC-11 | 旧レポートとの比較表で、追加した事実・残した未承認事項・表示した不整合が一覧できる |

---

## 13. 実装順序と変更を避ける対象

| フェーズ | 作業 | 主成果物 |
|---|---|---|
| P0（監査） | 現行固定DBのField inventory／既知矛盾をテスト化 | `tests/test_vv_report_facts.py` |
| P1（中核） | Fact抽出・選択Run結合・出典管理・状態分類 | `VVReportFacts.py`相当 |
| P2（第4章） | 4.1～4.4のRequested/Effective/Generated表と基準Case表 | Markdown・Audit JSON |
| P3（第1～3章） | 範囲・QoI・試験実績・欠測一覧・Rule状態を自動入力 | Markdown第1～3章 |
| P4（表示） | 長い表を分割、脚注・出典リンク・GFMプレビュー確認 | 新テンプレート |
| P5（QA） | 全オフライン回帰テスト、実データ再生成、Manifest、版管理 | 修正版ZIPと検証ログ |

**変更しない対象**：`vv_study`の元DB、原本Run、解析モデル、Geometry ID、Test ID、既存判定テーブルの内容、ソルバーバイナリ。必要な警告や説明は**出力側で解決**し、データ修正が必要なら別途DB移行・レビュー手続として扱う。

---

## 14. 現時点で人間が決める事項

1. **正式なIntended Use**：どの用途（概念設計・性能・安定微係数等）に、どの精度で使いたいか。
2. **基準Caseの採用**：報告書で代表例に使うCase ID（提案候補GEO-03-G3）。これはNumerical Acceptanceの承認と別。
3. **技術規約**：VSPAEROの係数・符号・座標系の外部比較用規約。
4. **形状メタデータの不整合の扱い**：G1/G2、GEO-04-THICKの記録差異を誰がどの証拠で審査するか。
5. **Rule・Accuracy要求**：工学的PASS/FAILの閾値と承認者。
6. **文書統制**：機密区分、作成・確認・承認、発行版と日付。

これらの回答がなくても、**草稿の自動入力機能自体は実装可能**。未確定事項は適切に`MANUAL_DECISION`などと記録し、処理停止の原因にしない。

---

## 15. 最終技術方針

**自動化できるのは「証拠として記録された事実の取得・分類・整合確認・明瞭な表示」であり、工学的判断や過去の欠測データの補完ではない。**

特に第4章では、Requested／Effective／GeneratedとCase選択を分けた表を生成することが最も有効である。第1章はそれを要約し、第2章は検証した離散条件と未認定領域を分け、第3章は定義・出典の実績と未承認基準を分ける。こうすることで、現行SQLiteスキーマや解析原本を変更せず、レポートの情報価値と人間の審査効率を高められる。
