# VSPAERO V&V評価報告書とSQLite項目の対応表

**対象DDL**：`reference/schema.sql`（計画書v1.1／DB仕様v1.0.1に対応する17テーブル）

**留意**：本表はDBの既存カラム名に基づく報告書生成仕様である。現在のDBプログラムv0.1が全項目の自動生成・集計に対応済みという意味ではない。DB正本への書込みは既存VVDatabaseのトランザクション管理を用い、報告書ジェネレータは原則読み取り専用とする。

## 1. 識別子・データフロー

| 粒度 | テーブル | 主キー | 主なリンク |
|---|---|---|---|
| Study | `study` | `study_id` | `plan_version`, `plan_sha256`, `config_sha256` |
| 試験定義 | `test_definition` | `test_id` | `study_id`、試験25件 |
| 予定解析 | `planned_case` | `case_id` | `test_id`,`geometry_id`,`condition_id`,`setting_id` |
| 実行 | `run` | `run_id` | `geometry_id`,`condition_id`,`setting_id` |
| 計算点 | `case_point` | `point_id` | `run_id`,`target_condition_id` |
| Caseと結果の結合 | `case_run_link` | `(case_id,point_id)` | `planned_case`と`case_point` |
| Dataset | `dataset` | `dataset_id` | `point_id`,`condition_id`,`geometry_id` |
| Quantity | `quantity` | `quantity_id` | `dataset_id`, `source_artifact_id` |
| 比較 | `comparison` | `comparison_id` | `test_id`,2つのQuantity、`rule_id` |
| 評価 | `assessment` | `assessment_id` | `test_id`,`case_id`、用途領域 |
| 原本・生成物 | `artifact` | `artifact_id` | `run_id`,`dataset_id`,`relative_path`,`sha256` |
| 実行監査 | `run_event` | `event_id` | `run_id`と状態遷移 |

**重要な制約**：`assessment`は**全てが必ず1件のRunに直接結び付く構造ではない**。試験・QoI・適用領域の根拠として`rule_ids_json`と`evidence_ids_json`を利用する。`comparison`も1つのPredicted／Reference値を対象とするため、数条件にまたがる統合の検証結果は報告書作成時に明示的に集約する。

## 2. 報告書の各章との対応

| 報告書 | 正本テーブル・カラム | 導出処理／注意事項 | 現行v0.1の範囲 |
|---|---|---|---|
| 表紙・対象Study | `study.study_id/name/plan_version/plan_sha256/schema_version` | 抽出時の計画書とハッシュ照合 | 主な登録可能 |
| 1章 Executive Summary | `assessment`, `comparison`, `acceptance_rule` | 使用認定は人間のレビューと正式承認 | 要手動レビュー |
| 2章 用途・領域 | `flight_condition`,`geometry`,`assessment.applicability_domain_json` | Intended Use全文はDB外文書にも記録 | 部分的 |
| 3章 評価規則 | `acceptance_rule.policy_json/approval/reviewer/approved_utc` | 暫定Ruleを正式PASSに使用しない | 部分的 |
| 4章 Geometry | `geometry.physical_model_json/representation_json/source_vsp3_sha256` | 最終メッシュの品質は別証拠が必要 | 登録可能／実機未検証 |
| 4章 Flow | `flight_condition`のMach／Re／α／β／controls／基準量等 | 軸・座標・単位を必ず合わせる | 登録可能 |
| 4章 Numerics | `numerical_setting.requested_json`; `run.*snapshot_json` | Requested／Effective／Generatedを区別 | 部分的／実機未検証 |
| 5章 Test catalog | `test_definition`,`planned_case`,`case_run_link`,`run` | Testの実施状態は**派生集計**。run statusの直訳ではない | 25件登録済み |
| 5章 逸脱 | `run_event`,`planned_case`,`run`＋レビュー文書 | 逸脱管理用の専用テーブルは現行DDLにない | 報告書に手動記録 |
| 6章 VER | `dataset.kind='THEORY'`, `quantity`, `comparison` | ベンチマーク理論・VSPAERO比較を確認 | 取込基盤／自動評価未完成 |
| 7章 GRID／WAKE | `quantity`,`planned_case`,`comparison`,`assessment` | 実メッシュ・Wake感度は結果値と追加証拠を区別 | 3点暫定スクリーニングあり |
| 7章 GEO | `geometry`,`dataset`,`quantity` | モデル化差と数値誤差を別に解釈 | 自動差分は未完成 |
| 8章 XCD | `dataset.kind IN ('PANEL','NS')`,`quantity`,`comparison` | `DIRECT/CONVERTED/INTERPOLATED/ILLUSTRATIVE` | CSV取込・量の比較あり |
| 9章 VAL | `dataset.kind='WIND_TUNNEL'`, `usage_class`, `comparison` | 不確かさは定義・水準・計算根拠を検証 | 正式UQ・認定未完成 |
| 10章 STA/CTRL | `quantity.qoi/variant`,`quantity_dependency`, `dataset.metadata_json` | 定常FD／Adjoint／PQR・舵面Gain・単位の分離 | `.stab`原値取込中心 |
| 11章 DRAG | `dataset.kind IN ('EXTERNAL_CD0','DERIVED')`, `quantity`, `quantity_dependency` | 外部CD0方法A/Bは元値と別Dataset | 派生計算は実装済み |
| 12章 ACC | `assessment.*`, `acceptance_rule.*` | 認定は自動生成しない | 手動レビュー必要 |
| 13章 審査・承認 | `assessment.reviewer/approved_utc`＋発行文書 | 報告書全体の署名・承認は別管理 | 文書管理で対応 |
| 付録B 原本 | `artifact.relative_path/sha256/size_bytes/integrity_status` | Report Snapshot単位にmanifest作成 | 原本SHA-256照合あり |

## 3. 主な実データカラム

### 3.1 `run`の要求・実効・生成スナップショット

- `run.requested_snapshot_json`：要求した設定・環境
- `run.effective_snapshot_json`：OpenVSP APIから読んだ実効値
- `run.generated_snapshot_json`：`.vspaero`や生成メッシュ等、確認可能な実設定
- `run.binary_sha256`, `run.runner_sha256`：実行ファイル・実行アダプタの同定
- `run.status`, `run.created_utc`, `run.started_utc`, `run.ended_utc`：実行状況・時刻
- `run.fingerprint_sha256`, `run.environment_json`：比較・再現の環境情報

### 3.2 `quantity`の原値／正規化値

- `quantity.qoi`, `component`, `variant`：区別すべき物理量と算出方法
- `native_name/value/unit`：ネイティブ出力の値・名称
- `canonical_value/unit`：比較用に正規化した値・単位
- `convention_id`, `normalization_method_id`：変換規約
- `source_artifact_id`, `source_locator`：出典ファイルと列・行・抽出箇所
- `value_status`：`AVAILABLE`等。未計算・不正な値を0で埋めない。

### 3.3 `comparison`の不確かさ

- `delta_value`, `relative_delta`：同一の単位・規約に整合後の差
- `standard_uncertainty`, `expanded_uncertainty`, `uncertainty_metadata_json`：不確かさの値・定義
- `rule_id`, `verdict`, `rationale`：判定の規則・結果・理由
- `comparison_kind`：直接／換算／補間／参考の比較を区別

### 3.4 `assessment`の用途判断

- `numerical_status`, `engineering_status`：各評価の判定
- `disposition`：`APPLICABLE`／`CONDITIONAL`／`NOT_APPLICABLE`／`INSUFFICIENT_EVIDENCE`
- `applicability_domain_json`：用途・飛行条件・QoIの範囲
- `rule_ids_json`, `evidence_ids_json`, `restriction_json`：根拠と制限
- `reviewer`, `approved_utc`, `rationale`：技術承認と説明

## 4. DB Snapshot・発行版の推奨管理

本DDLには`report`、`report_snapshot`、`approval_log`、`deviation`という専用テーブルは**存在しない**。そのため報告書の発行にあたり、SQLiteを編集せずに次の方法を用いる。

1. 一貫した時点でSQLite DBのバックアップを取得し、SHA-256を記録する（実行中DBを単純コピーせずSQLite backup APIを推奨）。
2. 当該DBへ関連する元ファイル・設定・実行バイナリのArtifact一覧／ハッシュを書き出す。
3. Report ID、報告書版、発行時UTC、DB snapshotハッシュ、使用した後処理・図表生成ソフトウェアの版を`report_manifest.json`へ保管する。
4. 報告書の承認サイン・議事録等を報告書の付録として保持する。
5. DBを更新しても過去の発行版Snapshotを置換しない。

## 5. SQL照合の例（読み取り専用）

```sql
-- 25件のTest IDが揃っているか
SELECT test_id, category, requirement
FROM test_definition
ORDER BY test_id;

-- 各Test IDの予定ケース数
SELECT t.test_id, COUNT(c.case_id) AS planned_cases
FROM test_definition AS t
LEFT JOIN planned_case AS c ON c.test_id=t.test_id
GROUP BY t.test_id ORDER BY t.test_id;

-- Runのステータス一覧（SUCCEEDED はV&V PASSではない）
SELECT r.run_id, r.status, r.geometry_id, r.condition_id, r.setting_id,
       r.created_utc, r.ended_utc
FROM run AS r ORDER BY r.created_utc, r.run_id;

-- 指定したTest IDの比較結果
SELECT c.comparison_id, c.test_id, c.comparison_kind,
       c.delta_value, c.verdict, c.rule_id, c.rationale
FROM comparison AS c WHERE c.test_id = :test_id;

-- 使用可否の判断とレビュー者
SELECT a.assessment_id, a.test_id, a.qoi,
       a.numerical_status, a.engineering_status, a.disposition,
       a.reviewer, a.approved_utc
FROM assessment AS a ORDER BY a.created_utc;
```

**注意**：上記クエリはカラムの存在と関連付けを示すものであり、評価方法の物理的妥当性を自動保証するものではない。複数の試験・飛行条件点が存在する場合、意図に応じたJoin／フィルタを明示する。
