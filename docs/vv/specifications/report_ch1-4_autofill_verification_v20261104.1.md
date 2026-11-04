# 第1～4章自動入力：実装・検証記録（v20261104.1）

対象：`OpenVSP-main_v20261104.1`。2026-11-04。検証はLinux/Python 3.13、非ネイティブ実行。

## 検証範囲

- `python -m pytest -q -rs`: **174 passed, 13 deselected**（OpenVSP/VSPAERO/GUIを要するテストは未実行）。
- 添付の過去`vv_study`と固定DBを用い、`tools/regenerate_vv_report.py`を実行（通常版と明示profile版）。両方正常終了。
- 69計画Case、42採用成功Case、1,937採用Quantity、35図、25試験詳細Markdownを再生成。
- 本文13章、Markdownの表構文（markdown-it-py）、125件の相対リンクに欠損0。Manifest 87ファイルのSHA-256一致。
- 元SQLite SHA-256 `a5637d79db70eb48c7836a4fbd44877ff6e5cf9746d13881148f2880be964166` は前後一致。
- `report_profile.json`から`GEO-03-G3`を代表表示に指定して再生成済み。認定・工学PASSは付与していない。

## 詳細仕様 受入基準

| ID | 検証結果 | 証拠・残課題 |
|---|---|---|
| AC-01 | 部分充足 | 採用Run事実・欠測・区分・出典を自動生成。正式用途、要求精度等は引き続き手動。 |
| AC-02 | 充足 | 13章構成・既存方法論を維持し、対象の旧固定表を動的表へ差替え。 |
| AC-03 | 充足 | `report_field_provenance.json`、各監査CSV。 |
| AC-04 | 充足 | G1/G2 Stage、GEO-04-THICKの不一致を `CONFLICT` として表示。 |
| AC-05 | 充足 | BENCHMARK/G103A、採用成功/予定、代表/感度を区別。 |
| AC-06 | 充足 | Engineering／Applicabilityの判定は自動承認しない。 |
| AC-07 | 充足 | 元DBのSHA-256一致、解析原本は読み取りのみ。 |
| AC-08 | 充足 | 13章・25試験・既存CLI・Manifest・相対図リンク。 |
| AC-09 | 部分充足 | オフライン174件成功。全28件を独立テストで網羅していない（下記参照）。ネイティブ13件未実行。 |
| AC-10 | 部分充足 | GFMパーサーによる構文検証成功。WindowsのVS Code/GitHub実画面での目視確認は未実施。 |
| AC-11 | 部分充足 | 本文とAuditに自動事実・未承認・CONFLICTを列挙。旧版との差分表の自動生成は未実装。 |

## AF-01～28の観点

| テスト観点 | 結果 | 備考 |
|---|---|---|
| AF-01～04 | 確認 | 実データの章立て、母集団の区別、BENCHMARK/G103Aの分離。 |
| AF-05～06 | 確認 | 既知Geometry不一致を表示し、DBを修正しない。 |
| AF-07～11 | 確認 | Requested/Effective/Generated整合、不一致・欠測・False/0・Tess非補完。 |
| AF-12 | 確認 | 採用成功RunがないCaseから実効値を作らない。 |
| AF-13～14 | 部分確認 | 選択Run・Condition IDの結合処理あり。個別の多Run切替／重複IDテストは未追加。 |
| AF-15～16 | 確認 | Rule／Comparison／Assessment欠如と、承認・仮Ruleの区別。 |
| AF-17～20 | 確認 | DB不変、Markdown表・リンク・Manifest・出典。 |
| AF-21 | 部分確認 | UTF-8明示とLinuxでテストPASS。Windows cp932での再実行は未確認。 |
| AF-22～27 | 確認 | ネイティブAPI不要・Git欠測・指紋区別・Profile・Study不一致・数値照合許容差。 |
| AF-28 | 部分確認 | Case/Run/Point/Dataset/Quantity結合を検査。複数Pointの専用合成テスト未実施。 |

## 実装した主要ファイル

- `src/openvsppy/vv/VVReportFacts.py`（読取専用のFact/出典抽出、整合診断、監査出力）
- `src/openvsppy/vv/VVReport.py`（Marker置換、出力追加、Profile受取り）
- `docs/vv/templates/evaluation_report.md`（第1～4章AUTO Markerと重複表の整理）
- `docs/vv/templates/report_profile.json`（任意の明示的レビュー用入力）
- `tools/regenerate_vv_report.py`（旧DB読み取り専用の完全再生成）
- `tests/test_vv_report_autofill.py`（新規オフライン試験）

## 残る技術上の注意

- `CROSSCHECKED`は保存された設定記録の一致であって、ソルバーの物理的正しさやV&VのPASSを保証しない。
- G1/G2・THICK構成の不整合は、報告書に警告表示するだけで元のStudy記録は変更しない。
- `report_profile`は表示用代表Caseと文書管理草稿の情報であり、設計用途・数値精度・使用認定の正式な承認記録ではない。
- Windows上のOpenVSP 3.50.5新規ソルバー実行、Windows Markdownプレビューの目視確認は行っていない。
