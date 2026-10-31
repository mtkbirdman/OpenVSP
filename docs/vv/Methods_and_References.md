# VSPAERO V&V評価報告書：参考文献・採用手法対応

## 1. 資料の分類と適用上の注意

一次資料を優先して調査し、以下の役割に分けた。**本プロジェクトは公開規格の全部を実装することも、規格準拠証明を取得することも目的としない。** 有償規格は公開要約・文献情報により確認できる範囲があり、規格本文の条文適合性は別途正式に審査する。

| Reference ID | 文献／資料 | 本プロジェクトへの採用内容 | 限界・注意点 | 報告書章 |
|---|---|---|---|---|
| R-01 | ASME V&V 20-2009 (R2021) | 特定Validation点・QoIにおける解析と実験の差・不確かさ | 包絡内の内挿／外挿は原則、規格の点評価の直接対象外。一般的PASS閾値なし | 3、9、12 |
| R-02 | AIAA G-077-1998 (R2002) | Verification／Validation用語、CFD評価の段階性 | ソルバー全体の無条件な認証ではない | 3、5～9 |
| R-03 | NASA NPARC Documentation Guidelines | 形状、格子、境界条件、計算・収束、比較・感度・出典の報告 | WIND固有のファイル形式・乱流モデル欄は今回に強制しない | 4～11、付録 |
| R-04 | NASA-STD-7009B (2024) | 予定用途、モデル信頼性、プロジェクト側の承認 | NASA制度全体をプロジェクトに強制適用しない | 1、12、13 |
| R-05 | NASA-HDBK-7009B (2026) | 審査・利用可能性・制約の表現方法 | ハンドブックは規則そのものとは別の実務支援 | 1、12、13 |
| R-06 | Oberkampf & Trucano (2002) | Code／Solution Verification、Validation、誤差と不確かさの区別 | 一律の合格値を与えない | 3、6～9 |
| R-07 | Eça & Hoekstra (2014) | 格子細分化に基づく数値不確かさ評価の条件付き適用 | 適切な系統的格子列が必要 | 7 |
| R-08 | Celik et al. (2008) | 離散化誤差の報告、GCI等の条件と根拠 | 強制GCIでなく、成立性を検討 | 7 |
| R-09 | OpenVSP official V&V Studies at 3.50.5 tag | Hershey Bar／Ellipsoidの基本ベンチマーク | 新版・旧版の混同、形状設定と出力の相違に注意 | 6 |
| P-01 | V&V計画書v1.1 | 25 Test ID、対象QoI、仮の試験開始設定、数値・工学判定 | 計画書自体は空力的な精度を証明しない | 全章 |
| P-02 | DB詳細仕様v1.0.1とDDL | Run／Case／Point、元データ、比較、評価、証拠の追跡 | DBの実装済み範囲と仕様上の目標を区別 | 4、5、付録 |

## 2. 参考URL／書誌

- **R-01** ASME V&V 20: https://www.asme.org/codes-standards/find-codes-standards/standard-for-verification-and-validation-in-computational-fluid-dynamics-and-heat-transfer
- **R-02** AIAA G-077: https://webstore.ansi.org/standards/aiaa/AIAA0771998R2002
- **R-03** NASA NPARC: https://www.grc.nasa.gov/www/wind/valid/document.html
- **R-04** NASA-STD-7009B: https://standards.nasa.gov/standard/nasa/nasa-std-7009
- **R-05** NASA-HDBK-7009B: https://standards.nasa.gov/standard/NASA/NASA-HDBK-7009
- **R-06** Oberkampf & Trucano, *Progress in Aerospace Sciences* 38(3), 209–272 (2002): https://doi.org/10.1016/S0376-0421(02)00005-2
- **R-07** Eça & Hoekstra, *Journal of Computational Physics* 262, 104–130 (2014): https://doi.org/10.1016/j.jcp.2014.01.006
- **R-08** Celik et al., *Journal of Fluids Engineering* 130(7), 078001 (2008): https://doi.org/10.1115/1.2960953
- **R-09** OpenVSP 3.50.5 official V&V: https://github.com/OpenVSP/OpenVSP/blob/OpenVSP_3.50.5/examples/scripts/python_scripts/report.md
- **P-01** 同梱`reference/VSPAERO_3.50.5_VV_Plan_v1.1.md`
- **P-02** 同梱`reference/VSPAERO_VV_DB_Program_Spec_v1.0.1.md`、`schema.sql`、`test_catalog.csv`

## 3. 重要な用語と本文での使い分け

| 用語 | 本報告書での扱い |
|---|---|
| Code Verification（限定的） | 公開ベンチマークで特定機能の再現を確認。完全なソースコードの正しさの証明ではない |
| Solution Verification | メッシュ・後流等の数値依存性と、可能な範囲での数値不確かさの評価 |
| Validation | 風洞等の実験に対応する条件点・QoIで物理現実の再現性を評価 |
| Cross-code Comparison | パネル法／NSとの比較。独立した実験Validationとは別 |
| Numerical Screening | 試験開始値・暫定の差を用いた傾向確認。Engineering Acceptanceへ昇格させない |
| Engineering Acceptance | 設計・使用用途から承認された要求値に対する判定 |
| Applicability | 用途、QoI、機体形態、飛行条件、制限を明示した技術的な使用判断 |
| Independent Validation | 調整・チューニングに使わなかった、用途に対応する実験データでの評価 |
| External CD0 | VSPAEROの内蔵`CDo`とは区別して別Datasetに保存した積み上げ寄生抗力 |

## 4. 本プロジェクト固有の方法であり、規格の直接的な義務ではない項目

1. 主翼・尾翼をThin、胴体をThickとしたG3構成。
2. `WakeNumIter=5`、`RootWakeNodes=8`等の暫定開始値。
3. `CLi`を基本揚力、`CDiw`を基準の誘導抗力評価候補として扱う運用。
4. `CD_A=CDiw+CD0_build`および`CD_B=CDiw+(CDo-CDo_baseline)+CD0_build`の比較設計。
5. `APPLICABLE/CONDITIONAL/NOT_APPLICABLE/INSUFFICIENT_EVIDENCE`という内部の用途判定カテゴリ。
6. `|E_Q|+U_Q<=T_Q`を候補とする、正式に承認された場合の内部の保守的判定式。
7. 全25 Test IDの分割、SQLiteの17テーブル、Run・Artifact単位の証拠追跡。

**注意**：上記はプロジェクト固有の合理的な仮説・運用案である。実行後の感度・Validation結果に基づき承認されてはじめて標準設定や設計用途への正式適用となる。
