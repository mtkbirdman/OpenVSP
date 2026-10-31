# VSPAERO 3.50.5 試験別V&V評価記録テンプレート

> 1 Test IDにつき原則1記録を作る。複数ケース・評価量・参照データを含めてもよい。正式報告書の付録Aから参照する。空力計算の結果を入力するまでは`NOT_EVALUATED`を維持する。

| 管理項目 | 記録値 |
|---|---|
| Record ID | 《未設定：例 TEST-REC-GRID-01-001》 |
| Test ID | 《既存25 Test IDのいずれか》 |
| Test Title | 《test_catalog.csvで確認》 |
| Report ID／版 | 《未設定》 |
| 計画書版 | V&V Plan v1.1 |
| Study ID | 《未設定》 |
| 実施分類 | `REQUIRED`／`CONDITIONAL` |
| 試験実施状況 | `NOT_STARTED`（報告用の表示値。SQLite列の列挙値ではない） |
| 数値判定 | `NOT_EVALUATED` |
| 工学判定 | `NOT_EVALUATED` |
| 使用適格性 | `INSUFFICIENT_EVIDENCE` |
| 作成・レビュー日 | 《未設定》 |

## A. 試験目的・対象と制約

- 目的・前提：《計画書第7章とTest Catalogに照らして記入》
- 主要QoI：《単位・変数定義・出典を併記》
- 対象Geometry／飛行条件：《明示》
- 変更を許可する数値・形状パラメータ：《列挙》
- その他の固定条件：《明示》
- 適用不可・不足データ：《記入》
- 参照先と用途区分：`THEORY`／`PANEL`／`NS`／`WIND_TUNNEL`等。`CALIBRATION`／`VALIDATION`／`REFERENCE`の区別。

## B. CaseとRun一覧

| Case ID | Case Role | Geometry ID | Condition ID | Setting ID | Run ID／Point ID | Run status | Requested / Effective / Generated | Artifact ID |
|---|---|---|---|---|---|---|---|---|
| 《未設定》 | 《例 Medium》 | — | — | — | — | 未実行 | 《未確認》 | — |

差分として設定した項目以外が意図どおり一致しているかを必ず確認する。

## C. 数値結果

| Dataset ID | Quantity ID | QoI・成分 | Variant | 原値・単位 | 正規化値・単位 | 座標・基準点／処理方法 | 値状態 |
|---|---|---|---|---|---|---|---|
| — | — | — | — | — | — | — | `NOT_CALCULATED` |

### C.1 主要結果の可視化

- 図ID：《未設定》
- 軸・スケール・凡例：《記入》
- 作図プログラム・版：《記入》
- 図の元Quantity ID・条件・SHA-256：《記入》

## D. 数値Verification・感度評価

- 解析成否、実行の実効値と生成値の照合：《未評価》
- 反復履歴・残差・後流状態：《未評価》
- Coarse／Medium／Fineの実際のメッシュ情報・形状同等性：《未評価》
- 隣接水準差・収束傾向・非単調性：《未評価》
- GCI等を適用したか、適用しない場合の理由：《未評価》
- 数値不確かさ、方法、仮定、Method ID：《未算定》

| QoI | 方法 | 評価量／スクリーニング差 | 規則ID | Rule approval | 判定 | 理由 |
|---|---|---|---|---|---|---|
| — | — | — | — | 未承認 | `NOT_EVALUATED` | 《未評価》 |

## E. 参照値との比較・Validation（該当試験のみ）

| Comparison ID | Predicted quantity | Reference quantity | 比較種類 | 符号付き差 | 標準／拡張不確かさ | 条件の差 | Verdict |
|---|---|---|---|---|---|---|---|
| — | — | — | — | — | 《未算定》 | 《未確認》 | `NOT_EVALUATED` |

- Calibrationに使った資料と独立Validation資料の区分：《未設定》
- 参照データの品質・独立性・不確かさの出典：《未設定》
- データ変換・補間・補正とそれぞれの根拠：《未設定》
- 不一致原因の候補と、その裏付け：《未評価》
- 本試験が証明できない事項：《明記》

## F. Test ID別の特別な注意事項

| 試験群 | 追加記載事項 |
|---|---|
| VER | 解析／理論の定義、ベンチマークモデル版、参考解の近似誤差 |
| GEO | G1/G2/G3構成差、交線処理、モデル化差と数値差の識別 |
| GRID | Coarse/Medium/Fine、要素数、形状変化、細分化比、GCI適否 |
| WAKE | Wake反復・Wake節点・方式を別試験として実施 |
| XCD | 他ソルバーのモデル仮定、基準量、粘性寄与、定義差 |
| VAL | 風洞模型、装置・壁面・支持具・補正、測定不確かさ |
| STA | Steady FD、Adjoint、P/Q/R、無次元化・組合せ微係数 |
| CTRL | Control Group・Gain・偏向方向、摂動幅rad、非線形性 |
| DRAG | `CDiw`、`CDo`、外部CD0、基準CDo、重複算入防止、入力依存 |
| ACC | 用途、QoI、飛行条件領域、承認制限、承認記録 |

## G. 判定・考察・制限事項

**結果と考察は分けて記述する。**

1. 観測された数値事実：《未評価》
2. 根拠を伴う解釈：《未評価》
3. 未確認の原因仮説：《未評価》
4. 数値判定：`NOT_EVALUATED` ／Rule ID《未設定》
5. 工学判定：`NOT_EVALUATED` ／Rule ID《未設定》
6. 用途認定：`INSUFFICIENT_EVIDENCE`／Assessment ID《未設定》
7. 追加試験・アクションID：《未設定》
8. 上位報告書へ引用する結論文：《未設定》

## H. 検証・承認

| 項目 | 証拠 |
|---|---|
| 全Raw ArtifactのSHA-256照合 | 《未実施》 |
| DB Snapshotの一致 | 《未実施》 |
| ソルバー・実行環境の版照合 | 《未実施》 |
| 書込禁止の原本に改変がない | 《未実施》 |
| 作成者／レビュー者の確認 | 《未実施》 |
| 承認者と承認日時 | 《未承認》 |
