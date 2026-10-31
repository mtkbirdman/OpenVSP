# VSPAERO 3.50.5 Verification & Validation（V&V）評価報告書

> **正式な報告様式／未評価テンプレート v0.1**
>
> このファイルは実試験結果を記入して審査・発行するための報告書テンプレートである。**評価結果や使用認定を示す発行済み報告書ではない。** 未入力欄は証拠がそろうまで未入力・`NOT_EVALUATED`・`INSUFFICIENT_EVIDENCE`のままとし、架空の数値で埋めない。

| 文書管理項目 | 内容 |
|---|---|
| 文書番号 | `VSP-VV-REPORT-001`（提案番号。正式な管理番号は発行時に確認） |
| 文書名 | VSPAERO 3.50.5 V&V評価報告書 |
| テンプレート版 | `0.1` |
| 報告書の発行版 | **未発行** |
| テンプレート作成日 | 2026-10-23 |
| 報告対象Study ID | 《未設定：Study ID》 |
| 対象機体・形態 | 《未設定：機体型式、幾何版、G1/G2/G3の定義》 |
| 対象ソフトウェア | **OpenVSP／VSPAERO 3.50.5**（実行バイナリの同定は未確認） |
| 上位V&V計画 | `plan.md` |
| DB詳細仕様 | `database_spec.md`／同梱`schema.sql` |
| 対応試験 | 上位計画の全25 Test ID（変更不可） |
| 対象DBスナップショット | 《未作成：パス、SHA-256、抽出日時UTC》 |
| プログラムGit Commit／版 | 《未設定：commit／パーサ版／後処理版》 |
| 情報区分 | 《未設定：機密区分・配布範囲》 |
| 作成者・検証者・承認者 | 《未設定》 |
| 現在の技術的状態 | **空力V&Vは未完了。使用認定は行っていない。** |

## 記載規則と文書間の優先順位

1. **正本**：計画目的・試験定義・判定原則はV&V計画書v1.1、DB構造・実行履歴はDB詳細仕様v1.0.1とSQLite、発行時の技術判断は本評価報告書で管理する。矛盾があれば自動的に優先して採用せず、逸脱・改訂記録を作成する。
2. **入力欄**：`《未設定：...》`、`《未算定》`、`《該当するID》`等は編集対象。証拠がない場合に数字やPASSを作らない。
3. **実施状態／評価結果／用途判断は別物**。Runの`SUCCEEDED`は解析ファイルの取得成功を表し、数値Verification、実験Validation、設計への使用認定を意味しない。
4. 図表には必ず **図表ID、評価量の定義と単位、条件、出典（Run／Dataset／Quantity／Artifact ID）、生成手法／版**を付ける。
5. 本文には設計上の判断に必要な要点を残し、計算点ごとの全行・コンソールログ・生データはSQLiteと付録で追跡する。
6. 公表文献は採用した考え方の根拠であり、規格への全面準拠や第三者認証を主張しない。方法の変更・逸脱を記録する。
7. `CALIBRATION`に使った資料を独立`VALIDATION`として流用しない。既存パネル法・NS解析はCross-code Comparisonであり、原則として風洞Validationの代替ではない。
8. 未実施・未評価・比較不能な試験を、空欄のまま黙って省略しない。状態、理由、影響、追加作業を明記する。

---

## 1. 総括（Executive Summary）

### 1.1 審査対象と設計用途

| 項目 | 要約 |
|---|---|
| Intended use（予定用途） | 《未設定：例として全機定常6分力・初期安定性評価等を列記。実際の予定用途のみ入力》 |
| 機体・構成 | 《未設定：形状版、翼・尾翼Thin／胴体Thickの具体的対応》 |
| 飛行条件の評価範囲 | 《未設定：Mach、Re、α、β、舵角、トリム・質量特性》 |
| 評価対象QoI | 《未設定：評価量・定義・用途別重要度》 |
| 評価に使用した主たる根拠 | 《未設定：理論、風洞、他ソルバー、分割・後流試験》 |
| 結果の要約 | **NOT_EVALUATED：試験・比較結果は未記入** |
| 使用の総合判断 | **INSUFFICIENT_EVIDENCE：現在は認定証拠がない** |
| 主たる制限事項 | 《未設定。資料不足や適用限界があれば列記》 |

### 1.2 主要結論・決定事項

| Decision ID | 設計上の判断・結論 | 対象QoI／条件 | 判定 | 根拠Assessment ID | 承認者／日付 |
|---|---|---|---|---|---|
| DEC-01 | 《未評価：採用する基準メッシュ》 | 《未設定》 | NOT_EVALUATED | — | — |
| DEC-02 | 《未評価：採用するWake設定》 | 《未設定》 | NOT_EVALUATED | — | — |
| DEC-03 | 《未評価：定常空力の使用可能範囲》 | 《未設定》 | INSUFFICIENT_EVIDENCE | — | — |
| DEC-04 | 《未評価：安定微係数／舵効き》 | 《未設定》 | INSUFFICIENT_EVIDENCE | — | — |
| DEC-05 | 《未評価：CD0補正方法A/B》 | 《未設定》 | INSUFFICIENT_EVIDENCE | — | — |

### 1.3 結論の強さと未解決事項

《記入：評価済み事実、条件付き結論、未確認の仮説を区別。未解決事項は13章のAction IDを引用する。》

## 2. 評価目的・対象・適用範囲

### 2.1 使用シナリオと評価要求

| Use Case ID | 使用目的・意思決定 | 必要QoI | 条件領域 | 設計要求の出典 | 未達時の影響 |
|---|---|---|---|---|---|
| UC-01 | 《未設定》 | 《未設定》 | 《未設定》 | 《文書ID／頁／要求番号》 | 《未設定》 |

要求精度は**設計上の下流影響（性能、安定性、荷重、トリム、飛行力学、制御）から導く**ことを基本とする。他ソルバーとの一致率や風洞データの散らばりを、そのまま設計上の許容誤差へ転用しない。

### 2.2 対象形状と評価段階

| 段階 | 構成 | 解析表現 | 本報告書での役割 | 形状ID・版 |
|---|---|---|---|---|
| G1 | 主翼単体 | Thin | 基本翼特性 | 《未設定》 |
| G2 | 主翼＋胴体 | 翼Thin／胴体Thick | 翼胴干渉 | 《未設定》 |
| G3 | 主翼＋胴体＋水平・垂直尾翼 | 主翼・尾翼Thin／胴体Thick | 全機モデル（基準形態） | 《未設定》 |
| M1/M2 | 条件付きThick比較 | 《未設定》 | モデル化感度 | 《未設定》 |

解析用形状、風洞模型、既存パネル法、NS形状の差を表`TAB-02-02`として記録する。翼根・尾翼取付、胴体交線、端部Cap、支持具、推進装置、舵面の物理的なモデル化差を確認する。

### 2.3 包絡・除外範囲

| 評価軸 | 計画上の範囲 | 参照証拠の範囲 | 検証後の使用範囲 | 除外理由・残課題 |
|---|---|---|---|---|
| Mach | 《未設定》 | 《未設定》 | 未認定 | 《未設定》 |
| Re（基準長明記） | 《未設定》 | 《未設定》 | 未認定 | 《未設定》 |
| α | 《未設定》 | 《未設定》 | 未認定 | 《未設定》 |
| β | 《未設定》 | 《未設定》 | 未認定 | 《未設定》 |
| 舵角／トリム | 《未設定》 | 《未設定》 | 未認定 | 《未設定》 |
| 動的周波数・振幅 | 《未設定》 | 《未設定》 | 未認定 | 《未設定》 |

強い剥離・失速、遷音速衝撃波と境界層の干渉、強い非定常推進干渉は、固有の物理モデルと独立Validationなしに認定しない。境界条件を超える内挿・外挿の可否は個別審査とする。

## 3. 評価方法・指標・判定基準

### 3.1 V&Vの分解

| 区分 | 評価内容 | 主なTest ID | 主張できる範囲 |
|---|---|---|---|
| Limited Code Verification | 公式ベンチマーク・理論結果の再現 | VER-01/02 | 選定機能の再現確認。全実装の正当性を保証しない |
| Solution Verification | 分割、後流、反復、幾何表現の感度 | GRID・WAKE・GEO | 指定条件・QoIの数値依存性の評価 |
| Cross-code Comparison | パネル法、NSと比較 | XCD-01/02 | 物理モデル・設定差の診断。実験Validationとは区別 |
| Validation | 独立した風洞データと比較 | VAL-01/02 | 対応条件点・評価量での現実再現性の評価 |
| Special Methods | 微係数・舵効き・抗力構築 | STA・CTRL・DRAG | 定義一致、微係数、補正と参照整合性 |
| Applicability | 指定した設計用途に対する使用判断 | ACC-01 | 用途別・条件別の承認と制約 |

### 3.2 評価量と定義・符号・正規化

| QoI | ネイティブ名称 | 正規化・座標系 | 参照面積・長さ・基準点 | 比較基準 |
|---|---|---|---|---|
| 揚力・誘導抗力 | `CLi`, `CDi`, `CLiw`, `CDiw` | 《未設定：出力定義・風軸／機体軸》 | 《未設定》 | 《該当資料》 |
| 3軸モーメント | `CMix`, `CMiy`, `CMiz` 等 | 《軸・符号変換を確認》 | 《未設定》 | 《該当資料》 |
| 静的微係数 | `CL_alpha` 等 | rad⁻¹／その他を明記 | 《未設定》 | 《該当資料》 |
| 動安定微係数 | `Cm_q`, `Cl_p`, `Cn_r` 等 | 無次元角速度定義・定常／非定常区別 | 《未設定》 | 《該当資料》 |
| 舵効き | `Cm_delta_e` 等 | rad⁻¹、Gain・正方向明記 | 《未設定》 | 《該当資料》 |
| 抗力ポーラー | `CDiw`, `CDo`, 外部CD0, `CD_A/B` | 内蔵抗力・外部抗力の包含範囲を識別 | 《未設定》 | 《該当資料》 |

`CLi/CDi`と`CLiw/CDiw`は**別の評価系列**として保管する。`w`はwakeであり造波抗力ではない。両系列を加算しない。`CLtot/CDtot`等に含まれる`CDo`を、外部CD0と重複加算しない。`CDiw`を常に厳密な誘導抗力として扱えるとは限らないため、使用条件と参照比較を記載する。

### 3.3 誤差・不確かさの表示

同じQoI・単位・座標系・基準値へ整合させた値について、符号付き比較差を次で定義する。

$$
E_Q=Q_{\mathrm{VSPAERO}}-Q_{\mathrm{reference}}
$$

相対差`100 × E_Q / Q_reference`は**ゼロ近傍や符号反転のある参照値に対して強制的に報告しない**。必要に応じて、固定の代表スケール、データ系列の振幅、または絶対差で評価する。固定スケールは事前に承認する。

不確かさは利用可能な成分・定義・包含水準を分けて記載する。例として、互いに独立な標準不確かさとして正当に扱える場合に限り、以下を利用できる。

$$
u_{\mathrm{val},Q}^{2}=u_{\mathrm{num},Q}^{2}+u_{\mathrm{input},Q}^{2}+u_{\mathrm{exp},Q}^{2}
$$

相関がある場合は共分散項を加える。モデル形式に起因する偏りや、未検証領域への外挿不確かさまで自動的にこの式で説明できるわけではない。標準不確かさ・拡張不確かさ・信頼水準・包含係数`k`を混在させない。

**採用した評価式・集計方法・数値不確かさ推定方法は、Method ID・版・根拠文献・適用条件とともに記録する。** GCIは系統的な細分化系列等の成立条件を確認した場合だけ評価し、不成立なら感度結果として報告する。

### 3.4 数値・工学・使用認定の判定規則

| 判定対象 | 事前に必要な条件 | 結果 |
|---|---|---|
| Run実行健全性 | 実行完了、出力・設定・版の一致、原本の検証 | `run.status`（V&V合格とは独立） |
| Numerical Acceptance | 反復・分割・後流等、承認された数値的管理基準 | PASS／FAIL／INCONCLUSIVE／NOT_EVALUATED |
| Engineering Acceptance | 指定用途のQoI、承認済み要求精度、妥当な参照比較 | PASS／FAIL／INCONCLUSIVE／NOT_EVALUATED |
| Applicability | 上記評価、証拠の独立性・代表性、レビュー、使用制限 | APPLICABLE／CONDITIONAL／NOT_APPLICABLE／INSUFFICIENT_EVIDENCE |

**仮の数値スクリーニング基準を正式なEngineering PASSへ流用しない。** 例えば`|E_Q|+U_Q ≤ T_Q`を採用する場合も、`E_Q`の性質、`U_Q`の意味・包含水準、設計要求`T_Q`、適用する意思決定を明示し、プロジェクトの承認済み基準として管理する。これはASME V&V 20の一律のPASS式ではない。

| Rule ID | QoI／対象領域 | Numerics／Engineering／Applicability | 閾値・定義 | 承認状態 | 出典・承認者 |
|---|---|---|---|---|---|
| 《未設定》 | 《未設定》 | 《未設定》 | 《未設定》 | `PROVISIONAL`／未設定 | 《未設定》 |

### 3.5 Validation参照資料の資格

`dataset.usage_class`の`CALIBRATION`／`VALIDATION`／`REFERENCE`と、`evidence_class`（A/B/C）を記録する。風洞条件と解析条件の差、実験測定不確かさ、試験装置・支持具・壁面補正、補間の有無を明記する。参照資料が不十分でも参考比較は可能だが、不確かさに基づく正式認定へ自動昇格させない。

## 4. 解析対象・環境・モデル設定

### 4.1 実行ソフトウェアと再現性

| 項目 | 記録値 | 証拠ファイル／ハッシュ |
|---|---|---|
| OpenVSP GUI/API版 | 《3.50.5との一致を実測》 | 《未設定》 |
| VSPAERO実行バイナリ版 | 《3.50.5との一致を実測》 | 《未設定》 |
| APIファイル・実行バイナリSHA-256 | 《未設定》 | 《未設定》 |
| Python・OS・CPU・スレッド数 | 《未設定》 | 《未設定》 |
| 解析ユーティリティGit Commit | 《未設定》 | 《未設定》 |
| VVDBソフトウェア版 | 《未設定：現在の試験済範囲を併記》 | 《未設定》 |
| 計画書・設定ファイルのSHA-256 | 《未設定》 | 《未設定》 |
| データベースSnapshotのSHA-256 | 《未設定》 | 《未設定》 |

### 4.2 Geometryと最終メッシュ

| Geometry ID | モデル段階 | 形状ファイル | Thin Set | Thick Set | 機体／模型との同等性 | メッシュ証拠 |
|---|---|---|---|---|---|---|
| 《未設定》 | G3 | 《*.vsp3》 | 《主翼・尾翼》 | 《胴体》 | UNVERIFIED | 《*.vspgeom／図》 |

モデルの埋め込み、胴体・翼・尾翼交線、後縁、Cap、舵面のGain・物理偏向方向、左右対称性の違いを確認した結果を記録する。設定側の`Tess_W/SectTess_U`だけでは最終的なVSPAERO交線後の要素状態を保証しない。

### 4.3 代表数値設定

| 項目 | 計画書の試験開始候補 | 実際に使用した値 | 合理的な選定根拠／Test ID |
|---|---|---|---|
| 主翼`Tess_W` | 31 | 未確認 | GRID-01/02 |
| 主翼`SectTess_U` | 12 | 未確認 | GRID-01/02 |
| 胴体`Tess_W` | 25 | 未確認 | GRID-01/03 |
| `WakeNumIter` | 5（感度：5/10/20） | 未確認 | WAKE-01 |
| `RootWakeNodes` | 8（感度：8/12/16） | 未確認 | WAKE-02 |
| `FixedWakeFlag` | false | 未確認 | WAKE-03 |
| `ImplicitWake` | false | 未確認 | WAKE-03 |
| `Symmetry`／`RefFlag` | OFF／`MANUAL_REF` | 未確認 | GEO・全比較 |
| `Sref`,`bref`,`cref`、基準点 | 未確定（実行前必須） | 未確認 | O-03 |

**計画書の開始値をそのまま「検証済み標準値」と表現しない。** 実行ごとのRequested／Effective／Generatedを添付記録へ保存し、差があれば比較可能性を審査する。

### 4.4 基準飛行条件・実行条件

| Condition ID | Mach | Re（基準長） | α[deg] | β[deg] | 舵角・トリム | 参照量・座標系 | 用途 |
|---|---|---|---|---|---|---|---|
| 《未設定》 | — | — | — | — | — | — | 《未設定》 |

計画書O-01～O-09の未確定事項と解決状況は付録Dの管理表へ記録する。

## 5. 試験実施状況・逸脱管理

### 5.1 全25 Test IDの実施・評価状況

> この表のTest ID・名称・必須／条件付きの区分は、`test_catalog.csv`を正本として転記する。`run.status`／各テストの判定／最終使用認定を混同しない。

| Test ID | 試験名 | 要求区分 | 実施状況 | 評価状態 | 主なRun／参照Dataset | 備考・逸脱 |
|---|---|---|---|---|---|---|
| `VER-01` | Hershey Bar理論照合 | REQUIRED | 未実施 | NOT_EVALUATED | — | — |
| `VER-02` | Ellipsoid理論Cp | REQUIRED | 未実施 | NOT_EVALUATED | — | — |
| `GEO-01` | 主翼単体G1 | REQUIRED | 未実施 | NOT_EVALUATED | — | — |
| `GEO-02` | 主翼＋胴体G2 | REQUIRED | 未実施 | NOT_EVALUATED | — | — |
| `GEO-03` | 全機G3 | REQUIRED | 未実施 | NOT_EVALUATED | — | — |
| `GEO-04` | Thin/Thick形態比較 | CONDITIONAL | 未実施 | NOT_EVALUATED | — | — |
| `GRID-01` | 全機代表点の格子感度 | REQUIRED | 未実施 | NOT_EVALUATED | — | — |
| `GRID-02` | 第二飛行点の格子感度 | REQUIRED | 未実施 | NOT_EVALUATED | — | — |
| `GRID-03` | 局所分割・NRef等 | CONDITIONAL | 未実施 | NOT_EVALUATED | — | — |
| `WAKE-01` | 後流反復5/10/20 | REQUIRED | 未実施 | NOT_EVALUATED | — | — |
| `WAKE-02` | Wake節点8/12/16 | REQUIRED | 未実施 | NOT_EVALUATED | — | — |
| `WAKE-03` | Wake方式の比較 | CONDITIONAL | 未実施 | NOT_EVALUATED | — | — |
| `XCD-01` | 既存パネル法との照合 | REQUIRED | 未実施 | NOT_EVALUATED | — | — |
| `XCD-02` | NS解析との照合 | REQUIRED | 未実施 | NOT_EVALUATED | — | — |
| `VAL-01` | 風洞縦空力 | REQUIRED | 未実施 | NOT_EVALUATED | — | — |
| `VAL-02` | 風洞横方向 | REQUIRED | 未実施 | NOT_EVALUATED | — | — |
| `STA-01` | Steady/Adjoint比較 | REQUIRED | 未実施 | NOT_EVALUATED | — | — |
| `STA-02` | 微係数の分割/Wake感度 | REQUIRED | 未実施 | NOT_EVALUATED | — | — |
| `STA-03` | P/Q/R非定常 | CONDITIONAL | 未実施 | NOT_EVALUATED | — | — |
| `CTRL-01` | エレベータ舵効き | REQUIRED | 未実施 | NOT_EVALUATED | — | — |
| `CTRL-02` | エルロン舵効き | REQUIRED | 未実施 | NOT_EVALUATED | — | — |
| `CTRL-03` | ラダー舵効き | REQUIRED | 未実施 | NOT_EVALUATED | — | — |
| `DRAG-01` | 外部CD0方法A | REQUIRED | 未実施 | NOT_EVALUATED | — | — |
| `DRAG-02` | 外部CD0方法B | REQUIRED | 未実施 | NOT_EVALUATED | — | — |
| `ACC-01` | 使用認定 | REQUIRED | 未実施 | NOT_EVALUATED | — | — |

#### ステータス集計（入力時に算定）

| 集計対象 | 値 |
|---|---|
| Test ID総数 | 25（固定） |
| `REQUIRED` | 《入力：Catalogから集計》 |
| `CONDITIONAL` | 《入力：Catalogから集計》 |
| 実施済み | 《未集計》 |
| 未実施／実施失敗／対象外 | 《未集計》 |
| PASS／FAIL／INCONCLUSIVE／NOT_EVALUATED | 《未集計》 |

未実施の条件付き試験は`NOT_APPLICABLE`または未実施理由を別列へ記録する。`NOT_APPLICABLE`を判定`PASS`に換算しない。

### 5.2 試験計画からの逸脱

| Deviation ID | Test ID | 計画内容 | 実際の内容 | 理由 | QoI・認定への影響 | 承認記録 |
|---|---|---|---|---|---|---|
| DEV-001 | 《未設定》 | 《未設定》 | 《未設定》 | 《未設定》 | 《未評価》 | 《未設定》 |

### 5.3 解析ファイル・ソフトウェアの健全性

《記入：3.50.5での実行確認、状態遷移・原本ハッシュ照合、missing／changed artifact、解析エラーの除外基準。ソフトウェア単体テスト（SW）は空力V&Vとは別節で明示。》

## 6. 基本Verification（VER）

### 6.1 VER-01：Hershey Bar

- **目的**：有限翼の揚力傾斜・誘導抗力・分布荷重等について、明示された理論・公開参照値との対応を確認する。
- **解析形状・理論の版**：《未設定：公式モデル、使用した式・定数、適用領域》
- **条件・Run ID**：《未設定》
- **結果表 TAB-06-01**：`QoI | Theory/Reference | VSPAERO | Difference | 数値依存性 | Rule ID | Verdict | 参照ID`。
- **図 FIG-06-01**：《挿入予定：揚力特性・翼幅分布》
- **考察**：《未評価：理論の近似誤差、VLM/Panel相違、分割・後流感度を区別》
- **判定**：`NOT_EVALUATED`

### 6.2 VER-02：Ellipsoid

- **目的**：Thick体の面圧分布について理論解との一致性を確認する。
- **座標点とCpの定義**：《未設定：解析点の取出し、各軸・面上の比較点、補間方式》
- **条件・Run ID**：《未設定》
- **結果表 TAB-06-02**：`Axis/Surface | Point count | RMS difference | Maximum difference | Numerical check | Verdict`。
- **図 FIG-06-02**：《挿入予定：軸別Cp分布、差分》
- **判定**：`NOT_EVALUATED`

公式ベンチマークが成功しても、機体G3全体のValidation完了とは解釈しない。

## 7. 幾何モデル・メッシュ・後流のSolution Verification（GEO／GRID／WAKE）

### 7.1 GEO-01～GEO-04：構成要素と解析表現

| Test ID | 比較対象 | 主要QoI | 差の解釈・留意点 | 結果・証拠 |
|---|---|---|---|---|
| GEO-01 | 主翼G1単体 | `CLi`, `CMiy`, `CDiw`、荷重 | 基本モデルの健全性。単独でValidationではない | 《未評価》 |
| GEO-02 | G2−G1 | 同上、寄与差 | 翼胴干渉と形状・参照値差を分離 | 《未評価》 |
| GEO-03 | G3−G2 | 同上、尾翼寄与 | 翼・尾翼干渉、ダウンウォッシュの影響 | 《未評価》 |
| GEO-04 | Thin／Thick比較 | 同上、Cp | 解析表現そのものの変更。単純なメッシュ収束とは別 | 《未評価》 |

**FIG-07-01**：《G1/G2/G3の形状比較図》、**TAB-07-01**：《同一基準量での係数差》。

### 7.2 GRID-01～GRID-03：表面分割収束

| Case ID | Geometry ID | Condition ID | Setting ID | 実メッシュ数／品質 | CLi | CMiy | CDiw | CLiw | CDi | Run ID |
|---|---|---|---|---|---|---|---|---|---|---|
| 《Coarse》 | — | — | — | — | — | — | — | — | — | — |
| 《Medium》 | — | — | — | — | — | — | — | — | — | — |
| 《Fine》 | — | — | — | — | — | — | — | — | — | — |

- **TAB-07-02**：《隣接水準差、代表格子寸法、解析時間、残差等》
- **FIG-07-02**：《QoI別格子感度、必要なら特定接合部を拡大》
- **評価方法**：《各QoIの中・細差、単調／振動的傾向、形状・トポロジー変更の有無を説明》
- **GCI／外挿**：《適用可否、採用／不採用理由。成立条件を確認できた場合のみ数値を示す》
- **結果**：`NOT_EVALUATED`

GRID-03は条件付きの局所分割・`NRef`等の一因子感度であり、全機格子系列とは区別する。

### 7.3 WAKE-01～WAKE-03：後流の反復・節点・方式

| Test ID | 変更する唯一の主要因 | 基準候補 | 必要な比較 | 判定 |
|---|---|---|---|---|
| WAKE-01 | `WakeNumIter` | 5/10/20 | QoI差、反復／後流健全性 | NOT_EVALUATED |
| WAKE-02 | `RootWakeNodes`／`NumWakeNodes` | 8/12/16 | QoI差、後流離散化の影響 | NOT_EVALUATED |
| WAKE-03 | Fixed／Implicitなど方式 | 必要時 | 方式によるモデル差・収束差 | NOT_EVALUATED |

- **TAB-07-03**：《各RunのRequested／Effective／Generated、残差、`CLi`・`CLiw`・`CDi`・`CDiw`》
- **FIG-07-03**：《Wake反復または節点数とQoI》
- `CLi`と`CLiw`、`CDi`と`CDiw`の差がゼロになることを、単独の合格条件として扱わない。
- **基準設定への採否**：《未決。選定根拠をDEC-02へ参照》

### 7.4 数値Verificationの結論

| QoI | 対象条件 | メッシュ感度 | Wake感度 | 残差・反復 | Numerical Verdict | Assessment ID |
|---|---|---|---|---|---|---|
| 《未設定》 | 《未設定》 | 《未算定》 | 《未算定》 | 《未算定》 | NOT_EVALUATED | — |

## 8. 既存パネル法・NS解析との比較（XCD）

### 8.1 XCD-01：既存パネル法

《記入：モデル化（翼面、胴体、尾翼）、後流評価、座標・基準点・抗力定義、計算条件・分割・収束状態を一致表で示す。》

### 8.2 XCD-02：Navier–Stokes解析

《記入：NS解析形状、境界条件、粘性・乱流モデル、遷移、格子依存性、収束、支柱・風洞補正、抗力成分を明示する。NS解析を無条件の真値としない。》

| Comparison ID | Reference Dataset ID | QoI | FlightCondition | Comparison kind | VSPAERO | Reference | Difference | 計算／定義差の解釈 |
|---|---|---|---|---|---|---|---|---|
| 《未設定》 | — | — | — | — | — | — | — | 《未評価》 |

**FIG-08-01**：《縦空力・横方向・局所分布のコード間比較》。

直接比較ができない場合は`DIRECT / CONVERTED / INTERPOLATED / ILLUSTRATIVE`を区別し、使用した変換式と不確かさを記録する。

## 9. 風洞Validation（VAL）

### 9.1 試験装置・計測方法・風洞模型

| 項目 | 記録内容 | 出典・証拠 |
|---|---|---|
| 風洞名／実施日／試験報告書版 | 《未設定》 | 《未設定》 |
| 模型幾何・舵角・基準点 | 《未設定》 | 《未設定》 |
| 迎角・横滑り角測定／支持具／壁面補正 | 《未設定》 | 《未設定》 |
| 動圧・Mach・Re・大気条件 | 《未設定》 | 《未設定》 |
| 力・モーメント測定法／較正履歴 | 《未設定》 | 《未設定》 |
| 不確かさ成分／相関／包含水準 | 《未設定》 | 《未設定》 |
| `CALIBRATION`と独立`VALIDATION`の区分 | 《未設定》 | 《未設定》 |

### 9.2 VAL-01：縦方向空力

| Comparison ID | Condition ID | α | QoI | 計算値 | 風洞値 | 差 `E_Q` | 不確かさ（定義付き） | Rule ID | Verdict |
|---|---|---|---|---|---|---|---|---|---|
| 《未設定》 | — | — | CL／Cm／CD | — | — | — | 《未算定》 | — | NOT_EVALUATED |

**FIG-09-01**：《CL–α、Cm–α、CD–CL、必要に応じて差分・不確かさ》。

### 9.3 VAL-02：横・方向空力

| Comparison ID | Condition ID | β | QoI | 計算値 | 風洞値 | 差 | 不確かさ | Rule ID | Verdict |
|---|---|---|---|---|---|---|---|---|---|
| 《未設定》 | — | — | CY／Cl／Cn | — | — | — | 《未算定》 | — | NOT_EVALUATED |

**FIG-09-02**：《CY、Cl、Cnとβの関係》。静的風洞試験のデータを`Cm_q`等の動的Validationとして扱わない。

### 9.4 不確かさと比較条件の妥当性

《記入：直接比較可能性、換算・補間、測定不確かさと数値不確かさの推定根拠、モデル形式差を区別。量の定義や誤差構造が不明なら、適用判定の不足証拠とする。》

### 9.5 Validation結論

| QoI | Validation点 | 観測差の範囲 | 利用可能な不確かさ | 設計要求 | 判定 | 限界 |
|---|---|---|---|---|---|---|
| 《未設定》 | 《未設定》 | 《未算定》 | 《未算定》 | 《未承認》 | INCONCLUSIVE | 《未設定》 |

## 10. 安定微係数・舵効き（STA／CTRL）

### 10.1 安定微係数の定義

角速度は、比較対象と一致することを条件に、たとえば

$$
\hat p=\frac{pb}{2V},\qquad \hat q=\frac{q\bar c}{2V},\qquad \hat r=\frac{rb}{2V}
$$

のように無次元化する。原データが角速度のrad/s当たり、無次元角速度当たり、または強制振動から得た結合係数なのかを**必ず確認**する。`C_{m_q}`と`C_{m_{q+\dot\alpha}}`等を無条件に同一視しない。

### 10.2 STA-01／02：Steady、Adjoint、感度

| Test ID | 方法 | QoI | 値・単位 | Mesh／Wake | 参照／比較 | Verdict | 備考 |
|---|---|---|---|---|---|---|---|
| STA-01 | Steady FD vs Adjoint | `Cm_q`, `Cl_p`, `Cn_r` 等 | 《未算定》 | 《未設定》 | 《定義一致後》 | NOT_EVALUATED | 《未設定》 |
| STA-02 | 格子・Wake感度 | 同上 | 《未算定》 | 《未設定》 | 《未設定》 | NOT_EVALUATED | 《未設定》 |

### 10.3 STA-03：P/Q/R非定常（条件付き）

| 項目 | 記録値 |
|---|---|
| 物理・数値的方法 | 《未設定》 |
| Forced oscillation／P/Q/Rのモード | 《未設定》 |
| Reduced frequencyの定義・値 | 《未設定》 |
| 摂動振幅・角度の単位 | 《未設定》 |
| 時間刻み・周期数・後流設定 | 《未設定》 |
| 位相・組合せ微係数 | 《未設定》 |
| 動的参照資料の有無・比較可能性 | 《未設定》 |
| 判定・理由 | NOT_EVALUATED／《未設定》 |

非定常データがなくても定常微係数の感度調査はできるが、動的Validation済みとは主張しない。

### 10.4 CTRL-01～03：舵角特性と中央差分

| Test ID | Control Group／Gainの検証 | δ点・使用Run | QoI | rad当たり微係数 | 非線形性 | 参照比較・判定 |
|---|---|---|---|---|---|---|
| CTRL-01 | 《Elevatorの物理偏向確認》 | 《未設定》 | `Cm_delta_e`, `CL_delta_e` | 《未算定》 | 《未評価》 | NOT_EVALUATED |
| CTRL-02 | 《Aileronの左右・Gain確認》 | 《未設定》 | `Cl_delta_a`, `Cn_delta_a` | 《未算定》 | 《未評価》 | NOT_EVALUATED |
| CTRL-03 | 《Rudderの偏向方向確認》 | 《未設定》 | `Cn_delta_r`, `CY_delta_r` | 《未算定》 | 《未評価》 | NOT_EVALUATED |

中央差分は舵角幅を記録し、radを分母とする。例：

$$
C_{m_{\delta_e}}\approx\frac{C_m(+\Delta\delta_e)-C_m(-\Delta\delta_e)}{2\Delta\delta_e},
\qquad \Delta\delta_e\;[\mathrm{rad}]
$$

**FIG-10-01～03**：《操舵角と空力係数、局所傾き、非線形性・正負対称性》。大舵角での剥離・隙間・実際のヒンジ形状の影響は、解析手法の限界として別記する。

## 11. 外部CD0による抗力ポーラー（DRAG）

### 11.1 使用する抗力成分の宣言

| Quantity / Dataset ID | 成分 | 内包される物理 | Mach・Re・形状の対応 | 原本 |
|---|---|---|---|---|
| 《未設定》 | `CDiw` | 後流側評価の抗力 | 《未検証》 | `.polar` |
| 《未設定》 | `CDo` | VSPAERO内蔵寄生抗力項 | 《未検証》 | `.polar` |
| 《未設定》 | 外部`CD0_build` | Parasite Drag Toolの積み上げ成分 | 《未検証》 | 《未設定》 |

### 11.2 方法A（基準）

$$
C_{D,A}(\alpha)=CDiw(\alpha)+C_{D0,\mathrm{build}}
$$

VSPAEROの`CDtot`や`CDwtot`を元に外部CD0を追加しない。CD0の基準翼面積と、寄生抗力に含まれる項目を記録する。`CDiw`の誘導抗力としての妥当性は対応する検証範囲に限定する。

### 11.3 方法B（比較）

$$
C_{D,B}(\alpha)=CDiw(\alpha)
 +\left[CDo(\alpha)-CDo_{\mathrm{baseline}}\right]
 + C_{D0,\mathrm{build}}
$$

`CDo_baseline`の基準条件と、その増分を加える物理的理由を明記する。方法Bは必ずしも方法Aより正確とは限らない。`DRAG-01`と`DRAG-02`を別の評価とし、派生値は元データを上書きせず`DERIVED` Datasetへ記録する。

| FlightCondition ID | CL | CDiw | 外部CD0 | CDo | Baseline CDo | CD_A | CD_B | 風洞等の参照値 | Verdict |
|---|---|---|---|---|---|---|---|---|---|
| 《未設定》 | — | — | — | — | — | — | — | — | NOT_EVALUATED |

**FIG-11-01**：《CD–CL抗力ポーラー、A/B、風洞・NS参照値、許容精度》。

## 12. 最終的な使用適格性評価（ACC-01）

### 12.1 用途・QoI・領域ごとの判断

> **本章が設計上の最終成果である。** 判断はVSPAERO 3.50.5全体への一括認証ではなく、明示された用途・形態・QoI・飛行条件に限る。

| Use Case ID | QoI | Geometry／飛行条件・舵角 | Numerical Status | Engineering Status | 独立Validation証拠 | 使用判断 | 使用条件・制限 | Assessment ID |
|---|---|---|---|---|---|---|---|---|
| UC-01 | 《未設定》 | 《未設定》 | NOT_EVALUATED | NOT_EVALUATED | 未確認 | INSUFFICIENT_EVIDENCE | 使用未認定 | — |

許容される判定値は `APPLICABLE`／`CONDITIONAL`／`NOT_APPLICABLE`／`INSUFFICIENT_EVIDENCE`。`CONDITIONAL`は審査によって制限・追加確認・責任者を明示してはじめて採用する。未実施や曖昧な試験結果を条件付き適用へ自動移行しない。

### 12.2 使用許可文（発行時の記入欄）

| 項目 | 内容 |
|---|---|
| 許可対象 | 《未承認：用途、QoI、機体構成》 |
| 許可する飛行条件 | 《未承認：Mach、Re、α、β、舵角、トリム》 |
| 必須の形状・メッシュ・後流設定 | 《未承認：設定ID、改訂版》 |
| 使用する派生係数・CD0方法 | 《未承認》 |
| 使用上の禁止事項 | 《未承認：除外領域・物理現象・外挿》 |
| 再V&Vを要する変更 | 《未承認：Geometry、Solver、後処理、飛行条件変更》 |
| 審査・承認情報 | 《未承認》 |

### 12.3 予測領域への外挿と追加試験

《記入：Validation点の間で使う領域、検証点から離れた領域、物理現象が変わる領域を識別し、外挿の判断根拠・追加試験の優先度を記す。》

## 13. 結論・残課題・文書承認

### 13.1 結論

《実試験後に記入。検証済みの事項、認定できない事項、将来の適用へ残るリスクを短く記述。試験未実施であれば「V&V結果未取得・使用認定不可」と明記。》

### 13.2 要対応事項（Action Log）

| Action ID | 課題・不足証拠 | 関連Test ID／QoI | 重要度 | 担当 | 期限 | 完了証拠 | 状態 |
|---|---|---|---|---|---|---|---|
| ACT-001 | 《機体固有の飛行条件・承認済み要求精度の確認》 | ACC-01 | 高 | 《未設定》 | 《未設定》 | — | OPEN |
| ACT-002 | 《VSPAERO 3.50.5実機結合試験》 | VER／全Run | 高 | 《未設定》 | 《未設定》 | — | OPEN |
| ACT-003 | 《風洞の条件・不確かさ・独立Validation区分の確認》 | VAL-01/02 | 高 | 《未設定》 | 《未設定》 | — | OPEN |

### 13.3 承認欄

| 役割 | 氏名／所属 | 審査対象 | 判定・意見 | 日付 |
|---|---|---|---|---|
| 作成者 | 《未設定》 | 計算・出典・転記 | 未実施 | — |
| 解析検証者 | 《未設定》 | 再現性、式、整合性 | 未実施 | — |
| 独立レビュー者 | 《未設定》 | 根拠、欠落、判断 | 未実施 | — |
| 技術承認者 | 《未設定》 | 用途ごとの使用認定 | **未承認** | — |

### 13.4 改訂・発行履歴

| 版 | 日付 | 変更概要 | 作成者／承認者 | 対応DB Snapshot |
|---|---|---|---|---|
| Template v0.1 | 2026-10-23 | 既存計画書v1.1およびDB仕様v1.0.1に対応した未評価テンプレート | 未承認 | なし |

---

# 付録

## 付録A：試験別詳細記録

`test_record_template.md`を試験ごとに複製し、`Test ID`、目的、予定Case、実施Run、原本、変数定義、比較、不確かさ、判定、考察を記録する。25件の各Test IDの独立記録をDBから参照できるようにする。必要に応じて1試験内に複数のQoI・比較条件点を設ける。

## 付録B：図表・原本一覧（Artifact Register）

| 図表ID | 内容 | 生成方法・版 | 対応Test／Case／Run | Dataset／Quantity ID | 原本相対パス | SHA-256 |
|---|---|---|---|---|---|---|
| FIG-07-02 | 《格子収束》 | 《未設定》 | 《未設定》 | 《未設定》 | 《未設定》 | 《未設定》 |

画像ファイルを手動で差し替えた場合も新たなハッシュを取得し、報告書Snapshotの証拠集合を更新する。配布可能な根拠データの範囲は機密管理に従う。

## 付録C：報告書とDB項目の対応

詳細は`db_report_field_mapping.md`を参照する。出力済みのCSVや集計表であってもSQLiteの`dataset`／`quantity`／`comparison`／`assessment`へ辿れることを確認する。DBに保存されないレビュー文章と正式承認記録は、報告書自体の版管理で保全する。

## 付録D：計画書の未確定事項O-01～O-09

| O ID | 未確定事項の概要 | 確定値／根拠 | 審査状態 | 完了日 | 影響するTest ID |
|---|---|---|---|---|---|
| O-01 | Mach/Re/α/βの対象範囲 | 未設定 | OPEN | — | GRID／WAKE／XCD／VAL |
| O-02 | 舵角範囲・基準舵角 | 未設定 | OPEN | — | CTRL |
| O-03 | Sref/bref/cref、モーメント基準点・単位 | 未設定 | OPEN | — | 全試験 |
| O-04 | QoI別の正式な要求精度 | 未設定 | OPEN | — | VAL／STA／CTRL／DRAG／ACC |
| O-05 | 既存パネル法／NS／風洞データ | 未設定 | OPEN | — | XCD／VAL |
| O-06 | 非定常参照データの有無 | 未設定 | OPEN | — | STA-03 |
| O-07 | 実施・レビュー・承認体制 | 未設定 | OPEN | — | 全試験／ACC |
| O-08 | 3.50.5実行バイナリ／API一致 | 未設定 | OPEN | — | 全実機Run |
| O-09 | Case／Run計画とCalibration・Validationの分離 | 未設定 | OPEN | — | 全試験 |

## 付録E：証拠の固定と再現手順

1. 採用した`study.json`、V&V計画書、`test_catalog.csv`、DBスキーマの版・SHA-256を記録する。
2. 解析Runの`requested_snapshot_json`、`effective_snapshot_json`、`generated_snapshot_json`を照合する。
3. 各Datasetの`source_description`、`metadata_json`、`method_id/version`を確認する。
4. 各Quantityについて`native_name/value`と`canonical_value/unit`、軸・基準値・換算定義を確認する。
5. `comparison`のPredicted／Reference Quantity ID、比較種類、不確かさ、判定を追跡する。
6. `acceptance_rule.approval='APPROVED'`の根拠が発行版に対応するか確認する。
7. `assessment`の技術レビュー者、承認日、使用制約、根拠IDを確認する。
8. 原本と生成物のSHA-256を検査し、発行時のDB Snapshot（原則読み取り専用）を保存する。
9. 本文と試験別詳細記録のTest ID・図表ID・Resultの値・引用元に不整合がないことを確認する。

DBの`run.status`、`comparison.verdict`、`assessment.disposition`の値は、本文の自由記述や表の空欄から推測で自動生成しない。レビューコメント・承認はDB外で記録する場合、文書版と署名を追跡する。

## 付録F：主要参考文献

下記の資料を方法論の参考とする。**規格の全文を無条件に実装・準拠したという主張はしない。** 本プロジェクトに実際に採用した事項は`methods_and_references.md`の対応表を参照する。

1. ASME, *V&V 20-2009 (R2021), Standard for Verification and Validation in Computational Fluid Dynamics and Heat Transfer*. https://www.asme.org/codes-standards/find-codes-standards/standard-for-verification-and-validation-in-computational-fluid-dynamics-and-heat-transfer
2. AIAA, *G-077-1998 (R2002), Guide for the Verification and Validation of Computational Fluid Dynamics Simulations*. https://webstore.ansi.org/standards/aiaa/AIAA0771998R2002
3. NASA, *NASA-STD-7009B, Standard for Models and Simulations* (2024). https://standards.nasa.gov/standard/nasa/nasa-std-7009
4. NASA, *NASA-HDBK-7009B, Implementation Guide for NASA-STD-7009B* (2026). https://standards.nasa.gov/standard/NASA/NASA-HDBK-7009
5. NASA NPARC Alliance, *CFD Verification and Validation Documentation Guidelines*. https://www.grc.nasa.gov/www/wind/valid/document.html
6. Oberkampf, W. L. and Trucano, T. G. (2002), *Verification and validation in computational fluid dynamics*, Progress in Aerospace Sciences 38(3), 209–272. https://doi.org/10.1016/S0376-0421(02)00005-2
7. Eça, L. and Hoekstra, M. (2014), *A procedure for the estimation of the numerical uncertainty of CFD calculations based on grid refinement studies*, Journal of Computational Physics 262, 104–130. https://doi.org/10.1016/j.jcp.2014.01.006
8. Celik, I. B. et al. (2008), *Procedure for estimation and reporting of uncertainty due to discretization in CFD applications*, Journal of Fluids Engineering 130(7), 078001. https://doi.org/10.1115/1.2960953
9. OpenVSP, *OpenVSP Verification & Validation Studies*, official source tag `OpenVSP_3.50.5`. https://github.com/OpenVSP/OpenVSP/blob/OpenVSP_3.50.5/examples/scripts/python_scripts/report.md
10. `plan.md`（本プロジェクト上位計画）
11. `database_spec.md`と`schema.sql`、`test_catalog.csv`（本プロジェクトのDB定義）

---

## 発行前の必須チェック

- [ ] 25件のTest IDが全て存在し、必要試験と条件付き試験の扱いが明示されている。
- [ ] 解析実行のRun成功と、数値・工学・用途の判定が区別されている。
- [ ] 使用した3.50.5バイナリ・API・入力形状・計算設定が固定されている。
- [ ] 基準量・座標軸・微係数の定義・CD0の包含項が確認されている。
- [ ] 数値収束と実験Validationの根拠、未評価項目の影響が明示されている。
- [ ] Calibratedデータを独立Validationへ転用していない。
- [ ] 承認済み要求精度と、正式な使用認定の手続が確認できる。
- [ ] 図表・比較値をDB SnapshotとArtifact SHA-256まで追跡できる。
- [ ] 最終的な適用可能範囲／使用制限がQoIと飛行条件ごとに明確である。
- [ ] レビュー者・承認者・発行版・改訂履歴が埋まっている。

**未完了項目が残る場合**：その影響が理解・承認できる段階的な報告書として公開することはあり得るが、VSPAEROの特定用途への正式な使用認定が成立したかのように表現してはならない。
