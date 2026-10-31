# VSPAERO 3.50.5 V&V — 現行実装状況（v20261031.11）

**時点：2026-10-31**。本書は開発時点のソース実装と検証範囲を記述し、数値的・工学的V&Vの合格を認定しません。正規の要求事項はV&V計画書v1.1、詳細DB仕様v1.0.1、レビュー済み設定で確認します。過去の各版での主張は`CHANGELOG.md`に分離します。

## 1. 共通機能（現行コード）

| 領域 | 実装状況 | 実証境界 |
|---|---|---|
| Python APIによるComputeGeometry・VSPAEROSweep | 実装済み | Windows実機成功はv20261031.3のGEO-03-G3とv20261031.10のGEO-03-G3等で観測。v11修正の実機検証は未実施 |
| 解析モード別のネイティブ結果確認 | 実装済み | `.polar`、`.stab`、`.adjoint.stab`、`.pstab/.qstab/.rstab`を分離。STA系は実機未検証 |
| Case生成・登録・再利用・逐次Run | 実装済み（設定・Case依存） | 全予定Caseを自動生成できるわけではない |
| SQLite Run/Point/Dataset/Artifact追跡 | 実装済み | 新規DBのみ、旧DB Migrationなし |
| SHA-256、実行ソース指紋、結果監査 | 実装済み | Native終了コードはPython APIに非公開 |
| native履歴/安定微係数・CD0取込、比較、GCI候補 | 実装済み（一部限定） | 実際の単位・符号・メッシュ幅や基準条件の承認が必要 |
| CLI`execute-all`・CSV/PNG/13章Markdownレビュー報告 | 実装済み | 不足証拠を表示し人間の合否審査は行わない |
| VSPAERO NCPU指定（`--ncpu`／Case JSON） | 実装済み | API要求・読戻し・SQLite指紋・Run再利用をオフライン検証。v20261031.10 Windows実機でNCPU=8のGEO-03-G3成功、NCPU=4の一括実行実績あり |
| CLI進捗表示（標準・`--verbose`） | 実装済み | Case進捗・60秒生存確認・Wake観測、JSONL記録。v20261031.10 Windows実機で動作確認、v11の追加観測項目は未検証 |
| 新しい`src/openvsppy`レイアウト | 実装済み | Notebook 13件は`sys.path`方式、Editable Install任意 |
| ドキュメントの破壊的統合 | v20261031.6で実装 | 核心6文書の内容・SHA-256保全、計35→9文書 |

v11ではControl Group/Native Mesh事前検査・Native異常分類・Adjoint詳細進捗・Case別タイムアウトを追加。Windows実機での修正後成功は未検証。本版のオフライン回帰試験は別途`pytest`で確認します。前版v20261031.9では105 PASS / 13 DESELECTED（標準実行からNative API、VSPAEROソルバー、GUIを分離）。新しいLinux版のスタンドアロンVSPAERO診断スモークは**1 PASS**ですが、正式なGEO-03-G3やWindows Python API経由のNative試験は未実施です。

## 2. Test ID別の実装・証拠充足状況

以下の「共通経路あり」は対象Caseが承認・生成された場合の汎用パイプラインを指し、各試験の公式比較量・実機完走を意味しません。

| Test ID | 現行の実装範囲 | 未完了または要審査 |
|---|---|---|
| VER-01 | 候補ベンチマーク生成・実行経路 | Hershey Bar理論分布の比較、実機証拠 |
| VER-02 | 候補楕円体・格子設定・実行経路 | 理論Cp分布とネイティブ出力の比較 |
| GEO-01 | Thin/Thick構成比較の共通経路 | 実効Setの証拠・試験固有評価 |
| GEO-02 | 構成要素Caseの共通経路 | 干渉寄与の定量評価 |
| GEO-03 | G3を含むCase生成と共通経路 | v3のG3のみ実機成功、他Case未検証 |
| GEO-04 | モデル化切替の共通経路 | 実際のThin/Thick差分の検証 |
| GRID-01 | 3格子比較・GCI候補 | 実効格子幅と漸近範囲の審査 |
| GRID-02 | 第二条件の格子感度経路 | 実機比較・品質確認 |
| GRID-03 | LECluster/NRef等のCase候補 | 局所分割・胴体局所メッシュ検証 |
| WAKE-01 | Wake反復感度と履歴処理 | 収束性の工学的判断 |
| WAKE-02 | Wake節点数感度 | 全Caseでの実機検証 |
| WAKE-03 | Fixed/Implicit候補経路 | ソルバー方式の実効検証 |
| XCD-01 | 既存パネル法の参照Dataset取込・比較 | 参照データ、分布・座標変換 |
| XCD-02 | NS参照Dataset取込・比較 | Cp分布比較・位置補間・物理条件一致 |
| VAL-01 | 風洞参照Datasetの条件照合・比較 | 対象風洞データと不確かさの審査 |
| VAL-02 | 横・方向風洞データとの比較経路 | 参照・軸系・符号の承認 |
| STA-01 | `.stab`および`.adjoint.stab`取込 | 物理微係数への変換・実機検証 |
| STA-02 | 微係数の格子/Wake変化取込 | 数値感度の実機比較 |
| STA-03 | `.pstab/.qstab/.rstab`固有Parser | 実効ステップ・振幅・周期性・実機検証 |
| CTRL-01 | Elevator中心差分の共通経路 | Control Surface Gainと物理符号の実機比較 |
| CTRL-02 | Aileron差分の共通経路 | CMx符号と反転の回帰・実機比較 |
| CTRL-03 | Rudder差分の共通経路 | 横・方向係数の符号と偏向条件の実機比較 |
| DRAG-01 | ParasiteDragと方法A | CD0参照条件・二重計上審査 |
| DRAG-02 | 方法Bによる補正ポーラー | ベースCDoと適用範囲の承認 |
| ACC-01 | 25件の証拠一覧・13章報告書 | 工学的使用認定は人間が行う |

## 3. 状態と検証証拠

- `SUCCEEDED`：所定のファイル・数値・SQLite連結条件を満たしたワークフロー成功であり、V&V合格ではない。
- `NATIVE_ONLY`：ネイティブ値を保管したが、航空機軸系・符号・正規化が未承認。勝手に規範的な値へ変換しない。
- `BLOCKED`/不足：必要な条件・参照データ・Caseがない場合。補完値を仮定しない。
- 旧版のテスト報告は**その版の証拠**でありv6の新規実機成功ではない。直近のオフライン回帰は本版のテスト結果を`CHANGELOG.md`に追記する。

## 4. 主な残課題

1. **実機検証**：v6のGEO-03-G3、STA-01～03、全`execute-all`。
2. **数値Verification**：実メッシュ幅、後流反復収束、理論との対応、Cp/翼幅荷重の分布Parser。
3. **Validation**：XCD/VAL参照データ、基準系・条件・不確かさを用いた比較。
4. **安定・操舵**：P/Q/R実効数値設定、符号・正規化・物理軸の人間による承認。
5. **報告・承認**：Markdownはレビュー草稿。Word/PDF/XLSXの直接生成や自動的な使用認定は提供しない。

## 5. 正本と実装の対応

| 対象 | 場所 |
|---|---|
| 計画書 | `docs/vv/plan.md` |
| DB仕様 | `docs/vv/database_spec.md` |
| 報告書テンプレート | `docs/vv/evaluation_report_template.md` |
| 解析条件 | `vv_study/*.reviewed.json`、`vv_study/execution.reviewed.json` |
| ソース | `src/openvsppy/vv/` |
| SQL/Catalog | `vv_study/01_database/` |
| 状態の実記録 | `vv_study/01_database/*.sqlite`、各Run証拠 |
