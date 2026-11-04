# CHANGELOG — OpenVSP/VSPAERO V&V（統合版）

本書は旧`docs/vv/`内のRelease Notes、実装・仕様Matrix、Windows実行ガイド、過去のテスト結果など**29件を破壊的に統合した開発記録**です。過去の操作コマンドは当時の構成を示す記録であり、現行CLIではありません。現行操作は`cli_guide.md`、現行実装状況は`implementation_status.md`を参照してください。

## v20261104.0 — Markdown評価報告書の正確性・可読性・可搬性改善（2026-11-04）

- 報告書に明示マーカーを導入し、5.1節の静的「全件未実施」表を固定SQLite／採用Runの状態表へ置換。Run成功、数値判定、使用認定は分離したままとする。
- 正式13章の`evaluation_report.md`を索引・代表図中心に整理し、`test_reports/<Test ID>.md`に全Quantity、Dataset/Method、メッシュ/履歴、参照比較を保持。元CSVは従来どおり全Quantityを出力。
- 図の横軸をCTRL指令角・実測メッシュ数・Wake条件等の物理量とし、GEO等のカテゴリは連続線を使用しない。同名QoIのMethod／Dataset違いを上書きせず分離。主たるQoIを試験群ごとに限定し、L/Dと微小係数を同一棒図へ混合しない。
- 図の相対パス・Figure IDと出典の`figure_manifest.csv`を追加。`results_snapshot.json`を`REPORT_SNAPSHOT_V2_RELATIVE_PATHS`に拡張。
- ライブDBからの生成データを登録した後にDBを固定し、**本文と表CSVは固定DB**から生成。DB自身と報告書自身のハッシュは自己参照を避け`evidence_manifest.json`で管理。
- Markdown GFMテーブル構文、物理横軸、複数Dataset、全25 Test、固定DB、移動後の画像パスなどのオフライン回帰テストを追加。
- **既存Run、SQLiteスキーマ、V&V計画書の数値要求は変更していない。** 旧版の出力とは表示形式・図の点数が異なるため、既存報告書や解析原本を上書きしない。Windows/OpenVSP 3.50.5での新規ソルバー解析は本版では未確認。

## 文書体系の整理（2026-10-31）

- `docs/vv/`の運用ファイル名を固定し、計画書・DB仕様・評価報告書雛形もファイル名から版数を除去。文書内の正式な版数と既存の技術的な試験条件は維持。
- 現行パスの参照をコード・設定・テスト・文書で更新。`README.md`を日本語で再構成し、一度限りの比較スクリプトを現行ツリーから削除（過去コミットには残存）。
- 上位計画書は文書間参照のファイル名だけを修正したため、元のSHA-256 `31a30e103fbdbc9268b9a390cf2b4d8b29df6f611f25fb938a8b08fd150b7130` から `d699486564082163c199ebb5c030832727cef07524d36e9f41c6dee5df90143d` に変化。新規Studyの`plan.sha256`は更新済み。**過去のStudy・Runが参照した旧ハッシュを改変しないこと。**
- WindowsのGitチェックアウトで計画書とネイティブ証拠のバイト列が書き換わらないよう、`.gitattributes`で改行変換を制限。旧版の表記・記録は履歴として保持。

## 最近の変更

### v20261031.11 — 一気通貫のV&V実行を優先したNative障害対応

- 実際のThin/Thick Geometry Setに合わせてControl Groupを無効化し、部分的に欠落したGroup・不明な親形状・非ゼロ偏向の欠落Groupを明確に拒否。承認済みGain・偏向符号を保全。
- `VSPAEROComputeGeometry`後に`.vspgeom`・`.vkey`の構造検査を追加。異常なNRef等では入力を勝手に修正せず、`mesh_diagnostics.json`に証拠を保存してソルバー起動前に停止。
- 結果ファイル未生成時は監視ログの最終観測段階と異常分類を`api_diagnostics.json`へ保存。Python APIで得られないソルバー終了コードは推測しない。
- AdjointのForce/Moment段階・GMRES反復を`--verbose`で可視化。Mode/Case別のタイムアウトを`execute-all`のマニフェストから指定可能（既定`ADJOINT=10800s`）。
- 本体のVSPAERO呼出しは引き続きPython API。Wake反復数・NCPU・制御舵面Gain・SQLite schema・V&V計画書は変更しない。
- Windowsで観測済みの障害に対応する回帰テストを追加。オフラインの回帰結果とWindows実機でのNative再実行結果は区別する。

### v20261031.10 — VSPAERO NCPU実行時指定

- `run`、`run-test`、`execute-all`に`--ncpu`を追加（正の整数のみ）。Case JSONの`numerical_setting.requested.vspaero.NCPU`も認識し、省略時は従来どおり1。
- 明示したCLI値はRun実行時のみ優先され、レビュー済みCase/setting_id、G103Aモデル、Study/DB Schema、V&V計画書は変更しない。
- OpenVSP API `VSPAEROSweep/NCPU`設定・読戻し、診断JSON、標準進捗、Run指紋・SQLite環境記録を一致させる。
- 再利用時はNCPU一致を必須化。過去版Run（NCPU記録なし）は安全側に再利用しない。旧Runを消去/変更しない。
- APIスタブ・CLI・DB・再利用制御のオフライン回帰試験を追加。**123 PASS / 13 DESELECTED**、Linux/Windows 3.50.5実機での新NCPU値は未検証。

### v20261031.8 — テスト構成の整理（2026-10-31）

- Release-date based `test_vv_20261031*.py` and monolithic V&V test files were reorganized into 9 feature-based files.
- `tests/conftest.py` now provides a single isolated Study/SQLite fixture factory and seeded Case fixture; equivalent duplicate test setup was removed.
- `test_vv_repository_layout.py` checks required documentation and notebook bootstrap semantics without freezing expandable file counts.
- Test coverage and native evidence fixtures are retained. No production Python code, approved G103A input, SQLite schema or V&V plan is changed.
- Offline baseline: 103 PASS / 8 SKIP (same test population as v20261031.7). Native Windows VSPAERO execution NOT performed in this release.


## リリース別の要点

| 版 | 変更の概要 | 検証範囲・留保 |
|---|---|---|
| 2026-10-30 | P9–P14に係るRun再利用・指紋・理論/格子/Wake/安定・参照比較/DRAG・報告基盤、Windows winfixを整理 | 旧DB移行・実行手順は現行版では廃止。実機未確認項目が多い |
| v20261031.0 | `execute-all`、一括集計、PNG図、13章Markdownレポートを追加 | 25 Test ID全件の条件が揃うわけではない |
| v20261031.1 | VER候補、GEO-04、GRID-03、STA-03、DRAG-01/02接続とレビュー済み条件の強化 | VER理論分布/実機検証は未完 |
| v20261031.2 | WindowsのControlSurface Group読戻し、Hershey Barスパン、VER-02格子、非有限値取込を修正 | GEO-03-G3のAPI実行が `.polar`なしで失敗した実機記録あり |
| v20261031.3 | ワーカーのRun-local cwd、VSPAEROSweep `RedirectFile`、API診断JSON、厳密なネイティブファイル検証 | **Windows実機でGEO-03-G3が20/20反復完了、SQLite `SUCCEEDED`登録**（下記） |
| v20261031.4 | 解析モード別`.stab`/Adjoint/P-Q-R取込、ネイティブ値区分、DB原子的登録、完全なControlSurface補助ファイルハッシュ | オフライン試験。STA系Windows実機未検証 |
| v20261031.5 | `src/openvsppy`へ破壊的パッケージ移行、CLI/worker新モジュール名、Notebook `sys.path`、新規DB専用化、`docs/vv`集約 | 94 PASS / 8 SKIPのオフライン試験、Windowsネイティブ未検証 |
| v20261031.6 | `docs/vv` 35→9文書に統合。6件の仕様・雛形をそのまま保持し、運用・実装状況・履歴の3件を再構成。README参照を修正 | 文書・参照・回帰テスト（98 PASS / 8 SKIP）、CLI新規DB 25 Test ID・68 Case、13章報告書生成。Windowsネイティブ実行はしない |
| v20261031.9 | pytest既定のオフライン分離、Native API/VSPAERO/GUIの明示マーカー、共通fixtureと作業領域、個別プロセス時間制限ランナー、LinuxスタンドアロンSolver Smokeを追加。pytest.ini廃止 | オフライン105 PASS/13 DESELECTED、Linux vspaero実機スモーク1 PASS。Python APIを要するNative試験とWindows API経由Runは未実施 |
| v20261031.8 | 版番号付きテストの機能別整理、Study fixture共通化 | オフライン103 PASS/8 SKIP。Windows APIネイティブ未検証 |
| v20261031.7 | `execute-all`に標準stderr進捗、`--verbose`のWake・出力・DB情報、60秒heartbeat、共通JSONLイベントを追加。API実行条件は維持 | 103 PASS / 8 SKIP。子プロセス実行でWake・heartbeat・timeoutを検証し、CLIでJSON/進捗分離・68 BLOCKEDを確認。Windows 3.50.5ネイティブ解析は未検証 |

## 2. 2026-10-31 Windows API成功の確定証拠（v20261031.3）

- Case：`GEO-03-G3`、Run：`RUN-9691c15aea824527`、Point：`PT-043e039f51f14960`。
- Windows実測子プロセス：`vspaero.exe -omp 1 <Run raw model_effective>`、親プロセスは`python.exe -m vv.worker`（v3の旧モジュール名）。同じソルバーSHA-256：`f8cc76fc7f14f7ce284bf70ca410039dbe8749b7c155d26c3b75b4c05ebfc98d`。
- Wake反復20/20、`Total setup and solve time: 457.932000 seconds`、監視ログ末尾`Done`。
- 出力：`.polar` 3,025 B、`.history` 20,587 B、`.lod` 185,339 B、`.adb` 1,895,204 B。SQLite `SUCCEEDED`、全Artifact `VERIFIED`、レポート生成を実機コマンド出力で確認。
- API Run第20反復（`.history`）：CLtot=0.509707363631、CDtot=0.018307937794、CMytot=-0.05487592527、CLiw=0.499464650875、CDiw=0.005219332193。
- 直接起動対照（`.polar`）：CLtot=0.509706605733、CDtot=0.018312180664。入力の`.vspgeom`/`.vspaero`/`.csf`はSHA-256で一致。数値差は小さいが、 `.vsp3`/`.vkey`のバイト同一性は確認されなかった。
- 差異はv3にて統合した作業ディレクトリと`RedirectFile`などの変更のどれが決定要因だったかは**未特定**。再現性・ネイティブ終了コード取得は別問題。
- 本成功は**解析実行/DB登録の成功**であり、V&V工学合格・後続版の実機成功を証明しない。

## 3. 3.50.5ネイティブ制御・不具合の経緯

- 公式3.50.5ソース：`VSPAEROMgr.cpp`は`VSPAEROSweep`を通じ`ProcessUtil`でプロセスを管理。Windowsの`CreateProcess`、ブロッキング`ReadFile`と`WaitForSingleObject`を使用する。ソルバーのネイティブ終了コードはPython APIに返されない。`ERROR 7: Could not open Polar file`だけで終了原因は特定不能。
- v2の失敗：`Updating wake grids...`付近で監視ログが終わり、`.polar/.history/.lod/.adb`が見つからなかった。**直接実行**では20反復・正常結果を取得し、形状・Solver入力だけでは説明できないことを確認。
- v3の出力監視ファイル`raw/vspaero_api_stdout.log`を使用すると途中進捗が観測できる。解析実効モードはGEO-03-G3=`STEADY_SWEEP`/`STABILITY_OFF`、STA-01 Steady=`STABILITY_DEFAULT`、取り違えない。
- Control Surface Gainは承認済み設定（Aileron +1/+1、Elevator +1/-1、Rudder -1）を維持。変更を伴う障害回避は採らない。
- `-pstab`/`-qstab`/`-rstab`にはC++実装由来の128 time steps、5°振幅、25 Wake反復、32ノード等の設定がある。**設定要求値とソルバー内部実効値は区別する**。
- 非有限なネイティブ値は勝手に0に変換しない。微係数は軸系・正規化の承認がない限り`NATIVE_ONLY`で保管する。

## 4. 旧テスト・実機診断の扱い

- v20261031.6はオフライン98 PASS / 8 SKIP（文書整合性テスト4件を追加）。新規SQLiteは25 Test ID・68 Case、13章報告書生成を確認。
- v20261031.5はオフライン94 PASS / 8 SKIP。8 SKIPはOpenVSPネイティブAPIやGUI条件等による。STA/P/Q/R実機成功の根拠ではない。
- v20261031.4は91 PASS / 8 SKIP、v20261031.3は83 PASS / 8 SKIP（各版固有のテスト）。**版間でテスト母集団は異なる**ため数だけでソフト品質を比較できない。
- 2026-10-30のWindows回帰では13件の失敗修正（Group読戻し、スタブ期待値など）を試みた履歴がある。旧版の失敗事象と対策は現在版の再現テスト証拠ではない。
- GEO-03-G3以外の25試験に対し、必要なXCD/VAL参照データ、理論分布、各モードの現物出力などは依然欠けており、レポートが生成できてもValidation PASSを意味しない。

## 5. 統合元の文書台帳

以下はv20261031.5の同梱文書から統合したものです。**本文をそのまま連結しておらず、現行コードと整合する事項および歴史的な事実に整理しました**。当時は正本6件のハッシュを保って残しました（現在の計画書は上記のとおり参照先の整理でSHA-256が変化しています）。元の全内容が必要な場合は配布元`OpenVSP-main_v20261031.5.zip`を参照してください。

| 統合元（旧`docs/vv/`基準） | 統合先 |
|---|---|
| `RELEASE_NOTES_v20261031.5.md` | `CHANGELOG.md` |
| `REVIEWED_20261031_V2_README.md` | `cli_guide.md` |
| `VV_CLI_Guide.md` | `cli_guide.md` |
| `VV_FULLVV_SPEC_MATRIX.md` | `implementation_status.md` |
| `VV_IMPLEMENTATION_MATRIX_20261031.md` | `implementation_status.md` |
| `VV_IMPLEMENTATION_STATUS.md` | `implementation_status.md` |
| `VV_LAYOUT_SPEC_MATRIX_v20261031.5.md` | `implementation_status.md` |
| `VV_SPEC_IMPLEMENTATION_MATRIX_20261030.md` | `implementation_status.md` |
| `VV_SPEC_MATRIX_20261031.1.md` | `implementation_status.md` |
| `VV_SPEC_MATRIX_20261031.2.md` | `implementation_status.md` |
| `VV_SPEC_MATRIX_20261031.4.md` | `implementation_status.md` |
| `history/EXECUTE_ALL_GUIDE_20261031.md` | `cli_guide.md` |
| `history/FULLVV_CHANGELOG.md` | `CHANGELOG.md` |
| `history/RELEASE_NOTES_v20261030.md` | `CHANGELOG.md` |
| `history/RELEASE_NOTES_v20261031.0.md` | `CHANGELOG.md` |
| `history/RELEASE_NOTES_v20261031.1.md` | `CHANGELOG.md` |
| `history/RELEASE_NOTES_v20261031.2.md` | `CHANGELOG.md` |
| `history/RELEASE_NOTES_v20261031.3.md` | `CHANGELOG.md` |
| `history/RELEASE_NOTES_v20261031.4.md` | `CHANGELOG.md` |
| `history/TEST_RESULTS_20261031.md` | `CHANGELOG.md` |
| `history/TEST_RESULTS_FULLVV.md` | `CHANGELOG.md` |
| `history/TEST_RESULTS_v20261031.1.md` | `CHANGELOG.md` |
| `history/TEST_RESULTS_v20261031.2.md` | `CHANGELOG.md` |
| `history/VV_CLI_Guide_pre_v20261031.5.md` | `CHANGELOG.md` |
| `history/WINDOWS_API_DIAGNOSTICS_v20261031.3.md` | `CHANGELOG.md` |
| `history/WINDOWS_FIX_13_FAILURES.md` | `CHANGELOG.md` |
| `history/WINDOWS_FULLVV_GUIDE.md` | `cli_guide.md` |
| `history/WINDOWS_GUIDE_v20261031.1.md` | `cli_guide.md` |
| `history/WINDOWS_GUIDE_v20261031.2.md` | `cli_guide.md` |

> この台帳は文書移動の追跡を目的とします。統合先の記述は、最新仕様と履歴を明確に区別しています。削除前の元文書は今回のZIPには再収録していません。
