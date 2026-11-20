# Albano 2025 原本ZIP：96原産集団の座標と実験単位を実データから初認証

**2026-11-19 | Draft PR #161 | ORIGINAL SOURCE-ONLY | ユーザーアップロードの真正ZIPをローカル監査、原本CSV／適応度データはGitHubに未登録**

## 1. 公式原本とのバイト単位照合（確認済み）

今回のアップロードファイルは `doi_10_5061_dryad_bk3j9kdqk__v20250905.zip`（ZIP自体451,871バイト、SHA256 `855451829407ebbed5560aa01f2b02cdee78fa1ca529cd10d383a1b5e320d80a`）。ZIP内部の真正原本 `Albano_et_al_2025_EM_Experiment_Data.csv` は**421,805バイト**。バイト列からローカルで計算したSHA256が

`b3796db37145dd32ab148f5c2105ca0646923feab153d8857d4f90a5539f080c`

であり、事前凍結していた**Dryad公式APIのdataset `10.5061/dryad.bk3j9kdqk`、version `394175`、file ID `4337649`の掲載SHA256と完全一致**した。前回のCI上の401/403アクセス失敗と違い、**元CSVの真正バイトデータを取得・検証できた**。

ZIPにはこれ以外に `README.md`、現代気候の `WorldClim_z-standardized.csv`、気候PCAスコアCSVの計4ファイルがある。今回、他の気候CSVから21kaの値を推定していない。READMEは気候データの行順を記述するが、PCA表の各原産集団行には直接IDラベルがなく、そのまま新しい座標ごとの歴史気候値とみなさない。

## 2. 実データsource列の監査結果（計算・確認済み）

元CSVの35列のうち、`Garden, Position, Population, Rep, ContOrigin, ContGarden, Latitude, Longitude` の**8つの原産地／実験ID列だけ**を抽出・集計した。生存率、成長率、結実、花数、HCN遺伝子型、温度適応度等の応答列の値は統計解析、出力、モデリングをしていない。

| 原著source-only指標 | 実CSVの値 | 解釈 |
|---|---:|---|
| 元個体レコード | **1,858行** | 元々の共通圃場実験の植物レコード。研究数ではない |
| `Population` 個別ID | **96** | 96原産集団ラベル |
| 緯度経度の完全一致した独立したペア数 | **96** | 全96集団でIDごとの1座標ペア、ID間の座標完全重複0 |
| 集団内の座標不一致 | **0** | 同じ集団の別の移植先でも座標一致 |
| 原産大陸 | **NAM 47 / EUR 49集団** | 原産地の構成 |
| 移植先 | **4園** | Lafayette 465行、Mississauga 465行、Montpellier 470行、Uppsala 458行 |
| ID×gardenの実測組合せ | **372 / 384** | 存在しない12組合せを「0適応度」として補完しない |
| 4園すべてに現れる原産集団 | **91** | 他3集団は2園、他2集団は1園のみ |
| `Garden + Position` | **全1,858行一意** | 個体を識別できる合成キー |
| `Garden + Population + Rep` | **107重複キー（余剰131行）** | `Rep` だけでは個体を一意に識別できない |
| 同大陸内の最近接原産地2点 | **0.388079 km** | 1km未満の別の集団IDが1組ある |

注意：**異なる位置座標が96個あるだけでは、96個の独立したLGM気候格子があると証明したことにはならない**。原産集団の「Latitude,Longitude」という真正のCSV記録と、著者が採集時に個別GPSを測量したかどうか・どの精度かという測地出典の認証を区別する。96の位置ペアから歴史気候に意味のある独立情報があるかは、21kaの気候ラスタを定義してから応答盲検で別途検定する。

またREADMEは最初の生育季節を北米で2020年・欧州で2021年と記しているので、**庭・大陸と年の交絡**を扱う必要がある。植栽個体の多さが大陸間の実験年を独立複製するわけではない。原産地のNAM/EURと庭大陸のAmerica/Europeは別フィールドであり混同しない。

## 3. 成果物の保管と再現

- **GitHubには原本CSV・ZIP・96座標の行データをコミットしない**。実データの座標ごとのsource-onlyレジストリは今回のユーザーセッション内に生成済み。
- 本PRの `source_registry/albano2025_verified_source_only_receipt_v01.json` は全96行の座標そのものを含まない**集計・原本ハッシュの受領書**。
- `albano2025_verified_archive_source_only_v01.py` は指定の元CSVをZIPから唯一名で抽出、サイズ/SHA256を厳密チェックし、source列だけで36等の高次元応答は解析せず、集団ID・大陸・庭の整合性、重複Rep ID、欠落組合せ、空間近接を検査する。
- `test_albano2025_verified_archive_source_only_v01.py` は原本を必要としない合成ZIP上でfail-closedな検査を実行する。`test_albano2025_frozen_source_receipt_v01.py` は公式Dryadハッシュ受領書と今回の集計を照合。CIでは**元ZIPを持たないため生データそのものを再ダウンロードしたと主張しない**。実ZIPは本セッションで実際に確認した。

## 4. 科学的判定と次の主線

**SOURCE_ID→COORDINATE MAPPING AUTHENTICATED IN ORIGINAL BYTES: PASS.** 旧10研究の最大の障害であった、元ZIPへのアクセスと元ID×位置の真正データ照合を**Albano 1研究で克服**した。

**NOT YET A NEW INDEPENDENT LGM EFFECT PANEL.** 1植物種・4共通圃場は4独立研究ではなく、独立LGM格子・気候履歴・種間一般性は未検証。現代気候・地理的近接性・大陸／導入史・庭と実験年の分離を応答盲検で検討すること。必要な今後の対比は最低でも`21ka历史気候 vs 現代気候 + 座標由来の地理項 + 大陸 + garden/年`であり、どの気候変数セット・景観距離・分割・停止判定を使うかは**応答を見る前に新規契約として固定する**。

Bontrager 44/50停止、旧TTF-Q、211種の遺伝予測力のnegativeは変更しない。投稿可能な「LGMが適応／可塑性／種分化を予測」とまだ主張しない。PR #161はDraft、main未変更。
