# TTF｜Raabová 2007–Pánková 2014 個体IDデータ探索：二次アーカイブまでの実測監査と終端判定

**2026-11-16 | Draft PR #155 | source-only、適応度未開封 | Bontrager 44/50 種停止判定不変**

## 今回の到達点

Raabová et al. 2007 と Pánková et al. 2014 は**著者自身が同じ2004年開始の成体相互移植実験だと明記している**。独立した物理実験として数える単位は論文2本に対して**成体移植1コホート**。個体番号の対応は別問題。

2014年に公開されたPLOS ONE Table S1の正規原本は旧Excel形式の147,456バイト（SHA256 `7a8de0c203aec2b2edf0a7b6bceea649663a17095471a475fdb8a0abd69665c5`）。説明シートは `mark` を「野外での植物個体番号」と明記。適応度値を読まずに、2014年内の`survival`と`performance`の`St × mark`キー**280件**を接続し、うち**279件**で原産集団`So`まで整合。原産地の食い違い1件とperformanceだけにある1キーは補正しない。[2014元Excel内部監査](./source_registry/pankova2014_original_plant_mark_internal_linkage_v01.json)参照。

**2007年側の対応する個体札番号列はまだ未取得**。しかし「個体が何本同一だったか」を知らなくても、同一物理実験を2論文で独立2回と誤計数しない判定は既に成立する。

## 新規の二次公開アーカイブ発見と検証

[Yates & Fraser (2014)](https://doi.org/10.1111/eva.12181) は、J. Raabováからデータの提供を受けたことをAcknowledgementsに明記している。彼らの公開データは[Dryad DOI 10.5061/dryad.81k4v](https://doi.org/10.5061/dryad.81k4v)であり、公式[Zenodo保存記録 #4963190](https://zenodo.org/records/4963190)に**DOI完全一致の保存コピー**がある。原本ファイルは`Source.Population.Size.Performance.In.Novel.Environments.xlsx`、**162,038バイト、公開MD5 `31753f59df54a13ca8fbc992a348702d`**。

ファイルリスト・DOI認証：[GitHub Actions #43874024889](https://github.com/zuizui0223/TTF/actions/runs/43874024889) **成功、合成7テスト成功**。原本Excelの完全な公開MD5照合・実測SHA256 `da3aa332c84fc0d0e2d79b8663563d51e79baee7d9927d87231bf8195550bd18` の後、**シート名と冒頭テキスト見出しだけ**を抽出した。実行：[GitHub Actions #43875572888](https://github.com/zuizui0223/TTF/actions/runs/43875572888) **成功、合成6テスト成功**。

| 元Excelの4シート | 初行のデータ範囲・項目 | 元の2007個体札番号の明示列 |
|---|---|---|
| `Home vs. Away` | Taxa, Population, Transplant.Site, Species, Home/New Survival, Effect.size, LOR, SE, Article | 未確認（冒頭見出しにはなし） |
| `Common Garden ` | Population, Transplant.site, Home, Population size, Survival, Dead, Sample.Size, Article | 未確認（冒頭見出しにはなし） |
| `Reciprocal Transplant` | Population, Transplant.site, Home, Population size, Survival, Dead, Sample.Size, Article | 未確認（冒頭見出しにはなし） |
| `References` | 元文献一覧 | 該当なし |

これは元研究者から提供された**集団・移植先単位の情報をメタ解析に変換した二次公開データ**である。元論文で解析された規模は13種・111集団・874生存率推定値。Excelの`Population`を2004年の植栽個体`mark`として扱わない。**新しい独立したRaabová成体移植実験ではない**。

原本由来の列見出しだけを[再現可能なJSON台帳](./source_registry/yatesfraser2014_raabova_secondary_original_xlsx_schema_v01.json)に保存。元Excel全体、数値的な生存率・効果量・分散、植物の生IDはGitHubに転載・解析していない。

## カレル大学公式記録の出典品質

- [2008年博士論文の公式記録](https://dspace.cuni.cz/handle/20.500.11956/94265?show=full)では添付資料として本文・要旨・審査記録の**PDF**のみ掲示されている。
- [2009年の学位論文の公式記録](https://dspace.cuni.cz/handle/20.500.11956/95341?show=full)でも本文・要旨・審査記録の**PDF**のみ掲示されている。
- GitHub Actionsから同サイトを実アクセスした際は2件ともタイムアウト、2007年のElsevier出版社ページはHTTP 403。[先行アクセス記録](./raabova2007_individual_tag_source_access_receipt_v01.json)参照。

**添付ファイル一覧にXLS/CSVがないことは、著者の元の手元の個体札データが存在しないという意味ではない**。学位論文PDFにも個体札の生データが表として載っているかどうかは、現時点では確認していない。

## TTFの次工程判定

**このRaabová2007→Pánková2014の「2007年個体番号」探索は、今のTTF主線ではここで停止が妥当。** 理由は、元個体表が回収できなくても**成体2004年コホートの重複排除は確定済み**で、2007個体を追加取得しても新しい独立した移植実験を増やせないからである。

原著提供者への非公開元データ依頼は実施しておらず、依頼が必要な場合は個体札番号・原産地ID・移植先ID・観測年の**識別子だけ**を対象とする。

次の科学的主目標は、歴史気候・局所適応に独立した情報を付加できる**新しい原著実験と原産集団の測地情報**を、応答を見る前に認証すること。これまでの35のHargreaves未照合研究ラベルのうち21研究には相互移植先由来の代理原産地座標があるが、真の原産地GPSとして認証したものではない。旧Bontragerの44/50種ゲートに後付けで加えて救済しない。

**最終source-only状態：`DERIVATIVE_META_ANALYSIS_RECOVERED_ORIGINAL_ROW_LEVEL_2007_TAG_NOT_CERTIFIED_2004_EXPERIMENT_ALREADY_DEDUPLICATED`**。Bontrager、PROCEED、TTF-Q・遺伝解析の凍結判定を保全。PR #155はDraft、`main`未変更。
