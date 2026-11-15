# TTF：Raabová 2007 / Pánková 2014の個体追跡IDの公開可能性を監査

**2026-11-15 | PR #155 (Draft) | 原著の公開Excelの識別列のみ検査 | fitness・適応度・遺伝・可塑性・種分化の応答値は未開封**

## 結論・現在地

**Pánková et al. (2014) のTable S1原本は実際に取得でき、野外植物の個体番号 `mark` が記録されている。** 同じExcelの`survival`と`performance`シートで、元の植物個体番号と移植先を使った記録照合が可能。一方、**Raabová et al. (2007) の個体番号を含む同じ型の原本をまだ特定できていないため、2007年対2014年の個体IDを1対1照合したとは主張しない**。

### 正規公開ソースとハッシュ

- 2014年原著：[Pánková, Raabová & Münzbergová (2014), PLOS ONE](https://doi.org/10.1371/journal.pone.0093967)
- 公式Excel：[Table S1, `pone.0093967.s001.xls`](https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0093967.s001&type=supplementary)
- ファイル形式：Excel旧OLE2/BIFF8、**147,456バイト**、取得原本SHA256 `7a8de0c203aec2b2edf0a7b6bceea649663a17095471a475fdb8a0abd69665c5`。公開ファイル自体はコミット・転載せず、ハッシュだけを固定。
- 研究データでなく**列の見出しと変数辞書**を読む段階：[Actions #43720495186](https://github.com/zuizui0223/TTF/actions/runs/43720495186)、原本メタデータJSON SHA256 `4a271585e9deb20f3d919e9f6a34a40a11ba1299b4bc90d1c79282a5678b8b2c`。
- 原本で確認できたシート名：`survival`、`performance`、`explanation`。説明シートが`mark`を**The number of the plant in the field**（現地の植物個体番号）、`nohelp`を**Ordering Numbers**（行の並び順）、`St`をTarget site、`So`をSite of originと定義する。なお、自動列名ヒューリスティックは`mark`を個体ID候補に認定しなかったが、**原著コードブックが直接その意味を確定したため、ヒューリスティックの陰性は「個体IDが存在しない」根拠にはならない**。

### 2014年の元Excel内だけでの個体番号照合

先に[許可列と対応キーを固定](./pankova2014_s1_original_plant_mark_join_contract_v01.json)したあと、`survival`と`performance`の`mark`・`St`・`So`列の**値だけ**を読み取った。行ごとの生存値・葉数・バイオマス・菌根感染率などの生物応答セルは読まず、原本の識別子・個体数・不整合数だけを計算した。

| メタデータ監査の単位 | 結果 |
|---|---:|
| `survival`の元個体行 | **547** |
| `performance`の元個体行 | **281** |
| `survival`内の`St×mark`重複 | **0** |
| `performance`内の`St×mark`重複 | **0** |
| 2シートで`St×mark`が一致 | **280** |
| そのうち`So`（原産集団）まで一致 | **279** |
| 2シートで`St×mark`一致だが`So`が食い違う | **1** |
| `performance`にしかない`St×mark` | **1** |
| `survival`にしかない`St×mark` | **267** |

元Excelの単純な`mark`だけでは`survival`に150種類・`performance`に144種類の値しかなく、**個体番号は移植先`St`とセットで扱う必要がある**。2014年の本体論文は生存判定母数547個体、最終の形質・バイオマスが測定可能だった281個体と報告しており、このメタデータの行数はMethodsの数と整合する。ただし、一部のID不整合は残るため**281/281完全追跡できたとは言わない**。原産集団の不一致1件とperformanceのみ1件は原データの差異として温存し、元の公表値やIDを勝手に修正しない。

- 正規のソースのみの実行：[Actions #43721024583](https://github.com/zuizui0223/TTF/actions/runs/43721024583)、原本ソースアクセス試験5件・XLS早期テキスト列検査8件・個体ID照合コード試験8件＝**21件全合格**。
- 本来のsource-only集計成果物SHA256：`7f30d8c34fc317090782cb6f3696e6166c4ad99746d49c315ba7d4e8b0fd1d13`。
- [保存した機械可読メタデータ集計](./source_registry/pankova2014_original_plant_mark_internal_linkage_v01.json)は上記Actionsの**2つの原始成果物ハッシュ**と個体の照合／不整合数を持つ。Excel本体、植物札番号の実値、元の個票・応答値はGitHubには保存していない。

### 実験コホート対個体の独立性

Raabová et al. 2007は2004年に開始した6原産集団×6移植先の成体相互移植実験を報告し、Pánková et al. 2014は原著Discussionで文献[21]（Raabová 2007）を「同じ相互移植実験の2年後の先行研究」と明示している。**成体移植の独立実験数は2論文合わせて1**。この確定は[著者の明示的引用に基づく受領書](./raabova2007_pankova2014_author_confirmed_adult_cohort_v02.json)と[論文×実験の重複排除表](./source_registry/raabova2007_pankova2014_experiment_report_dedup_v01.json)で記録済み。

しかし、本件の**2014年`survival`と`performance`の内部照合279件**を、**2007年と2014年の同一個体279件**と読み替えることは絶対にできない。両論文を跨ぐ個体タグの対応記録が公開資料でまだ確認できないからである。2007年の元個体番号テーブルが取得できれば、`St×mark`に年・実験の対応を足して**新しい事前固定照合**が可能だが、元情報がなければ重複人数・植物数はnull（未評価）のままにする。

**この段階での判定：`2014_EXACT_PLANT_MARK_CROSSSHEET_PARTIAL_LINKAGE_279_OF_281_NO_2007_TAG_CROSSWALK`。** Bontrager 44/50種の停止、PROCEED 56種、既存のTTF遺伝・TTF-Qアウトカムは全て不変。PR #155はDraft、`main`への反映はしない。
