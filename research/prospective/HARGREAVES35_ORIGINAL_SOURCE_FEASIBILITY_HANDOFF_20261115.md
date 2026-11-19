# HargreavesのBontrager未照合35研究：原産集団の地理情報を回復できるか

**2026-11-15 | TTF Draft PR #155 | response-blind source-only handoff | fitness・Genetic・speciation応答は未開封**

## 判定（今回の新規結果）

Bontragerの**全植物原本**とHargreavesの全植物原本を先に照合して確定した「研究ラベルが一致しない35件」について、Zenodoに保存された**真正Hargreaves原本のメタデータのみ**を再監査した。

- 35研究すべてに植物行が存在（合計**1,234行**）。これは独立した生物学的反復数ではない。
- 35研究すべてに原産集団IDが2種類以上、31研究には移植先IDが2種類以上存在。
- **30研究・167研究×種×実験ID×移植先×fitness_typeセル**で、同じ移植先・実験・指標内に異なる原産集団IDが2種類以上ある。
- その167セルのうち**165セルは移植先座標が一意で有限**、2セルは移植先座標が欠損または曖昧。後者はThompson et al. 1991のラベルに属する。
- 5研究は元の植物分類群ラベルが本プログラムの厳格な二名法検査を満たさず、「解析可能な植物種が生物学的に存在しない」とは言わない（Eranen & Kozlov 2009, Gellie 2016, Heimes 2016, Raabová 2008, Roy 1998）。事後のゆるい学名修正によってこのファミリーを救済しない。
- **原産集団そのもののlat_source/long_source列はHargreaves元CSVに存在しない**。166でも167でもなく、**新たに地理参照まで検証できた原産集団はゼロ**。観測したのは原産地IDの数だけであり、独立した座標の数ではない。
- **独立した原著実験の同一性の認証も未実施（認証0）**。研究ラベル非一致は「独立実験」の証明にならない。

本ファイルの137/1511/4947/6536といった前回の研究×種×原産地／移植先共有数は変更していない。既存のBontrager事前固定44/50種の停止も変更していない。

## 原著資料の回復順序（応答を見る前の資料取得優先度）

出典は[固定済み35研究の集計](./source_registry/hargreaves35_study_geometry_v01.json)と、[論文DOI照合台帳](./hargreaves35_primary_article_recovery_shortlist_v01.json)を参照。以下の「セル数」は元データ中の**重複する比較配置候補**であり、独立した論文数・種数・有効反復数ではない。

| 原著ラベル | 候補セル | 出典DOI | 原産地座標 |
|---|---:|---|---|
| Raabová et al., Biological Conservation 2007 | 24 | https://doi.org/10.1016/j.biocon.2007.07.007 | 原著要検証 |
| Ellis et al., New Phytologist 2007 | 18 | https://doi.org/10.1111/j.1469-8137.2007.02043.x | 原著要検証 |
| Tomiolo et al., Ecology 2015 | 12 | https://doi.org/10.1890/14-1445.1 | 原著要検証 |
| Melo et al., New Phytologist 2014 | 10 | https://doi.org/10.1111/nph.12779 | 原著要検証 |
| Volis et al., Oecologia 2002 | 10 | https://doi.org/10.1007/s00442-002-0999-0 | 原著要検証 |
| Baack & Stanton, Evolution 2005 | 8 | https://doi.org/10.1111/j.0014-3820.2005.tb01063.x | 原著要検証 |
| Brancalion et al., Ecology and Evolution 2018 | 8 | https://doi.org/10.1002/ece3.4248 | 原著要検証 |
| Arany et al., Plant Ecology 2009 | 6 | https://doi.org/10.1007/s11258-008-9530-y | 原著要検証 |
| Kindell et al., Journal of Ecology 1996 | 6 | https://doi.org/10.2307/2261336 | 原著要検証 |
| Kalske et al., Journal of Ecology 2012 | 4 | https://doi.org/10.1111/j.1365-2745.2012.02008.x | 公開HTMLのDMS表記を原本で再検証 |

上記10研究ラベルは**106/167セル**を含む（比較配置の集中による便宜的な原著資料回収の順序であり、fitness結果で選んでいない）。各書誌DOIは掲載誌・著者・年・論文タイトルの対応で確認した候補だが、**Hargreaves原本に収録された実験・集団と論文付録の実験・集団が1:1に対応することまでは認証していない**。Tomioloの論文には4植物種が登場する一方、Hargreavesの本件原本で厳格な二名法に合致するものは3種類であり、元論文の全対象をそのまま換算していない。Kalskeの公開HTMLにある集団座標のDMS表記は非標準値を含むため、文字列の見た目から座標を創作・補正しない。

## 同一性・観測・科学的予測を別々に扱う

1. **出版物の同定**：DOI・著者・年・タイトルの一致。これだけで元の庭IDや原産集団IDを解釈しない。
2. **原産地の地理証拠**：原論文のTable/Methods/Supplementにある原産集団に対応した原著WGS84座標、または原産地位置を独立に検証できる原著地名・地図。地名のみの場合は測地誤差を明示し、気候ラスタへの位置確定を自動許可しない。
3. **実験の独立性**：論文番号だけでなく集団・庭・実験年・実験IDを対応付け、同一実験の複数の報告を二重計上しない。
4. **新しい未来の解析**：適応度を開く**前**に別ファミリーとして計画・閾値・学習評価分割・帰無と検出力・地理不確実性を凍結する。本文の35、30、167、106を前回停止済みの44種へ足して50種を超えたことにしない。

## 再現可能性

- 原本：Hargreaves Dryad DOI https://doi.org/10.5061/dryad.0vt4b8gv1、Zenodo正規保存 https://zenodo.org/records/4409291。元CSVサイズ4,382,308バイト、SHA256 `f00eec59e5388da5a2243d391c4e639dc77eda8d98abbf9a8a6846b14da7bb71`。
- 固定ルール：[hargreaves35_source_geometry_contract_v01.json](./hargreaves35_source_geometry_contract_v01.json)。
- 新規スクリプト：[hargreaves35_source_geometry_v01.py](./hargreaves35_source_geometry_v01.py) は原本ハッシュと既存35研究リスト受領書ハッシュの両方を検証し、10個の許可された元CSVメタデータ列しかインデックスしない。
- 元の研究リストの受領書SHA256 `8268387cdaec3a835086418f7ce91f91777cbbc1026e148c6e6e0399647106c3`。
- 実原本の唯一の新規実行：[GitHub Actions #43708203840](https://github.com/zuizui0223/TTF/actions/runs/43708203840) 成功、既知真相の合成テスト7/7成功。
- [保存したsource-only集計](./source_registry/hargreaves35_study_geometry_v01.json) のSHA256 `7354b7fc559a16b0a2656ef251da84add8d4062d74d8e797092a2d3c76450dca`。35研究名と研究別の集計数は掲載するが、原産地／庭の生IDや座標、個体のfitness値は掲載しない。

**最終判定：`SOURCE_ONLY_SOURCE_GARDEN_GEOMETRY_AVAILABLE_ORIGIN_GEOGRAPHY_NOT_QUALIFIED`**。前回のPR #154（44<50）の停止、56植物PROCEEDのnullゲート停止、既存のHe/Fst/COI/Jensen/TTF-Qはすべてそのまま。PR #155はDraft、適応・可塑性・種分化について新たな経験的効果は報告しない。
