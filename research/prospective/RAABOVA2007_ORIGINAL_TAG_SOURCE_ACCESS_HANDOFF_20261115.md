# TTF｜Raabová 2007年成体移植の個体札番号表を探す source-only 実行記録

**2026-11-15 | Draft PR #155 | 2007原本取得ゲート | 適応度・遺伝・LGM効果は未開封**

## 現在の結論

**同じ2004年開始の成体相互移植コホートを2論文が報告した事実は著者自身の2014年論文で確認済み。しかし、2007年論文の元個体札番号表は現在の公開アクセス経路では取得できなかった。**

これは「元表が存在しない」という意味ではない。次の出典へ実際にアクセスした結果を、取得制限とデータ不在を区別して記録した。

| 一次資料・リポジトリ | 入口 | 今回のアクセス結果 | 個体タグの元行を取得 |
|---|---|---|---|
| Raabová et al. 2007 出版社 | [Biological Conservation DOI](https://doi.org/10.1016/j.biocon.2007.07.007) | **HTTP 403** | いいえ |
| Raabová 2008 博士論文 | [カレル大学公式記録](https://dspace.cuni.cz/handle/20.500.11956/94265) | **TimeoutError** | いいえ |
| Raabová 2009 同テーマ学位記録 | [カレル大学公式記録](https://dspace.cuni.cz/handle/20.500.11956/95341) | **TimeoutError** | いいえ |

博士論文リポジトリには本文や要旨のPDFが掲載されていること自体はWeb検索の公式記録から確認できるが、**元CSV／Excelの不在を全アーカイブについて確定することはできない**。本文PDFは個体別データ表の代替にならない。記事の本文だけから個体番号を推定したり、他の年のIDを補完しない。

出典探索の事前契約：[`raabova2007_original_data_location_preflight_contract_v01.json`](./raabova2007_original_data_location_preflight_contract_v01.json)。コード：[原データ出典索引の監査](./raabova2007_original_data_location_preflight_v01.py)。正規実行：[GitHub Actions #43723238370](https://github.com/zuizui0223/TTF/actions/runs/43723238370)、**6件の合成アクセス安全性テスト成功**。機械可読の状態：[2007元個体ID公開ソースアクセス受領書](./raabova2007_individual_tag_source_access_receipt_v01.json)。

## 既に取得できている2014年側の識別情報

Pánková et al. 2014 [PLOS ONE Table S1](https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0093967.s001&type=supplementary) は**原本XLS**を取得・ハッシュ照合済み。原本の説明シートに `mark` = "The number of the plant in the field"、`St` = "Target site"、`So` = "Site of origin" とある。適応度・生存率・葉数・バイオマスなどの数値値は分析していない。

- 2014年`survival`：547個の一意な`St×mark`。
- 2014年`performance`：281個の一意な`St×mark`。
- 同じ2014年の2シート間：280キーが一致し、**279キーは`So`も整合**。1キーは原産集団ラベルに矛盾、1個は`performance`だけに存在。勝手に修復しない。
- **2007対2014の個体タグ番号の一致数 = 未評価（null）**。同じ物理実験を追跡したことは、個体番号が全年度で変わっていないことの証明とは異なる。

[2014個体マーク内部監査結果](./source_registry/pankova2014_original_plant_mark_internal_linkage_v01.json) と [2014原本説明文書](./PANKOVA2014_ORIGINAL_PLANT_MARK_CROSSSHEET_HANDOFF_20261115.md)が確定情報。

## 2007年の未解決入力に必要なもの

個体重複を新たに評価するなら**2004年に設置された成体移植コホートの2005/2006年時点の元個体番号表**か、2014年の番号札との**対応表**が必要である。照合に最小限必要なのは元の`mark`または個体ID、移植先K1–K3・S1–S3、原産集団K1–K3・S1–S3、2004年成体移植であることの実験識別、観測年であり、適応度・葉長・開花数などの生物応答値は不要。

この原本が公開された場合に限り、以下を**応答未開封の別契約**で評価する：

1. 書誌／実験コホート同一性確認（既に完了）。
2. 2007年原本個体番号の空白・重複・同一IDに複数の原産地があるかの品質監査。
3. 2014年`St×mark×So`と年別IDのマッピング、1対1・多対1・矛盾の集計。
4. 原本間で同一タグ体系であるという証拠がなければ**「類似した番号」を個体同一と自動認定しない**。

**第三者著者へのデータ提供依頼は送っていない。** 研究者本人の承認・正規の公開経路なしに元データの個体列を取得できたと主張しない。

## 科学的停止境界

成体移植コホート重複排除は**2論文・1物理実験**で確定。2007年に別途実施された種子播種2003年と2004年は、この成体移植とは別の実験として扱う。2014年の279件の内部照合を**2007年との一致個体数に変換しない**。

Bontrager PR #154の**44/50種停止判定**、PROCEED、TTF-Q、既存のHe/Fst/COI系の凍結解析は維持。追加の適応度係数・可塑性・種分化・LGMによる進化応答は未推定。PR #155はDraft、`main`は未変更。

**終端判定：`ORIGINAL_2007_TAG_ROW_DATA_NOT_RETRIEVED_PUBLISHER_403_THESIS_TIMEOUT_CROSSPAPER_MATCH_UNKNOWN`**。
