# TTF：Raabová原著の年次・集団同定誤りを訂正し、正しい集団名で地理品質を調べ直した結果

**2026-11-15 | Draft PR #155 | source-only / exploratory QA | 適応度・遺伝・LGM効果の応答未開封**

## 結論と訂正

**旧「Sudová 2010 → Raabová 2007の4地点」比較は対象原著の引用先を取り違えていた。** Sudová (2010) Table 1の脚注にある「ref. 50」は **Raabová et al. 2008, Oecologia, DOI 10.1007/s00442-008-1156-1** を指し、**Raabová et al. 2007, Biological Conservation, DOI 10.1016/j.biocon.2007.07.007** ではない。2007論文は2倍体6集団の2地域（K1–K3, S1–S3）、2008論文は倍数性接触帯の2倍体3・6倍体6の計9集団。Sudová表の2倍体2地点・6倍体2地点を2007原著の「3,1,5,6」に対応させた旧監査は科学的な外部検証を構成しない。

**旧 1/4・距離1.93168/0.9278/3.76521/6.51907 km は、異なる対象集団を比較した数値の履歴として残すが、原産地位置精度が3/4地点で失敗したという解釈は撤回。** データ・距離閾値を後付けで変えない。原本JSONを改変せず、[訂正JSON](./raabova_2007_2008_wrong_primary_reference_erratum_v01.json)と旧[引き継ぎ文書](./RECIPROCAL_SOURCE_AS_GARDEN_EXTERNAL_GEO_AUDIT_20261115.md)冒頭の明示的撤回で記録した。旧外部照合ワークフローの自動実行を停止し、5件の科学的参照先回帰テストに成功（[Actions #43711184196](https://github.com/zuizui0223/TTF/actions/runs/43711184196)）。

確認できる原著：
- [Sudová et al. 2010](https://doi.org/10.3732/ajb.1000114)；出版社の引用文献番号**50**を確認。
- [Raabová et al. 2008](https://doi.org/10.1007/s00442-008-1156-1)；2008年の異なる倍数性接触帯パネル。
- [Raabová et al. 2007](https://doi.org/10.1016/j.biocon.2007.07.007)；2倍体のみのK1–K3/S1–S3。
- [Castro et al. 2011](https://doi.org/10.1007/s10682-010-9439-5)；チェコ北部のマルチ倍数性集団の原著地名・GPS表。

## 原本データの集団IDによる確認

Hargreaves 2020の公表元CSV（Dryad DOI 10.5061/dryad.0vt4b8gv1、原本SHA256 `f00eec59e5388da5a2243d391c4e639dc77eda8d98abbf9a8a6846b14da7bb71`）について、**原著研究名・分類群名・source/siteラベルだけ**を認証済み原本から読み取った。

| 元研究 | Hargreaves行数 | 異なる原産集団ID | 異なる庭ID | 原著コード |
|---|---:|---:|---:|---|
| Raabová 2007 | **180** | **6** | **6** | K1 K2 K3 S1 S2 S3、それぞれ元ID側・庭側に1個 |
| Raabová 2008 | **405** | **9** | **9** | 2007のK/Sコードは現れず、分類群表記に倍数性等の追加テキストあり |

元集団・庭の生ID、原座標、適応度、fitness、遺伝結果は出力・評価していない。出力は集団コードの分類集計のみ。[元原著ID形式の監査実行 #43711576389](https://github.com/zuizui0223/TTF/actions/runs/43711576389)は成功し、合成6テストに合格。受領書SHA256 `aa4b4482c327a7a5f9e777f1e3951471272ec432611239fb67bb66b2082d52cb`。

## 別の正しい参照に基づく探索的地理照合

「旧検定が失敗したので参照先を差し替えて成功させる」のではない。**新しい探索的な原著地名・コード照合**として、Raabová **2007年**原著で正しく名前が確認できるS1=Malíč、S2=Holý vrch、S3=Encovanyを、Castro et al. (2011) Table 1の同名**2倍体**の公表DMS座標（同表でのID 1/2/3）と対応付けた。Hargreaves側の原産集団IDはK/Sコードを字句上完全一致させ、同じ研究・種・IDが移植先として現れる場合の代理座標を使用。**空間的近接性からIDを割り当てていない。**

| Raabová 2007集団コード | 同名集団の2011年原著 | 座標記載差（km） | 旧1 km距離を参考とする記述 |
|---|---|---:|---|
| S1 | Malíč | **0.735240** | 1 km以内 |
| S2 | Holý vrch | **1.306310** | 1 km超 |
| S3 | Encovany | **1.918088** | 1 km超 |

3地点のうち1地点だけが従来の1 km以内。ただし**地名が同じでも実験区画・記録位置・採集範囲が違う可能性**があり、この探索的距離は「実際のGPS精度検証」「同一実験の認証」でも「LGM適応の帰無結果」でもない。新規の事前記述済みの探索的契約：[raabova2007_castro2011_named_sites_exploratory_contract_v01.json](./raabova2007_castro2011_named_sites_exploratory_contract_v01.json)。[Actions #43711907352](https://github.com/zuizui0223/TTF/actions/runs/43711907352)で原本と名寄せコードを固定して実行・合成7テスト成功。新規結果SHA256 `f1364cf928d00c36551aa2542e4fac5be71e1fd6ac739955e5fb687d48b23b59`。**地理的閾値は緩和していない**。

## 生物学的な現在地

- 肯定的：HargreavesとBontragerの原本は同一収集系統内212原産集団で非常に近い座標を保持。相互移植の原産集団IDを移植先として持つと、35研究中21研究について76個の**一意の代理座標**・142比較候補セルを回復可能。数値は元の凍結結果から変更なし。
- まだ不十分：外部原著から集団IDに紐づく**独立測地の正確な原産地緯度・経度**を確認した状態ではない。S1–S3の原著名対応は確認したが、異なる年の採集地・実験庭の微小位置は一致するとは限らない。
- 未実施／不可：新しいLGM/0BP差の計算、元のBontrager 44/50種停止ゲートの救済、適応度・可塑性・種分化の効果量解析や新規の確認的な解釈。

**最終診断：`WRONG_PAPER_REFERENCE_CORRECTED_2007_SITE_NAMES_EXPLORATORY_ONLY_SOURCE_GEO_UNCERTAIN`**。元のPR #154/旧TTF-Q・He/Fst/COI/Jensen/PROCEED凍結結果と`main`は変更しない。PR #155はDraftのまま。


## SHA固定後の永続保存と実行状況

- Raabová 2007/2008元データのK/S集団コード分類を[研究ID別の恒久JSON](./source_registry/raabova_2007_2008_original_ID_metadata_v01.json)に保存（正規の元原本監査 [Actions #43711576389](https://github.com/zuizui0223/TTF/actions/runs/43711576389)、SHA256 `aa4b4482c327a7a5f9e777f1e3951471272ec432611239fb67bb66b2082d52cb`）。
- 新しい探索的なS1–S3の公表地名座標比較を[同名集団3地点の恒久JSON](./source_registry/raabova2007_castro2011_3site_exploratory_v01.json)に保存（[Actions #43711907352](https://github.com/zuizui0223/TTF/actions/runs/43711907352)、SHA256 `f1364cf928d00c36551aa2542e4fac5be71e1fd6ac739955e5fb687d48b23b59`）。
- 2つともGitHub Actions [#43712276329](https://github.com/zuizui0223/TTF/actions/runs/43712276329) で元の成果物ZIPからハッシュ再照合し、PR #155の開発ブランチへ恒久化した。元CSV全体、集団の生ID／座標、fitness値は保存していない。
- **旧Sudová×Raabová2007の1/4は科学的に不適切な原著取り違えであり、現在のQAでは使用禁止。** 2011年公表地名へのS1–S3探索的比較は正しい研究名に対応するが、現地GPS測地精度の証明にはならない。
