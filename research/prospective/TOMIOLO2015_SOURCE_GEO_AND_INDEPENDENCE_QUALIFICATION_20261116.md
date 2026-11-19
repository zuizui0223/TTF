# TTF｜未照合35研究の独立性資格認定：Tomiolo 2015の原著Figshare・原産地位置検証

**2026-11-16 | Draft PR #155 | source-only pilot | 生物学的適応度・21ka LGM・遺伝・可塑性・種分化の応答未開封**

## 現在の結論

**Tomiolo et al. (2015) のHargreaves由来「原産集団ID 6個」は、一意な位置座標としては2地点である。** 2015年論文の公式Figshare原著付録Appendix AにあるM（地中海性）とSA（半乾燥）の地点座標と、同一研究・種・原産集団IDが移植先にも現れる場合だけに限定したHargreaves座標は一致（M:表示精度内で0km、SA:0.0004km）。しかし**別の独立GPS実測という証拠ではない**。両者が公表された分単位座標を共用・転記した可能性があり、位置の微小誤差は未評価。3つ目のA（乾燥）地点はAppendix Aの気候傾度表に登場するが、**「原産集団が3つに増えた」という根拠にはならない**。

## 原著・原本の照合

- 論文：[Tomiolo, Van der Putten & Tielbörger 2015, Ecology](https://doi.org/10.1890/14-1445.1)（原産地×土壌×近隣植物を組み合わせた相互移植）。
- 論文にリンクされた元資料：[Figshare collection 10.6084/m9.figshare.c.3307653.v1](https://doi.org/10.6084/m9.figshare.c.3307653.v1)。
- Figshare metadata-onlyの原本一覧：Appendix A `appendix-A.htm` 134,278バイト（MD5 `f2dc5a3c1715a72d297112434dac8f45`）、Appendix B `appendix-B.htm` 7,084バイト、Appendix C `appendix-C.htm` 104,739バイト。ファイル一覧の照合には[Actions #43882113357](https://github.com/zuizui0223/TTF/actions/runs/43882113357)、6合成テスト成功。
- Appendix Aの研究地点・実験配置説明だけを元HTMLから取得。本文ファイルSHA256 `3d9f27e87755c0d0c28f2e7d4798e73df20ad8bb70dc9236fb0302c7e4348447`。[Actions #43882985831](https://github.com/zuizui0223/TTF/actions/runs/43882985831)、7合成テスト成功。数値の生存率・結実・花数・繁殖等の応答値はモデル化していない。初回の限定的サイト文脈抽出に定性的な「mortality」記述が紛れ込んだため、後続実行では**出力文脈を遮断する技術的な除外語**を追加した。これを「完全な結果未閲覧の原著論文本文」とは誤称しない。
- 著者自身の2013年[テュービンゲン大学博士論文](https://publikationen.uni-tuebingen.de/xmlui/bitstream/handle/10900/49977/pdf/Sara_TomioloTHESIS_040913.pdf)では、地中海性サイトMがN31°42′ E35°03′、半乾燥サイトSAがN31°23′ E34°54′と記載され、両方とも**分単位の代表地点**。Appendix Aの表と一致する。学位論文と2015年原著を独立した移植実験2件として加算しない。実験IDの直接照合は未認証。
- 原本Hargreaves Dryad `10.5061/dryad.0vt4b8gv1` はSHA256 `f00eec59e5388da5a2243d391c4e639dc77eda8d98abbf9a8a6846b14da7bb71` で検証。保存済み資格資料の[Tomiolo論文の原産地ID/庭ID幾何](./source_registry/hargreaves_source_as_garden_georecovery_v01.json)は**6原産集団ID・2個の独立した位置点・同一庭比較12セル**。4種が原著で扱われた一方、現行のHargreaves厳格二名法対象は3種。学名を結果後に再採点して新しい種を足さない。

### 1. 分単位の元サイト位置との距離

既存のHargreaves原著ファイルにある`study×binomial taxon×source`と、`study×taxon×site`が文字列完全一致する場合にのみ移植先WGS84から「原産地代理座標」を取った。二地点を匿名P1/P2として緯度順に固定し、Figshare Appendix Aに公表されたM・SA・Aの三地点との**全2×3距離行列**を計算した。近いから元IDを書き換える処理はなし。

| 匿名の代理位置 | SA（半乾燥）との距離 | M（地中海性）との距離 | A（乾燥）との距離 |
|---|---:|---:|---:|
| P1 | **0.0004 km** | 37.9732 km | 58.8356 km |
| P2 | 37.9728 km | **0.0000 km** | 96.4949 km |

- 元データ資格結果：[Actions #43883825212](https://github.com/zuizui0223/TTF/actions/runs/43883825212)、合成7テスト成功。
- [オリジナルのGitHub Actions結果をバイト一致で保存したJSON](./source_registry/tomiolo2015_geo_proxy_matrix_v01.json) SHA256 `78f1096d0e3ba985a6bfc55f2d7f2911f8b8187322b57ef0243fa2c35bf8e98c`（git blob `140c98b13aecb296c26d30d804a8f00a0b7ce340` と元成果物が一致）。
- [3ファイルと位置情報のsource-only provenance](./source_registry/tomiolo2015_figshare_appendixA_source_provenance_v01.json) にファイルサイズ・MD5・地理点の独立性未認証を記載。

**原著表との座標の近さは、座標出典の再利用／整合性の証拠であって独立測地の裏付けにはならない。** 公表座標が分単位のため、緯度で±0.5分 ≈0.9km、経度で緯度31°付近なら±0.5分 ≈0.8km程度の表記上の丸め誤差がありうる。研究サイトそのものが100×250m規模なのも含め、**気候ラスタの1kmセルを確定的な原産地曝露として扱えない**。

## 2. Tomioloは歴史気候×局所適応の独立検証を今すぐ行えるか

**行えない。** 根拠：

1. **独立した歴史気候の入力地点数は、6 source IDではなく最大2サイト**。この違いを無視すると擬似反復になる。
2. 独立した元原産集団GPSはまだ確認できず、M/SAの代表地点だけが原著に載る。Hargreaves座標はその値と事実上一致する。
3. 2013年著者学位論文と2015年論文は同じ地域・植物・実験プログラムを扱い、別の独立実験のサンプルとして数えてはならない（個体／実験IDの直接重複数は未検証）。
4. 原著Figshare資料に独立した4種×全原産集団GPSのデータファイルは見つかっておらず、確認できたのは研究地・土壌・図表などの付録。これはファイル構成の確認であり、非公開GPSが存在しないという断定ではない。
5. LGM 21ka/現代0BPのソース別気候差は計算しておらず、適応度の正負・効果量も見ていない。旧Bontragerの**44/50**凍結ゲートにTomioloの元種・IDを後付けで足さない。

判定は`ORIGINAL_PUBLISHED_SITE_LOCATIONS_CONCORDANT_NOT_NEW_INDEPENDENT_GPS_OR_HISTORICAL_FITNESS_PANEL`。この研究は地理参照監査の事例として残し、確認的な旧パネルを救済しない。

## 3. 次の35研究資格審査で優先すること

- **Ellis et al. 2007, Argyroderma pearsonii**：[原著](https://doi.org/10.1111/j.1469-8137.2007.02043.x)に6集団、3流域×2石英被覆条件という明確な原著設計があり、Hargreavesに一意な庭由来代理位置は6点、比較候補18セル。ただし原論文のFig. 1は位置図で、元集団コード→高精度GPSの文字列テーブルは未確認。2006年のEllisら別論文にも同属の位置情報があるが、**同じ6原産集団と直接同定できるまで転記しない**。6集団が10km前後の範囲に分布するため、LGM/現在気候の独立歴史対比が残るとも仮定しない。
- **Brancalion et al. 2018, Euterpe edulis**：[原著](https://doi.org/10.1002/ece3.4248)はブラジルのCaetetus・Carlos Botelho・Ilha do Cardosoの3森林タイプで、既存代理位置は3点・比較候補8セル。ただし原著で確認できたのは主に公園・森林の位置と補遺の環境条件であり、**採種母集団に対応する精密GPS**は未認証。広大な保護区の代表位置を原産集団の真値に変換しない。

これまでのHargreavesの未照合35研究に存在する**同一移植先での複数原産集団ID比較**は167セル、相互移植先由来の座標代理を2点以上持つ研究ラベルは21。うちRaabová2007とその同一成体移植コホートの後続Pánková2014は重複計数してはならない。**残りは「研究ラベルとして候補」であって認証済みの独立原著実験ではない。**

既存PR #154 Bontrager 44/50、PROCEED、TTF-Q、He/Fst/COI/Jensenなどはすべて不変。PR #155はDraft、`main`へは統合しない。
