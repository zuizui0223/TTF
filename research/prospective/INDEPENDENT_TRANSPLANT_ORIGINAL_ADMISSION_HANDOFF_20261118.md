# TTF｜独立相互移植研究の原著資格審査：source-only 第1段

2026-11-18 | Draft PR #155 | 生物学的応答を開かない

## 現在地

Ellis 2007の実験配置は、元24アームと二次抽出18アームの差、設計行列の階数7→6、零空間、64通りの復元感度、63通りの等分散仮定下の精度感度まで検証した。これらは設計監査であり、LGMによる適応を示すものではない。さらなる配置解析より、独立した一次データの資格認定を優先する。

## 独立パネルの資格（新規データを開く前に要求）

1. 原著に対応する原産集団IDと、その集団の独立に検証可能な位置・測地精度。
2. 原著の実験コホート、移植先、年、個体・集団単位の重複照合。論文やメタ解析の違いは独立実験の証明ではない。
3. 原産地×移植先の比較が存在し、必要な環境対比を設計行列として識別可能。
4. 現在気候・距離・障壁を調整しても独立した歴史気候変位が残るか、応答を見る前に評価。
5. 独立種・独立研究での学習/検証分割、偽陽性と検出力、停止基準を別の新規契約として固定。

## 今回実装した source-only 第1段

凍結済みHargreaves原著回収候補10件（比較配置候補106セル）を機械的に再点検する。これら106セルは独立反復ではない。元の登録では原産地GPSが認証済みの研究0、独立した原著実験の認証0、LGM履歴対比の資格0。**これは実験や座標が存在しないという結論ではなく、必要な証拠がまだ揃っていないという意味。**

原著メタデータの次の確認対象として、Brancalion 2018（森林3地点の原産集団位置）、Baack & Stanton 2005（倍数性と原産集団ID）、Kalske 2012（非標準DMS座標の原本確認）を挙げる。Ellis 2007は元の欠落6アームと原産地GPSが必要。Raabová 2007 / Pánková 2014は同一成体コホートの重複計数を禁止。Tomiolo 2015は元6 source IDを6個の独立GPSと数えず、代表2地点のままHOLDする。これらは元の資料構造による資料回収優先度であり、適応度効果に基づくランキングではない。

`independent_original_transplant_admission_v01.py` は10研究の出典・研究ID・候補セル数・資格状態を再集計し、独立実験やGPSを推測してPASSにしない。7つのfail-closedテストをActionsで実行し、結果をsource-only artifactに出力する。応答列、個体fitness、AFLP、LGM気候値は一切参照しない。

## 科学的判定

**NO_NEW_INDEPENDENT_HISTORICAL_FITNESS_PANEL_QUALIFIED**。Bontrager 44/50停止を事後救済しない。既存TTF-Q、GenDivRange He、CaliPopGen Fst、Jensen break-location、211種TTFの凍結結論も変更しない。PR #155はDraft、mainは変更しない。


## 2026-11-19 source-only追記（新しい独立パネルの合格ではない）

- **Kalske 2012**：原著の4原産集団DMSは、厳密構文で3/4不正。構文のみ正常なUtöも通常変換で66.498333°N、15.357778°Eとなり、論文が指定する南西フィンランド群島と地理的に矛盾する。4/4原産地座標をQUARANTINE。自動繰り上げ・他論文の同名集団GPSの流用を禁止する。コード `kalske2012_dms_geographic_quarantine_v02.py`、記録 `KALSKE2012_DMS_GEOGRAPHIC_QUARANTINE_V02_20261119.md`。
- **Brancalion 2018**：原著Methodsによれば3森林原産集団×各10母樹＝30母樹、同じ採種ロットを共通圃場・種子移植・苗移植の3方式に使う。3方式を3独立母樹コホートとして計数しない。Hargreaves由来の8比較セルも独立試験8件ではない。コード `brancalion2018_maternal_cohort_gate_v01.py`、記録 `BRANCALION2018_SHARED_MATERNAL_COHORT_HANDOFF_20261119.md`。

対象候補10件の合否状態は変わらず、**独立原産地GPS認証0、独立実験認証0、歴史気候対比資格0**。この0は無存在の証明ではなく資格認定の保留。凍結済みLGM・適応・TTFの科学的応答は未開封・未変更。


## 2026-11-19 原著source-only追加監査

**Brancalion 2018**：原著Table S1取得をGitHub Actions上の公式WileyとPMCの4経路で実行したが、403/404/404/404で真正DOCXの内容は未取得。3原産集団・30母樹、24ブロック、72プロット、各苗／種子移植計2,160配置の実験設計を確認。植物の配置数は原産地の独立気候曝露数を増やさず、**3原産地だけでは切片＋現代気候2項＋LGM1項の4列モデルは必ずランク不足**。詳細は`BRANCALION2018_ORIGINAL_TABLE_S1_SOURCE_RANK_STOP_20261119.md`、テストとActionsはPASS。単独の新規LGM効果パネルへは昇格しない。

**Baack & Stanton 2005**：原著Table 1（2000年の苗移植）は3時期それぞれ4原産×4庭の異なる組合せを使用。**48 round-specific arms、時期を畳むと38 unique directed source→garden pairs、6原産・8庭、10 local comparisons**。Hargreavesの8候補セルとは集計粒度が違うため、差分を「欠損40」と呼ばない。元の原産集団GPSは依然未認証。倍数性（2倍体/4倍体）、融雪季節、地理的集団の混同を防ぐ設計判定を`BAACK2005_ORIGINAL_ROUND_STRATA_HANDOFF_20261119.md`に記録した。

**資格判定を維持**：真正の独立原産地GPSをもつ候補研究0、独立原著実験の資格認定0、historical contrast資格0。これは実験が無いことや結果が否定されたことではない。旧Bontrager 44/50と既存科学的結果の再調整・事後救済はしない。


## 2026-11-19 Melo 2014 original-source geography gate

**Melo et al. 2014** は、元論文 `The system` のCabarita Beachの2原産集団（Dune/Headland）のDMS記載が構文上有効である一方、両点の直線距離**3.505057 km**は原著の定性的説明**約0.8 km**と整合しない（差2.705 km、約4.38倍）。元位置点が採集時GPSか代表点か未認証で、出版誤記の原因も同定していない。原著DMSを上書き・推定修復せず **SOURCE-GEO HOLD**。Hargreavesの原産2ラベル・庭2ラベル・10候補セルをそのまま認証済み独立研究数としない。詳細は `MELO2014_ORIGINAL_COORDINATE_DISTANCE_CONFLICT_20261119.md`、再現コード `melo2014_original_coordinate_distance_gate_v01.py` と11件の専用テストを追加。元独立GPS認証0、LGM応答・新規遺伝・適応データ未開封。

### source-only資格プログラムの現在判定

KalskeはDMS不正3地点＋地理不一致1地点、Brancalionは3原産地点のsource-level design ceiling、Baackは3時期で異なるsource×garden比較、Meloは公開座標と本文距離の不一致。それぞれ**異なる資料品質・識別の失敗モード**であり、同じ研究を形式上の研究数や座標数で水増ししない。10原著候補のうち一次位置・実験ID・歴史気候対比まで認証された新規独立検証パネルは引き続き0。

**後続の主目的は、原産地GPSと実験コホートの真正記録が回復できる研究を得ること**。既存のBontrager 44/50停止を後付け再解析で救済しない。元の表現型・遺伝・局所適応・可塑性・種分化・LGM効果は新たに推定していない。


## 2026-11-19 Arany 2009 source/garden実験単位の追加検証

**Mosleh Arany et al. (2009), *Plant Ecology*, DOI 10.1007/s11258-008-9530-y** のオープンアクセス原著Methodsを照合。原産5集団（Meijendel砂丘3、Leiden/Noordwijk内陸2）、各集団10母植物、実験庭2（砂丘、Leidenの代替庭）、母植物×庭ごとの苗8、総配置**800苗**。原著は「元の内陸採種地点にはenclosureを設置できなかった」と明記し、Leidenの庭を実験代替として使った。砂丘原産2集団は**5m**間隔。したがってHargreavesの5原産ラベル、2庭ラベル、6比較セルは独立した原産地気候5点のGPSが認証されたことを意味せず、gardenをsource GPSに代用しない。

実装：`arany2009_original_source_garden_grain_v01.py`、9テスト、`ttf-arany2009-original-grain-v01.yml`。原著設計の内容を `ARANY2009_ORIGINAL_SOURCE_GARDEN_GRAIN_20261119.md` に固定した。**新たに原産地GPS認証0・独立歴史気候対比資格0を維持**。同研究の既刊生態学的適応結果を否定しない。

既に検証済みMeloの原著座標と著者関連博士論文の地点照合を14テスト版へ修正し、GitHub Actions Push run [#44482122147](https://github.com/zuizui0223/TTF/actions/runs/44482122147)でPASS。関連研究の同じ地点座標を新しい独立個体コホートと誤計数しない。
