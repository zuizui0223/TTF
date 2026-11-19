# TTF：Volis 2002／Kindell 1996の原著スコープを閉じる

**2026-11-19 | Draft PR #155 | Source-only / 応答盲検 / 10研究shortlistの第1段資料点検**

## Volis et al. 2002 (*Oecologia*)

**一次資料**：原著 DOI https://doi.org/10.1007/s00442-002-0999-0 。著者のベングリオン大学サーバー掲載の出版PDF https://tzin.bgu.ac.il/life/Faculty/Volis/pdf/pdf4.pdf 、Methods（原著p.132–133）、Figure 1の移植設計を確認。

- 1996年に採種した原産集団は**イスラエルの地中海型M（Beit Guvrin National Park）1集団と砂漠型D（Negevの実験用wadi）1集団**。**原産集団2**、庭M/Dの**2現地環境**。原著は保護区とNegev実験区の記述を載せるが、原産地集団そのものの独立した緯度経度を当該Methodsから認証していない。
- 植生除去移植では**各環境3母植物の全同胞から混合種子プール**を作成して使用。個体母樹がプール内で独立した原産地気候を代表するわけではない。密度2段階（16/64株m⁻²）×純植／混植があり、追加の独立原産集団ではない。
- 自然植生を残す別実験では各環境**10母植物**由来の2週齢苗を相互移植。植生除去実験の3母植物混合プールと同じ統計的単位・同一苗IDとはみなさない。D/Mのそれぞれで現地サブサイト2区分を使うが、これは**2つの原産地気候単位が4つに増えたわけではない**。
- 同じ研究チームが2002年に公表した別の*Evolution*論文（DOI 10.1111/j.0014-3820.2002.tb01453.x）は地域ごと5集団×2で形質比較を行う**別論文・別研究配置**。こちらの10集団をOecologiaの元移植source=2に後付けで足さない。Hargreaves候補にはsource=2, gardenラベル2, 同一庭候補セル10。

**結果：原著2原産集団では、LGMの独立寄与を現代気候・地理・砂漠／地中海性生境差から分ける多変量モデルを単独では資格認定できない。** 新しい生物学的効果を推定したわけではない。

## Kindell et al. 1996 (*Journal of Ecology*)

資料：DOI https://doi.org/10.2307/2261336 。今回直接確認したのは**出版された抄録の再掲記録**（https://ournativeplants.com/publications/surrounding-vegetation-transplant-age-detection-local-adaptation-perennial-grass-aristida-stricta-7206/）とJSTORの掲載情報であり、**原著Methodsの全アーム表は未検証**。確定できる構造と未確認の構造を区別する。

- アメリカの多年生草本*Aristida stricta*の**2原産集団を2移植地へ相互移植**し、植生除去／未処理の**2処理**、苗齢**3段階**を比較した。
- 設計因子の単純な組合せは2 source×2 gardens×2 vegetation×3 ages=24カテゴリー候補となるが、**各カテゴリーに実際の観測が存在したと立証したわけではない**。この24を「原著24アーム観測済み」と呼ばない。
- Hargreaves二次抽出は**2原産ラベル・6 gardenラベル・6比較候補セル**。二次garden IDが原著の地理的に独立した「6庭」である根拠はなく、**元の現地移植先は2地域**。二次IDの意味は元データとの詳細joinが必要。
- 元の個体・母系ID、元原産集団GPS、実験地点の座標精度、元の全アーム表は未認証。

## 10研究shortlistの資格に関する結論

Raabová、Ellis、Tomiolo、Melo、Baack、Brancalion、Kalske、Aranyに続き、**VolisとKindellをsource-onlyレベルで調査**。10候補すべての研究単位の初回監査が完了。重要な留保として、Kindellなど一部は本文全量未認証であり、初回資料調査の終了と**元データ・研究ID・独立GPSの資格認定が完了したことは別**。

`volis_kindell_original_study_scope_closure_v01.py` は元の10研究ショートリストのVolis/Kindell行と資料から判明した設計の上限を照合し、原産地GPS認証0、独立実験コホート認証0、独立歴史気候対比資格0を維持する。10件の合成テスト・専用Actions・source-only成果物で再現する。

**初回10研究の判定：NO_NEW_INDEPENDENT_HISTORICAL_FITNESS_PANEL_QUALIFIED**。これらを旧Bontrager 44/50停止の事後救済に使わず、LGM適応や可塑性・種分化の新しい肯定効果も主張しない。今後は別の独立原著コホート・真正GPSを得て初めて新規パネルの事前資格を設計する。
