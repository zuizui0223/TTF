# Baack & Stanton 2005：3回に層化された原著移植実験のsource-only監査

**2026-11-19 | Draft PR #155 | 生存・生長・遺伝・LGM応答の新規モデルは未開封**

原著：Baack & Stanton (2005), *Evolution* 59:1936–1944, DOI https://doi.org/10.1111/j.0014-3820.2005.tb01063.x . PDFのMaterials and MethodsとTable 1に原著の原産集団ID・移植先ID・時期（Early/Middle/Late）が明記されている。地図から座標を視覚推測しない。

## 実際の原著の単位

全研究ではHoosier Pass近傍の**10集団（2倍体D=4、4倍体T=6）**を対象とするが、用途ごとのsource IDとgarden IDが異なる。

種子移植実験は3個の2倍体原産集団と5個の4倍体原産集団の種子を2移植先へ配した設計で、原著Seedling Transplant ExperimentのIDセットではない。**この2実験を同じsource×garden行列とみなさない。**

苗移植実験（2000年）の原著Table 1：

| 時期 | 原産4集団 | 移植先4集団 | 自集団戻し |
|---|---|---|---:|
| Early | D-LP, D-HPNG, T-HPWM, T-PM | D-LP, D-HPNG, T-HPWM, T-PM | 4 |
| Middle | D-LP, D-HPNG, T-HPWM, T-WP | D-HPNG, D-HPWW, T-P7, T-WP | 2 |
| Late | D-BP, D-HPNG, T-PM, T-WP | D-BP, D-HPNG, T-PM, T-WP | 4 |

集計：**6ユニーク原産集団、8ユニーク移植先、3時期×4×4＝48の時期固有移植アーム、時期を無視すれば38ユニーク有向source→garden組、local比較10アーム、異倍数性の優占地への移植24アーム**。

- 原著は3時期を**別の実験として解析**する理由として、融雪時期と利用できる原産地・庭の違いを挙げる。全時期に均一な6×8完全交差が行われたわけではない。
- 中期の移植先T-P7、D-HPWWは原産地として苗を提供していない。したがって全sourceと全gardenを単純に掛け合わせ、存在しない移植アームを補完することは禁止。
- Hargreaves候補短冊の**6原産集団ラベル、8庭ラベル、8比較セル**は元研究の苗移植ID数と一致する部分があるが、**8二次セルと原著48時期固有アームは異なる集計単位**。差の40を「欠落アーム」とは扱えない。
- 異倍数性優占地への移植は送粉による遺伝子交流や系統分化の独立測定ではない。倍数性が強く地理的に偏っており、これを気候・LGM効果と読み替えない。

原著図には集団位置図があるが、GPS精度が保証された原産集団座標を回収していない。隣接するBaack 2004の35集団分布記録・地図を、同名集団IDと実験コホートの直接照合なしに2005年実験原産集団の真正位置へ代用しない。

## 判定

`baack2005_round_stratified_original_design_v01.py`と9テストで元Table 1の48/38/10/24、6 source、8 garden、時期の差をfail-closedで復元する。

**ROUND_SPECIFIC_EXPERIMENT_DESIGN_NOT_COMPLETE_6x8_MATRIX_OR_LGM_PANEL**。元論文の倍数性・種分化の結果を本研究で再検証したものではなく、historical fitness panel資格は引き続きHOLD。Bontrager 44/50等の凍結した科学的判定は不変。PR #155 Draft、main未変更。
