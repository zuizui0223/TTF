# Ellis 2007：移植比較グラフの識別可能性（source-only）

2026-11-17 | TTF Draft PR #155 | 記述的追加監査。fitness、AFLP、LGM適応度、可塑性、種分化の応答は一切開いていない。

## 結果

凍結済みの原著照合受領書（SHA256: 28cea27d59affd846c9ee490f20daedd36220e6e75041f429663a208dcd239ce）だけを入力とした。3流域（Ariz/Sout/QK）×2石英被覆（dense/sparse）の6集団をノードとし、原著・二次抽出に収録された source→garden の移植配置をエッジとする。

| 事前に報告された比較配置 | 原著の実験設計 | Hargreavesの二次抽出 |
|---|---:|---:|
| 有効なsource→garden比較 | 24 | 18 |
| 同流域・異石英被覆 | 6 | 0 |
| 比較グラフの無向連結成分 | 1 | 2 |
| 成分サイズ | 6 | 3と3 |
| 配置として存在するアーム区分 | 3 | 2 |

**新しい数学的診断：二次抽出は石英被覆の種類ごとに2つの連結成分へ分断されている。** 原著にあった6本の同流域・異石英被覆移植アームが全て抽出されていないため、二次データだけではその実験的な移植応答の比較を直接識別できない。局所移植と異流域・同石英被覆の比較配置は残るが、その差が地理障壁の因果効果であると認定したわけではない。

このグラフのエッジは**移植比較**であり、野生集団間の遺伝子流動でも第四紀の移動回廊でもない。二次抽出に比較が無いことは原著の実験が存在しなかったことを意味しない。異流域・異石英被覆の12組は原著設計で実施されていないため欠損に数えない。

## 再現性

- 入力：source_registry/ellis2007_24_to_18_origin_garden_arm_audit_v02.json（上記SHA固定）
- 実装：ellis2007_transfer_graph_scope_v01.py（集計から一意に復元できる完全区分のみグラフ化、部分区分ならfail closed）
- テスト：test_ellis2007_transfer_graph_scope_v01.py（元受領書SHA照合、観測成分、6本欠落、重複区分・破損区分・応答境界の否認）
- 出力：source_registry/ellis2007_transfer_graph_scope_v01.json（個体ID・生座標・fitness数値無し）

実行コマンド：

    python -m unittest discover -s research/prospective -p test_ellis2007_transfer_graph_scope_v01.py -v
    python research/prospective/ellis2007_transfer_graph_scope_v01.py --output _source_only/ellis2007_transfer_graph_scope_v01.json

2026-11-17のオフラインの同等ロジックで7テストが合格。GitHub側での再実行は別途Actions結果を確認すること。

## TTFとしての判定

COMPARISON_GRAPH_HAS_TWO_QUARTZ_COMPONENTS: **利用可能データの識別可能性の限界**が新たに明確化された。LGMの独自予測力・適応・可塑性・種分化の肯定結果ではない。

Bontrager 44/50種の停止、PROCEED、北米哺乳類80種の凍結ΔAUC、211種の遺伝TTF、mainの科学的結論を更新しない。Ellisの単一種を旧停止ゲートに追加しない。
