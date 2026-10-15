# TTF — Corrected NOAA ETOPO land-only sensitivity for the original 26 Costa Rica localities (2026-10-15)

**Source-pinned exploratory modern topography audit; NOT a historical paleoshoreline, not a new genetic or climate LCP hypothesis test.**

## Root cause of failed v0.1 ETOPO runs

The v0.1 contract had mistakenly selected `ETOPO_2022_v1_15s_N00W090_surface.tif`. NOAA's tile code gives the **northern tile boundary**, so `N00W090` spans **15°S–0°**, outside *Bungalotis erythus*'s 10.455°–11.05°N original localities. The 2026-10-13 runner confirmed bounds [-90,-15,-75,0] and compound CRS **EPSG:9518**, which combines horizontal WGS84 with vertical EGM2008 height; the EPSG:4326-only assertion was also false. Original CI runs #38407040635/#38407334076 therefore failed **correctly** and no fictitious mask result was accepted.

Source-only v0.2 contract was committed **before looking at the first correct-tile elevation value**:
`docs/supporting/historical_host_lcp_etopo_static_mask_contract_v0.2.json`
(SHA256 `af9a7c3aa0b6c9708d9fc26e13a1ca52ac8181f3254a4d4d92618792dd9431a6`).

The metadata-only workflow [#38673983044](https://github.com/zuizui0223/TTF/actions/runs/38673983044) fetched the correct NOAA `N15W090` tile (0–15°N, 90–75°W), EPSG:9518, 3600×3600 original native 15-arcsecond cells, 23,311,599 bytes, and exact SHA256
`2e0f33efebcbe6cc50f3a6dc7c13c78f46e04f648b5f6f3a5ca180c8a5c2f35c`.
This SHA and metadata were frozen **before** reading ROI elevation values in
`benchmarks/frozen/historical_host_lcp_etopo_N15W090_source_pin_v0.2.json`.
Source artifact ID `11944988218` includes the original public TIFF and metadata receipt.

## Real native topography result

After the source pin, [GitHub Actions #38674429445](https://github.com/zuizui0223/TTF/actions/runs/38674429445) ran the exact original 26 response-blind localities and the fixed 0.12° bounding-box margin. All **6 synthetic compound CRS/frozen threshold/source guard tests passed**.

- Native ETOPO window: **201 rows × 114 columns = 22,914 cells**.
- Valid original ETOPO elevation cells: **22,914/22,914**.
- Present land (>0m EGM2008 elevation): **22,914/22,914 (100%)**.
- Present submerged elevation cells: **0**.
- Additional cells for an *illustrative, static -120m sea-level threshold*: **0**.
- Under both present and static -120m rules: **22,914 accessible land pixels**.
- All **26/26 original localities** lie above present 0m; ETOPO sampled heights **147.6538–1008.3745m**.

The frozen output (byte-exact GitHub Actions JSON, SHA256 `9b793cf4d249f3447a66b20193a443553570fdb6d46e1086501060ba9ce2a8f3`) is
`benchmarks/exploratory/historical_host_lcp_etopo_N15W090_static_land_readiness_v0.2.json`. Runner artifact 11944919142.

**Ecological interpretation:** For this **specific** ROI, neither current ETOPO sea pixels nor a hypothetical uniform 120m sea-level fall creates/removes grid corridors. The previously reported `+1.666` host and `+17.794` self climate cost-km LGM changes cannot be attributed to presently submerged terrain in this ROI or to the tested *static* eustatic lowering on unmodified ETOPO elevations.

**Important limitations**:
1. This is NOT a 21ka coastal or ice-mass model: no glacio-isostasy, tectonics, sea level curve, historical topographic reconstruction, or time-sliced ice/glacier data.
2. The ETOPO tile is not downsampled or aligned to the earlier native CHELSA 101×57 grid; a later separate protocol would be needed to enforce land masks on original CHELSA LCP route cells.
3. Historic host distributions and connectivity not inferred; the climate costs remain a *modern* host-climate analogue projected into LGM CHELSA climate.
4. The original genetic one-shot 259 species remains **NOT_EVALUABLE_EMPIRICAL_RESPONSE**; no COI identity or genetic effects were opened.
5. The 665-point host climate reference versus the 26-point insect self reference remains a more substantial interpretive asymmetry than the now-checked static shoreline.

Official references: NOAA ETOPO2022 doi:10.25921/fd45-gt74 and ETOPO scientific dataset paper doi:10.5194/essd-17-1835-2025. EPSG:9518 WGS84 horizontal/EGM2008 vertical.
