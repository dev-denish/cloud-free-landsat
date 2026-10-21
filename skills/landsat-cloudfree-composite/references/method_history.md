# How the 2009 west-block (Puttalam) cloud-free image was built: history and lessons

This is the project record for `D:\Arockiaraj\LULC_2009\10_Scripts`, scripts 40–54, from 2026-10-03 to 2026-10-04. The input was 105 USGS Landsat 5 C2 L2SP bundles for 2009 (D2 package, `00_Source_USGS_bundles`, read through `/vsitar/` without extracting). The block is paths 141/142, rows 054/055, about 4,300 km².

## Attempts in order

| Script | Approach | Result / user verdict |
|---|---|---|
| 40 | Screen every 2009 scene over the block (QA) | No cloud-free 2009 pass exists. Best is 30 Jan 2009 (69 % USGS-clear). |
| 41 | One pass, cloud-masked | Too many holes |
| 42 | Best pass + greedy single-date fills + thermal cloud/shadow test | User: "even worse" (specks, haze) |
| 43 | Per-pixel median of all clear 2009 observations | Seamless, but rejected (flat, not real values) |
| 44 | All-month 3 km tiles, cleanest date, feather blend; defines `read_date` (QA + thermal) | Still blotchy and hazy. Contact sheet: every date has popcorn cumulus. |
| 45 | Dry season only (Jan–Apr), 1.5 km tiles, local Wallis match to the dry-season median | Much better |
| 46/47 | Planted polygons: each polygon from ONE real scene that is 100 % clear | 775 polygons from one scene, 15 from 2+, 0 unclear |
| 48 (v2) | Temporal shadow/cloud test; land/water matched separately; gaps = the real value closest to the reference; real values kept | Grey cloud-like patches by the lagoon |
| 49 (v3) | Dark-object haze test (5th pct blue over vegetation per 1.5 km tile, relative to the clearest date) | Haze fixed, but fragmented around the lagoon |
| 50 (v4) | Date pool widened to Jan–Sep (dry zone: 6 May, 25 Jul and 11 Sep clean), 0.03 off-season penalty | Clean lagoon; small specks left |
| 51 (v5) | USGS confidence ≥ medium masked; speck test (≤ 300 px blobs vs the clean median) | Specks gone; pale straight-edged polygons on the Kalpitiya spit |
| 52 (v6) | Wallis only where a real reference exists (≥ 3 % of the window) | Polygons gone; sea patches remain |
| 53 (v7) | Sea harmonised: one smooth sea-colour field + per-date fine detail | **ACCEPTED** |
| 54 | Burn v7 into the island mosaic: exact DN, provenance, replaced mask | Analysis mosaic (`Mosaic_new`) |

On 2026-10-04 the user chose the harmonised (colour-matched) west block for analysis. With real values, R7 "needs review" clusters went from 91 to 135; with the harmonised block they were 117.

## Lessons

- A thermal test alone flags cool water and mangrove as cloud. Require that the pixel also *looks* like cloud.
- Thin haze passes the USGS QA and thermal masks. Test haze relatively, across dates at the same place (dark object).
- A spatial-only speck detector flags bare fields. Compare against other dates of the same pixel instead.
- Wallis matching to a nearest-pixel-filled reference makes straight-edged Voronoi artefacts. Match only where the reference is real.
- Sea tone changes with tide and turbidity from date to date. Harmonise the water low-pass across dates and keep each date's detail.
- In dry-zone areas, "dry season only" can be too strict around lagoons. Widen the pool and penalise off-season dates instead.
- Keep two products: exact real values (provenance) and a harmonised copy (classification).
- Show the user zoomed tiles. They judge by eye, often on a phone.
