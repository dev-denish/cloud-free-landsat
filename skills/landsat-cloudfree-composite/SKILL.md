---
name: landsat-cloudfree-composite
description: Build a cloud-, shadow- and haze-free Landsat (TM/ETM+/OLI, USGS Collection-2 Level-2) image for an area from many dates, where every pixel is ONE real observation (no averaging) plus a colour-matched copy for classification. Use when the user asks for a cloud-free / clean / haze-free image, mosaic or composite, wants clouds, cumulus specks, shadows or haze removed, wants to patch or replace a cloudy part of a mosaic, or wants plots/polygons taken from a single clean scene. Works from USGS .tar bundles or extracted folders; runs with the QGIS Python.
---

# Landsat cloud-free composite (one real date per pixel)

This method came out of the 2009 LULC project (Puttalam west block, Sri Lanka, Oct 2026). The user rejected these approaches as hazy, blotchy or speckled: single best pass plus fills, full-year median, and all-month tile patching. What finally worked was **v7**, implemented in `scripts/cloudfree_composite.py`. The history and the reasons are in `references/method_history.md`; read it before changing thresholds.

## When to use what

| User wants | Use |
|---|---|
| Clean image for an AOI from a stack of Landsat bundles | `cloudfree_composite.py` (this skill) |
| Only to know which dates are usable | `--screen-only` (QA-only, seconds) |
| Plots / polygons each from ONE clean scene (field checks, plot statistics) | add `--keep-whole polygons.shp` |
| Composite for classification / K-Means without date seams | use `<prefix>_sr_harmonised.tif` |
| Exact USGS values (indices, reporting, provenance) | use `<prefix>_dn_real.tif` / `_sr_real.tif` + `_source.tif` |

## How to run

Python: the QGIS Python has GDAL + numpy + scipy. On this PC (PowerShell):

```powershell
& "C:\Program Files\QGIS 3.44.9\bin\python-qgis-ltr.bat" "$HOME\.claude\skills\landsat-cloudfree-composite\scripts\cloudfree_composite.py" `
  --bundles "D:\path\to\USGS_bundles" --aoi "D:\path\aoi.shp" --crs EPSG:32644 `
  --years 2009 --season 0101-0930 --prefer 0101-0430 --out "D:\path\out" --prefix myarea
```

The AOI can be given three ways: `--aoi` as a vector (needs `--crs`), `--aoi` as a raster mask (that raster's grid is used), or `--bounds XMIN YMIN XMAX YMAX --crs ...`.

Always do this in two steps:
1. Run `--screen-only` first, then show the user the `_screen.csv` table: date, coverage and USGS-clear share. Agree on `--season` (the date pool) and `--prefer` (the dry or least-cloudy season) with them.
2. Then run the full composite. Its speed is about 600 km² × 25 dates in 2.5 min, and the 4,300 km² west block took about 40 min. Memory is roughly n_dates × pixels × 15 bytes. Split AOIs that would need more than ~8 GB.

Key options (the defaults are the accepted v7 values):
- `--season MMDD-MMDD`: the date pool. It may wrap the year end.
- `--prefer MMDD-MMDD`: dates outside it get the `--off-pen 0.03` score penalty.
- `--min-date 0.02`: drops dates with less of the AOI clear than this.
- `--tile 50`: 1.5 km patches.
- `--min-share 0.90`: the clean share a tile needs.
- `--feather 6` / `--feather-water 25`: feathering in px.
- `--haze-exc 0.010`: haze test threshold, as dark-object blue above the clearest date.
- `--match-win 101`: Wallis window.
- `--no-qa-conf`: skip the USGS confidence ≥ medium mask.
- `--no-sea-harmonise`: skip the sea-colour harmonising.

## Method (what each step catches)

1. **USGS QA_PIXEL** masks fill, dilated cloud, cloud and shadow (plus cirrus on OLI). The cloud/shadow **confidence ≥ medium** bits (9 and 11) are also masked, grown by 1 px; these catch the faint cumulus the flag bits miss.
2. **Thermal cloud test.** A pixel counts as cloud only if it is colder than the 51 px local land mean **and** looks like cloud: SWIR1 > 0.12, NDVI < 0.25, not water. Without the "looks like cloud" part, cool lagoon water and wet mangrove were flagged in every scene. A shadow search follows within 25 px of clouds.
3. **Temporal outlier test**, 2 passes against the per-pixel median of the clear dates. Shadow = NIR and SWIR1 below 0.75 × median; cloud = blue above median + 0.025.
4. **Dark-object haze test.** For each 1.5 km tile, take the 5th percentile of blue over clear vegetation. A tile is hazy when that value is more than 0.010 above the clearest date at the same tile. The test is relative only, because land cover changes the level. Thin haze passes all QA, thermal and temporal masks and caused the "grey cloud-like patches".
5. **Speck test.** This needs at least 3 clean observations. It flags whitish or dark blobs of ≤ 300 px compared with the median of the other clean dates. Bigger blobs are real change.
6. **Patch choice per tile.** The score is clean share − haze − HOT rank − off-season penalty. The weights are feathered, and the real value comes from the date with the largest weight.
7. **Harmonised copy.** Each date gets a Wallis match to the clean median, with land and water matched separately. The match is only applied where a real reference covers at least 3 % of the window; otherwise the edges of the nearest-pixel fill show up as straight-edged false polygons. For water, a single smooth sea-colour field is used (the mean of the per-date low-pass) plus each date's own fine detail. This removes tide and turbidity patches.
8. **Gaps.** Tier 1 takes the least-hazy clear observation. Tier 2 takes the observation closest to the clean reference. Neither tier averages.
9. **`--keep-whole`.** Each polygon (+1 px ring) gets one date: the most clear pixels, then the least hazy.

## Outputs (`--out`, `--prefix`)

| File | Contents |
|---|---|
| `_sr_real.tif` | Float32 SR, Blue Green Red NIR SWIR1 SWIR2 |
| `_dn_real.tif` | UInt16 exact USGS DN (SR = DN × 0.0000275 − 0.2) |
| `_sr_harmonised.tif` | Colour-matched copy |
| `_source.tif` | Date index per pixel (category names = dates) |
| `_dates.csv` | Per-date clear share and share of the composite |
| `_truecolour.tif`, `_falsecolour543.tif`, `_preview.png` | Display images |
| `_screen.csv`, `_summary.txt` | Screening table and run summary |

## Checking the result (always, before reporting success)

- Read `_summary.txt`: "still empty" should be 0, and the tier-1/tier-2 gap counts should be small.
- Look at the result as **zoomed tiles**, not only the whole preview. Small cumulus specks, shadow dots and haze bands are only visible zoomed. The user often checks from a phone; send zoomed PNG tiles with SendUserFile.
- Look for: white/grey specks; dark dots next to them; grey haze bands; straight-edged pale polygons (a reference-fill artefact); rounded colour patches in sea or lagoon (date tone); and square tile steps over open water.
- If an older accepted product exists, compare against it (identical-value share, side-by-side PNG).

## Known limits / pitfalls

- Over **open lagoon or sea**, faint square tile steps and dark blobs can remain. Some are real (shoals, seagrass), some are date effects. Use the harmonised copy, or raise `--feather-water`. Do not treat water texture as cloud without checking several dates.
- A mosaic in which every pixel is a different real date gives **date-boundary artefacts in K-Means/R7**. In the project, flagged clusters went from 91 to 135; the user then chose the harmonised copy for analysis. Keep the real values for provenance.
- **QGIS file locks:** if an output raster is open in QGIS, a re-run fails with "Permission denied". Remove those layers first.
- The PC sleeps overnight and pauses long runs; keep the session active.
- The PowerShell tool can block inline scripts. Write scripts to a file with the Write tool, then run them. Run Python that reads extracted downloads with `-I`.
- To burn the composite into a bigger existing mosaic, use the same CRS and grid. Re-read exact DN from `_dn_real.tif` (or convert harmonised SR as DN = round((SR + 0.2) / 0.0000275), min 1). Replace only where the base mosaic is valid, and record the provenance (source date and a replaced mask). Project script for this: `D:\Arockiaraj\LULC_2009\10_Scripts\54_build_new_mosaic_2009.py`.
- Sentinel-2 is not supported: the bands, QA and thermal inputs differ (it has no thermal band). Adapting it would need SCL or s2cloudless masks in place of steps 1–2.
