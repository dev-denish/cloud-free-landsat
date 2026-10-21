#!/usr/bin/env python
"""Cloud-, shadow- and haze-free Landsat Collection-2 Level-2 composite for an area of interest.

Generalised from the 2009 LULC west-block (Puttalam) work, scripts 44 + 50-53 (method "v7", accepted by the user 2026-10-03).
Every output pixel is ONE real observation of ONE date (no averaging); a second, colour-matched copy is written for classifiers.

Pipeline
  1. per acquisition date (all path/rows of that date merged, first bundle with data wins):
       USGS QA_PIXEL fill / dilated cloud / cloud / shadow (+ cirrus on OLI) + cloud/shadow CONFIDENCE >= medium (grown 1 px)
       + thermal cloud test (cold vs 51-px local land mean AND looks like cloud: SWIR1 > 0.12, NDVI < 0.25, not water)
       + geometric shadow test near clouds (NIR and SWIR1 < 0.65 x local mean), all grown GROW px
     dates with < MIN_DATE of the AOI clear are dropped
  2. temporal outlier test (2 passes) against the per-pixel median of clear dates: shadow (NIR & SWIR1 < 0.75 x median, not
     water) / cloud (blue > median + 0.025)
  3. haze: dark-object blue (5th pct of SR blue over clear veg, per TILE px tile); a tile of a date is HAZY when it is > DO_EXC
     above the clearest date at that tile (relative only - the level depends on land cover)
  4. speck test against the clean median (needs >= 3 clean obs): small whitish / dark blobs <= 300 px
  5. per tile choose the date with the best score (clean share - haze - HOT rank - off-season penalty); feathered weights
     (land FEATHER px, water FEATHER_W px) x the date's clean mask; the REAL value = date with the largest weight
  6. display / harmonised copy: Wallis match (101 px) of each date to the clean median, land vs water separately, only where a
     real reference exists (>= 3 % of the window); water = one smooth sea-colour field (confidence-weighted mean of the per-date
     water low-pass, sigma ~1.2 km) + the date's own fine detail
  7. gaps: tier 1 = least-hazy clear observation; tier 2 = observation closest to the clean reference (one real value)
  8. optional --keep-whole polygons (+1 px ring): every polygon takes ONE date (most clear pixels, then least hazy), real values

Outputs in --out (prefix --prefix):
  _sr_real.tif         Float32 surface reflectance, one real date per pixel (6 bands Blue Green Red NIR SWIR1 SWIR2), 0 = no data
  _dn_real.tif         UInt16 exact USGS DN of the same values (SR = DN x 0.0000275 - 0.2)
  _sr_harmonised.tif   Float32 colour-matched copy (no date seams; use for classifiers / K-Means if seams cause artefacts)
  _source.tif          date index per pixel (category names = dates), _dates.csv (index, date, bundles, clear share, share used)
  _truecolour.tif / _falsecolour543.tif / _preview.png   display RGB (from the harmonised copy)
  _screen.csv, _summary.txt

Needs GDAL, numpy, scipy (e.g. the QGIS Python: "C:\\Program Files\\QGIS 3.xx\\bin\\python-qgis-ltr.bat").
Memory ~ n_dates x H x W x 15 bytes (30 dates x 5 M px ~ 2.3 GB). Split very large AOIs into blocks.
"""
import argparse
import csv
import glob
import os
import re
import sys
import time
from types import SimpleNamespace

import numpy as np
from osgeo import gdal, osr
from scipy import ndimage

gdal.UseExceptions()

SR_SCALE, SR_OFF = 0.0000275, -0.2
ST_SCALE, ST_OFF = 0.00341802, 149.0
BAND_NAMES = ("Blue", "Green", "Red", "NIR", "SWIR1", "SWIR2")
SENSOR_BANDS = {  # canonical order Blue Green Red NIR SWIR1 SWIR2, thermal band
    "TM": (("SR_B1", "SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B7"), "ST_B6"),      # LT04 / LT05 / LE07
    "OLI": (("SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7"), "ST_B10"),    # LC08 / LC09
}
PID_RE = re.compile(r"(L[CTEO]0[4-9]_L2S[PR]_(\d{6})_(\d{8})_\d{8}_\d{2}_T[12RT])")
PAL = [(31, 120, 180), (51, 160, 44), (227, 26, 28), (255, 127, 0), (106, 61, 154), (177, 89, 40), (166, 206, 227),
       (178, 223, 138), (251, 154, 153), (253, 191, 111), (202, 178, 214), (255, 255, 153), (0, 0, 0)]


def sensor_of(pid):
    return "OLI" if pid[:4] in ("LC08", "LC09", "LO08", "LO09") else "TM"


# ----------------------------------------------------------------------------------------------------------------- inputs
def find_bundles(root):
    """{product_id: (kind, location)} for USGS C2 L2 bundles: *.tar files or extracted folders (*_QA_PIXEL.TIF)."""
    out = {}
    for p in glob.glob(os.path.join(root, "**", "*"), recursive=True):
        name = os.path.basename(p)
        m = PID_RE.match(name)
        if not m:
            continue
        pid = m.group(1)
        if name.endswith("_QA_PIXEL.TIF"):
            out[pid] = ("dir", os.path.dirname(p))                 # extracted copy preferred (faster)
        elif name.lower().endswith(".tar") and pid not in out:
            out[pid] = ("tar", p)
    return out


def band_path(kind, loc, pid, name):
    return f"/vsitar/{loc}/{pid}_{name}.TIF" if kind == "tar" else os.path.join(loc, f"{pid}_{name}.TIF")


def in_window(mmdd, win):
    a, b = win
    return a <= mmdd <= b if a <= b else (mmdd >= a or mmdd <= b)


def make_grid(a, out):
    """AOI -> grid (bounds, wkt, res, H, W) + boolean mask."""
    res = a.res
    if a.aoi and a.aoi.lower().endswith((".tif", ".tiff", ".img", ".vrt")):
        ds = gdal.Open(a.aoi); gt = ds.GetGeoTransform(); b = ds.GetRasterBand(1)
        v = b.ReadAsArray(); nd = b.GetNoDataValue()
        mask = (v != 0) if nd is None else ((v != 0) & (v != nd))
        H, W = mask.shape
        return SimpleNamespace(bounds=(gt[0], gt[3] + H * gt[5], gt[0] + W * gt[1], gt[3]), wkt=ds.GetProjection(),
                               res=gt[1], H=H, W=W, mask=mask)
    if not a.crs:
        sys.exit("--crs is required with a vector AOI or --bounds (projected CRS in metres, e.g. EPSG:32644)")
    srs = osr.SpatialReference(); srs.SetFromUserInput(a.crs); wkt = srs.ExportToWkt()
    if a.bounds:
        x0, y0, x1, y1 = a.bounds
        vec = None
    else:
        vec = os.path.join(out, f"{a.prefix}_aoi_reprojected.gpkg")
        gdal.VectorTranslate(vec, a.aoi, format="GPKG", dstSRS=a.crs, accessMode="overwrite")
        lyr = gdal.OpenEx(vec).GetLayer(0); x0, x1, y0, y1 = lyr.GetExtent()
    x0 = np.floor(x0 / res) * res; y0 = np.floor(y0 / res) * res
    x1 = np.ceil(x1 / res) * res; y1 = np.ceil(y1 / res) * res
    W, H = int(round((x1 - x0) / res)), int(round((y1 - y0) / res))
    if vec:
        md = gdal.GetDriverByName("MEM").Create("", W, H, 1, gdal.GDT_Byte)
        md.SetGeoTransform((x0, res, 0, y1, 0, -res)); md.SetProjection(wkt)
        gdal.Rasterize(md, vec, burnValues=[1], allTouched=True)
        mask = md.ReadAsArray().astype(bool)
    else:
        mask = np.ones((H, W), bool)
    return SimpleNamespace(bounds=(x0, y0, x1, y1), wkt=wkt, res=res, H=H, W=W, mask=mask)


def rasterize_ids(path, g, out, prefix):
    """polygon id raster (1..N, 0 = none) on the grid, ALL_TOUCHED, + 1 px ring given to the neighbouring polygon."""
    vec = os.path.join(out, f"{prefix}_keepwhole_reprojected.gpkg")
    gdal.VectorTranslate(vec, path, format="GPKG", dstSRS=g.wkt, accessMode="overwrite", layerName="kw")
    md = gdal.GetDriverByName("MEM").Create("", g.W, g.H, 1, gdal.GDT_Int32)
    md.SetGeoTransform((g.bounds[0], g.res, 0, g.bounds[3], 0, -g.res)); md.SetProjection(g.wkt)
    gdal.Rasterize(md, vec, SQLStatement='SELECT fid AS rid, * FROM "kw"', attribute="rid", allTouched=True)
    ids = md.ReadAsArray()
    return np.where(ids > 0, ids, ndimage.grey_dilation(ids, size=(3, 3)))


# ------------------------------------------------------------------------------------------------------------ per date
def to_sr(dn):
    return np.where(dn > 0, dn.astype(np.float32) * SR_SCALE + SR_OFF, 0).astype(np.float32)


def local_mean(v, valid, win):
    f = valid.astype(np.float32)
    return ndimage.uniform_filter(np.where(valid, v, 0).astype(np.float32), win) / np.maximum(ndimage.uniform_filter(f, win), 1e-3)


def warp(kind, loc, pid, name, g, init):
    ds = gdal.Warp("", band_path(kind, loc, pid, name), format="MEM", outputBounds=g.bounds, xRes=g.res, yRes=g.res,
                   dstSRS=g.wkt, resampleAlg="near", srcNodata=None, dstNodata=None, warpOptions=[f"INIT_DEST={init}"])
    return ds.GetRasterBand(1).ReadAsArray()[:g.H, :g.W]


def qa_bad(q, pid, use_conf):
    flags = (((q >> 1) | (q >> 3) | (q >> 4)) & 1).astype(bool)      # dilated cloud, cloud, cloud shadow
    if sensor_of(pid) == "OLI":
        flags |= ((q >> 2) & 1).astype(bool)                          # cirrus
    conf = (((q >> 9) & 1) | ((q >> 11) & 1)).astype(bool) if use_conf else np.zeros_like(flags)
    return flags, conf


def screen_date(items, g, p):
    """QA-only clear share of the AOI (fast)."""
    have = np.zeros((g.H, g.W), bool); bad = np.zeros((g.H, g.W), bool)
    for kind, loc, pid in items:
        q = warp(kind, loc, pid, "QA_PIXEL", g, 1).astype(np.uint16)
        ok = ((q & 1) == 0) & ~have
        if not ok.any():
            continue
        f, c = qa_bad(q, pid, p.qa_conf)
        bad[ok] = (f | ndimage.binary_dilation(c, iterations=1))[ok]
        have |= ok
    m = g.mask
    return (have & m).sum() / m.sum(), (have & ~bad & m).sum() / m.sum()


def read_date(items, g, p):
    """exact DN (6,H,W) uint16, clear mask, have mask for one acquisition date (all bundles of that date merged)."""
    H, W = g.H, g.W
    dn = np.zeros((6, H, W), np.uint16); bt = np.zeros((H, W), np.float32)
    bad = np.zeros((H, W), bool); have = np.zeros((H, W), bool)
    for kind, loc, pid in items:
        bands, st = SENSOR_BANDS[sensor_of(pid)]
        q = warp(kind, loc, pid, "QA_PIXEL", g, 1).astype(np.uint16)
        ok = ((q & 1) == 0) & ~have
        if not ok.any():
            continue
        for i, b in enumerate(bands):
            dn[i][ok] = warp(kind, loc, pid, b, g, 0)[ok]
        try:
            t = warp(kind, loc, pid, st, g, 0).astype(np.float32)
            bt[ok] = np.where(t[ok] > 0, t[ok] * ST_SCALE + ST_OFF, 0)
        except RuntimeError:
            pass                                                      # L2SR product: no thermal band, QA only
        f, c = qa_bad(q, pid, p.qa_conf)
        bad[ok] = f[ok]
        bad[ok] |= ndimage.binary_dilation(c & ok, iterations=1)[ok]
        have |= ok
    if not have.any():
        return None
    have &= dn[3] > 0
    sr = to_sr(dn)
    blue, green, red, nir, sw1 = sr[0], sr[1], sr[2], sr[3], sr[4]
    water = (green - sw1) / np.maximum(green + sw1, 1e-6) > 0
    cloud = bad.copy()
    if (bt > 0).any():
        ref_bt = local_mean(bt, have & (bt > 0) & ~water & (blue < p.blue_cloud) & ~bad, p.local_win)
        ndvi = (nir - red) / np.maximum(nir + red, 1e-6)
        # thermal cloud only where it also LOOKS like cloud (cool lagoon water / wet mangrove were flagged in every scene otherwise)
        looks = ~water & (sw1 > p.cloud_swir) & (ndvi < p.cloud_ndvi)
        th = have & (bt > 0) & looks & (((bt < ref_bt - p.cold_k) & (blue >= p.blue_cloud)) | ((blue > 0.2) & (bt < ref_bt - 1.5)))
        cloud = ndimage.binary_opening(th, np.ones((2, 2), bool)) | bad
    near = ndimage.binary_dilation(cloud, iterations=p.shadow_dist)
    shadow = (have & near & ~cloud & ~water & (nir < p.shadow_frac * local_mean(nir, have & ~near, p.local_win))
              & (sw1 < p.shadow_frac * local_mean(sw1, have & ~near & ~water, p.local_win)))
    clear = have & ~ndimage.binary_dilation(cloud | shadow, iterations=p.grow)
    return dn, clear, have


# ------------------------------------------------------------------------------------------------------- temporal steps
def median_ref(DN, CL, H, W, chunk=256):
    ref = np.zeros((6, H, W), np.float32)
    for r0 in range(0, H, chunk):
        r1 = min(H, r0 + chunk)
        for i in range(6):
            st = np.stack([np.where(CL[k][r0:r1], to_sr(DN[k][i, r0:r1]), np.nan) for k in range(len(DN))])
            with np.errstate(all="ignore"):
                ref[i, r0:r1] = np.nan_to_num(np.nanmedian(st, axis=0), nan=0)
    return ref


def temporal_flags(x, ref, have, p):
    water = (ref[1] - ref[4]) / np.maximum(ref[1] + ref[4], 1e-6) > 0
    ok = have & (ref[3] > 0)
    shadow = ok & ~water & (x[3] < p.shadow_r * ref[3]) & (x[4] < p.shadow_r * ref[4])
    cloud = ok & (x[0] > ref[0] + p.cloud_db)
    return ndimage.binary_dilation(ndimage.binary_opening(shadow | cloud, np.ones((2, 2), bool)), iterations=p.tgrow)


def speck_flags(x, ref, have, nclean, cb=0.020, cv=0.020, sh=0.80, maxblob=300):
    ok = have & (nclean >= 3) & (ref[3] > 0)
    water = (ref[1] - ref[4]) / np.maximum(ref[1] + ref[4], 1e-6) > 0
    vis, refv = x[:3].mean(0), ref[:3].mean(0)
    cloud = ok & (x[0] > ref[0] + cb) & (vis > refv + cv) & ((x[2] - x[0]) < 0.08)
    shadow = ok & ~water & (x[3] < sh * ref[3]) & (x[4] < sh * ref[4]) & (x[1] < 0.88 * ref[1]) & (x[2] < 0.88 * ref[2])
    lab, _ = ndimage.label(cloud | shadow, np.ones((3, 3)))
    sz = np.bincount(lab.ravel()); keep = sz <= maxblob; keep[0] = False
    return ndimage.binary_dilation(keep[lab], iterations=2)


def water_lowpass(x, wm, f=4, sigma=10):
    H, W = wm.shape; h2, w2 = -(-H // f), -(-W // f)
    pad = ((0, h2 * f - H), (0, w2 * f - W))
    cover = np.pad(wm.astype(np.float32), pad).reshape(h2, f, w2, f).mean((1, 3))
    gw = ndimage.gaussian_filter(cover, sigma)
    out = []
    for i in range(6):
        mi = np.pad(np.where(wm, x[i], 0).astype(np.float32), pad).reshape(h2, f, w2, f).mean((1, 3))
        out.append(ndimage.zoom(ndimage.gaussian_filter(mi, sigma) / np.maximum(gw, 1e-6), f, order=1)[:H, :W])
    return np.stack(out), ndimage.zoom(gw, f, order=1)[:H, :W]


def lmean(v, w, win):
    return ndimage.uniform_filter(np.where(w, v, 0).astype(np.float32), win) / np.maximum(ndimage.uniform_filter(w.astype(np.float32), win), 1e-6)


def wallis(x, ref, w, win, minfrac=0.03):
    nb = ndimage.uniform_filter(w.astype(np.float32), win)
    mx, mr = lmean(x, w, win), lmean(ref, w, win)
    sx = np.sqrt(np.maximum(lmean(x * x, w, win) - mx * mx, 1e-8))
    sr_ = np.sqrt(np.maximum(lmean(ref * ref, w, win) - mr * mr, 1e-8))
    return np.where(nb >= minfrac, mr + np.clip(sr_ / sx, 0.5, 2.0) * (x - mx), x)


def dark_object(DN, CL, tid, ntile, ty, tx, p):
    n = len(DN); meas = np.full((n, ntile), np.nan, np.float32)
    for k in range(n):
        b1, red, nir = to_sr(DN[k][0]), to_sr(DN[k][2]), to_sr(DN[k][3])
        v = CL[k] & ((nir - red) / np.maximum(nir + red, 1e-6) > 0.25) & (b1 > 0)
        ids, vals = tid[v], b1[v]
        o = np.argsort(ids, kind="stable"); ids, vals = ids[o], vals[o]
        st = np.searchsorted(ids, np.arange(ntile)); en = np.searchsorted(ids, np.arange(ntile), side="right")
        for t in np.flatnonzero(en - st >= p.do_min):
            meas[k, t] = np.percentile(vals[st[t]:en[t]], 5)

    def fill(gr):
        gr = gr.reshape(ty, tx); bad = np.isnan(gr)
        if bad.all():
            return np.zeros(ntile, np.float32)
        if bad.any():
            gr = gr[tuple(ndimage.distance_transform_edt(bad, return_distances=False, return_indices=True))]
        return ndimage.median_filter(gr, 3).reshape(-1)
    with np.errstate(all="ignore"):
        base = np.nanmin(meas, axis=0)
    base = np.minimum(fill(base), p.do_base_max)
    return np.stack([fill(meas[k]) - base for k in range(n)]).astype(np.float32), meas


# --------------------------------------------------------------------------------------------------------------- output
def save(path, arr, dtype, names, g, meta, nodata=0, rgb=False, cats=None):
    nb = arr.shape[0] if arr.ndim == 3 else 1
    o = gdal.GetDriverByName("GTiff").Create(path, g.W, g.H, nb, dtype, ["COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=IF_SAFER"]
                                             + (["PHOTOMETRIC=RGB"] if rgb else []))
    o.SetGeoTransform((g.bounds[0], g.res, 0, g.bounds[3], 0, -g.res)); o.SetProjection(g.wkt); o.SetMetadata(meta)
    for i in range(nb):
        b = o.GetRasterBand(i + 1); b.SetNoDataValue(nodata); b.SetDescription(names[i])
        if cats:
            ct = gdal.ColorTable()
            for j in range(len(cats)):
                ct.SetColorEntry(j + 1, (*PAL[j % len(PAL)], 255))
            b.SetRasterColorTable(ct); b.SetRasterCategoryNames([""] + cats)
        b.WriteArray(arr[i] if arr.ndim == 3 else arr)
    o.BuildOverviews("NEAREST" if cats else "AVERAGE", [2, 4, 8, 16]); o = None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bundles", required=True, help="folder with USGS C2 L2 *.tar bundles or extracted folders (searched recursively)")
    ap.add_argument("--aoi", help="AOI: vector file (shp/gpkg/geojson) or a raster mask (non-zero = inside; its grid is used)")
    ap.add_argument("--bounds", nargs=4, type=float, metavar=("XMIN", "YMIN", "XMAX", "YMAX"), help="AOI box in --crs units")
    ap.add_argument("--crs", help="output CRS (projected, metres), e.g. EPSG:32644; not needed with a raster AOI")
    ap.add_argument("--res", type=float, default=30.0)
    ap.add_argument("--out", required=True)
    ap.add_argument("--prefix", default="cloudfree")
    ap.add_argument("--years", default="", help="comma list of years to use, e.g. 2009 (default: all)")
    ap.add_argument("--season", default="0101-1231", help="date pool MMDD-MMDD (may wrap the year end, e.g. 1101-0331)")
    ap.add_argument("--prefer", default="", help="preferred window MMDD-MMDD (dry season); other dates get --off-pen (default none)")
    ap.add_argument("--off-pen", type=float, default=0.03)
    ap.add_argument("--min-date", type=float, default=0.02, help="drop dates with less than this share of the AOI clear")
    ap.add_argument("--tile", type=int, default=50, help="patch / haze tile size in px (50 px = 1.5 km at 30 m)")
    ap.add_argument("--min-share", type=float, default=0.90)
    ap.add_argument("--feather", type=float, default=6); ap.add_argument("--feather-water", type=float, default=25)
    ap.add_argument("--match-win", type=int, default=101)
    ap.add_argument("--haze-exc", type=float, default=0.010, help="dark-object blue excess over the clearest date = hazy")
    ap.add_argument("--no-qa-conf", action="store_true", help="do not mask USGS cloud/shadow confidence >= medium")
    ap.add_argument("--no-sea-harmonise", action="store_true")
    ap.add_argument("--keep-whole", help="polygons (vector) that must each come from ONE date (e.g. plots / planted polygons)")
    ap.add_argument("--screen-only", action="store_true", help="only write the per-date QA clear share table and stop")
    a = ap.parse_args()
    if not a.aoi and not a.bounds:
        sys.exit("give --aoi or --bounds")
    p = SimpleNamespace(qa_conf=not a.no_qa_conf, blue_cloud=0.08, local_win=51, cloud_swir=0.12, cloud_ndvi=0.25, cold_k=3.0,
                        shadow_dist=25, shadow_frac=0.65, grow=3, shadow_r=0.75, cloud_db=0.025, tgrow=2,
                        do_min=150, do_base_max=0.038)
    os.makedirs(a.out, exist_ok=True)
    pre = os.path.join(a.out, a.prefix)
    t0 = time.time()
    log = lambda s: print(f"[{time.time() - t0:6.0f}s] {s}", flush=True)

    g = make_grid(a, a.out)
    m = g.mask; H, W = g.H, g.W
    log(f"AOI {m.sum() * g.res ** 2 / 1e6:,.0f} km2, grid {W} x {H} px at {g.res:g} m")
    season = tuple(a.season.split("-")); prefer = tuple(a.prefer.split("-")) if a.prefer else None
    years = {y.strip() for y in a.years.split(",") if y.strip()}
    bund = find_bundles(a.bundles)
    by_date = {}
    for pid, (kind, loc) in sorted(bund.items()):
        d = pid[17:25]
        if (years and d[:4] not in years) or not in_window(d[4:], season):
            continue
        by_date.setdefault(d, []).append((kind, loc, pid))
    log(f"{len(bund)} bundles found, {len(by_date)} acquisition dates in the year / season filter")
    if not by_date:
        sys.exit("no bundles match")

    # ---- screening (QA only, fast)
    screen = []
    for d, items in sorted(by_date.items()):
        cov, clr = screen_date(items, g, p)
        if cov > 0:
            screen.append((d, cov, clr, ";".join(i[2] for i in items)))
            log(f"  screen {d}: covers {cov:.1%}, USGS-clear {clr:.1%}")
    with open(pre + "_screen.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["date", "coverage", "usgs_clear_share", "bundles"]); w.writerows(screen)
    if a.screen_only:
        log(f"screen written: {pre}_screen.csv"); return
    use = [s[0] for s in screen if s[2] >= a.min_date]
    est = len(use) * H * W * 15 / 1e9
    log(f"{len(use)} dates with >= {a.min_date:.0%} USGS-clear; memory estimate ~{est:.1f} GB")

    # ---- full read
    DN, CL0, HAVE, dates = [], [], [], []
    for d in use:
        r = read_date(by_date[d], g, p)
        if r is None:
            continue
        dn, clear, have = r
        cm = clear & m
        if cm.sum() < a.min_date * m.sum():
            log(f"  skip {d}: clear {cm.sum() / m.sum():.1%} after thermal/shadow tests"); continue
        DN.append(dn); CL0.append(cm); HAVE.append(have & m); dates.append(d)
        log(f"  read {d}: clear {cm.sum() / m.sum():.1%}")
    n = len(DN)
    if n == 0:
        sys.exit("no usable dates")
    CL = np.stack(CL0)
    for it in range(2):
        ref = median_ref(DN, CL, H, W); added = 0
        for k in range(n):
            f = temporal_flags(to_sr(DN[k]), ref, HAVE[k], p)
            added += int((CL[k] & f).sum()); CL[k] = CL0[k] & ~f
        log(f"temporal test pass {it + 1}: {added:,} clear px flagged as shadow / cloud")
    ty, tx = -(-H // a.tile), -(-W // a.tile)
    tid = (np.arange(H)[:, None] // a.tile) * tx + (np.arange(W)[None, :] // a.tile)
    ntile = ty * tx
    do, do_meas = dark_object(DN, CL, tid, ntile, ty, tx, p)
    hazy = [(do[k] > a.haze_exc)[tid] for k in range(n)]
    CLN = np.stack([CL[k] & ~hazy[k] for k in range(n)])
    for k in range(n):
        log(f"  haze {dates[k]}: {1 - CLN[k].sum() / max(CL[k].sum(), 1):.0%} of its clear area hazy")
    ref = median_ref(DN, CLN, H, W); ncln = CLN.sum(0); nspeck = 0
    for k in range(n):
        x_ = to_sr(DN[k])
        f = temporal_flags(x_, ref, HAVE[k], p)
        sp = speck_flags(x_, ref, HAVE[k] & CL0[k], ncln)
        nspeck += int((sp & CL0[k] & ~f).sum())
        CL[k] = CL0[k] & ~(f | sp); CLN[k] = CL[k] & ~hazy[k]
    log(f"speck test: {nspeck:,} more px flagged")
    has_ref = CLN.any(0)
    ref = median_ref(DN, CLN, H, W)
    if (m & ~has_ref).any():
        idx = tuple(ndimage.distance_transform_edt(~has_ref, return_distances=False, return_indices=True))
        for i in range(6):
            ref[i] = np.where(has_ref, ref[i], ndimage.gaussian_filter(ref[i][idx], 6))
    water = m & ((ref[1] - ref[4]) / np.maximum(ref[1] + ref[4], 1e-6) > 0)
    hot = np.array([np.median((to_sr(DN[k][0]) - 0.5 * to_sr(DN[k][2]))[CL[k]]) if CL[k].any() else 9 for k in range(n)])

    Tfield, has_T = None, np.zeros((H, W), bool)
    if not a.no_sea_harmonise:
        Tsum = np.zeros((6, H, W), np.float32); Tw = np.zeros((H, W), np.float32)
        for k in range(n):
            x_ = to_sr(DN[k]); wm_ = CL[k] & ((x_[1] - x_[4]) / np.maximum(x_[1] + x_[4], 1e-6) > 0)
            if wm_.sum() < 2000:
                continue
            lp_, wt_ = water_lowpass(x_, wm_); c_ = np.clip(wt_ / 0.2, 0, 1)
            Tsum += lp_ * c_; Tw += c_
        Tfield = Tsum / np.maximum(Tw, 1e-6); has_T = Tw > 0.05
        del Tsum
        log(f"sea-colour field: water covered {(has_T & m).sum() / m.sum():.1%} of the AOI")

    def get_mt(k):
        x = to_sr(DN[k]); y = np.zeros_like(x)
        own_w = (x[1] - x[4]) / np.maximum(x[1] + x[4], 1e-6) > 0
        wk = np.where(has_ref, water, own_w)
        okm = CL[k] & has_ref
        for i in range(6):
            y[i] = np.where(x[i] != 0, np.where(wk, wallis(x[i], ref[i], okm & wk, a.match_win),
                                                wallis(x[i], ref[i], okm & ~wk, a.match_win)), 0)
        wm = CL[k] & own_w
        if Tfield is not None and wm.sum() >= 2000:
            lp, wt = water_lowpass(x, wm); sel = wm & has_T & (wt > 0.05)
            for i in range(6):
                y[i] = np.where(sel, np.maximum(Tfield[i] + (x[i] - lp[i]), 0.002), y[i])
        return y

    area = np.bincount(tid[m], minlength=ntile)
    share = np.stack([np.bincount(tid[CLN[k]], minlength=ntile) / np.maximum(area, 1) for k in range(n)])
    score = share - 5 * np.clip(do - 0.003, 0, None) - 0.01 * (hot.argsort().argsort()[:, None] / max(n - 1, 1))
    if prefer:
        score -= a.off_pen * np.array([not in_window(d[4:], prefer) for d in dates])[:, None]
    score[share < a.min_share] -= 1
    pa = np.argmax(score, 0)[tid]
    acc = np.zeros((6, H, W), np.float32); wsum = np.zeros((H, W), np.float32)
    best_w = np.zeros((H, W), np.float32); src = np.zeros((H, W), np.int16)
    for k in range(n):
        a_ = (pa == k).astype(np.float32)
        wt = np.where(water, ndimage.gaussian_filter(a_, a.feather_water), ndimage.gaussian_filter(a_, a.feather)) * CLN[k]
        if wt.any():
            acc += get_mt(k) * wt; wsum += wt
        better = wt > best_w; best_w[better] = wt[better]; src[better] = k + 1
    disp = acc / np.maximum(wsum, 1e-6); del acc
    log("patches blended")

    gap = m & (wsum < 1e-3)
    best_h = np.full((H, W), np.inf, np.float32); best_hk = np.full((H, W), -1, np.int16)
    for k in range(n):
        ok = gap & CL[k]; dpx = do[k][tid]; better = ok & (dpx < best_h)
        best_h[better] = dpx[better]; best_hk[better] = k
    tier1 = 0
    for k in range(n):
        take = gap & (best_hk == k)
        if take.any():
            disp[:, take] = get_mt(k)[:, take]; src[take] = k + 1; tier1 += int(take.sum())
    gap &= best_hk < 0
    best_d = np.full((H, W), np.inf, np.float32); best_k = np.full((H, W), -1, np.int16)
    for k in range(n):
        x = to_sr(DN[k]); ok = gap & HAVE[k]
        hz = np.clip(x[0] - ref[0], 0, None)
        dist = sum((x[i] - ref[i]) ** 2 for i in (1, 2, 3, 4)) + 9 * hz ** 2 + np.where(CL0[k], 0, 0.0004)
        better = ok & (dist < best_d); best_d[better] = dist[better]; best_k[better] = k
    tier2 = 0
    for k in range(n):
        take = gap & (best_k == k)
        if take.any():
            disp[:, take] = get_mt(k)[:, take]; src[take] = k + 1; tier2 += int(take.sum())
    empty = int((gap & (best_k < 0)).sum())
    log(f"gaps: hazy-but-clear tier {tier1:,} px, no-clear tier {tier2:,} px, still empty {empty:,} px")

    kw_note = ""
    if a.keep_whole:
        ids = rasterize_ids(a.keep_whole, g, a.out, a.prefix)
        pm = (ids > 0) & m; npoly = int(ids.max()) + 1
        tot = np.bincount(ids[pm], minlength=npoly)
        frac = np.stack([np.bincount(ids[pm & CL[k]], minlength=npoly) / np.maximum(tot, 1) for k in range(n)])
        rank = hot.argsort().argsort() / max(n - 1, 1)
        best = np.argmax(frac - 0.001 * rank[:, None], 0)
        full = int(((frac.max(0) >= 0.999) & (tot > 0)).sum())
        for k in range(n):
            sel = pm & (best[ids] == k) & CL[k]
            if sel.any():
                disp[:, sel] = get_mt(k)[:, sel]; src[sel] = k + 1
        kw_note = f"keep-whole polygons: {int((tot > 0).sum())} in AOI, {full} fully clear on one date ({pm.sum():,} px)\n"
        log(kw_note.strip())

    # real values from the source date (exact DN)
    dn_out = np.zeros((6, H, W), np.uint16)
    for k in range(n):
        s_ = src == k + 1
        dn_out[:, s_] = DN[k][:, s_]
    dn_out[:, ~m] = 0; disp[:, ~m] = 0; src[~m] = 0
    real = to_sr(dn_out)
    meta = {"method": "cloudfree_composite.py (QA + thermal + temporal + speck + dark-object haze masks; tile patches; one real date per px)",
            "dates": ",".join(dates)}
    save(pre + "_sr_real.tif", real, gdal.GDT_Float32, BAND_NAMES, g, meta)
    save(pre + "_dn_real.tif", dn_out, gdal.GDT_UInt16, BAND_NAMES, g, {**meta, "scale": "SR = DN x 0.0000275 - 0.2"})
    save(pre + "_sr_harmonised.tif", disp.astype(np.float32), gdal.GDT_Float32, BAND_NAMES, g, meta)
    save(pre + "_source.tif", src.astype(np.uint8 if n < 255 else np.uint16), gdal.GDT_Byte if n < 255 else gdal.GDT_UInt16,
         ["date index"], g, meta, cats=dates)
    land = m & ~water
    for name, bands in (("truecolour", (2, 1, 0)), ("falsecolour543", (4, 3, 2))):
        rgb = np.stack([np.clip(disp[b], 0, 1) for b in bands])
        sample = rgb[:, land][:, ::7] if land.any() else rgb[:, m][:, ::7]
        top = np.percentile(sample.max(0), 99.5)
        s = np.clip(rgb / max(top, 1e-3), 0, 1) ** 0.8
        s = np.clip(s + 0.5 * (s - ndimage.gaussian_filter(s, (0, 1, 1))), 0, 1)
        save(f"{pre}_{name}.tif", np.where(m, 1 + np.round(s * 254), 0).astype(np.uint8), gdal.GDT_Byte, ["red", "green", "blue"], g, meta, rgb=True)
    gdal.Translate(pre + "_preview.png", pre + "_truecolour.tif", format="PNG", width=min(W, 2400), height=0)
    used = np.bincount(src[m], minlength=n + 1) / m.sum()
    with open(pre + "_dates.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["index", "date", "clear_share_after_masks", "share_of_composite", "bundles"])
        for k, d in enumerate(dates):
            w.writerow([k + 1, d, round(CL[k].sum() / m.sum(), 4), round(used[k + 1], 4), ";".join(i[2] for i in by_date[d])])
    txt = (f"cloud-free composite {a.prefix}: {n} dates, AOI {m.sum() * g.res ** 2 / 1e6:,.0f} km2\n"
           + "".join(f"  {dates[k]}: {used[k + 1]:.1%}\n" for k in range(n) if used[k + 1] > 0)
           + f"averaged pixels: 0 (every pixel one real date)\ngap pixels: hazy-but-clear {tier1:,}, no clear obs {tier2:,}, still empty {empty:,}\n"
           + kw_note + f"run time {time.time() - t0:.0f} s\n")
    open(pre + "_summary.txt", "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
