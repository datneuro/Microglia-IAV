"""
═══════════════════════════════════════════════════════════════════════════════
03_Colocalization-Analysis.py — Post-hoc colocalisation analysis
═══════════════════════════════════════════════════════════════════════════════
"""

import argparse
import os
import sys
import glob
import warnings
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi
from scipy.stats import pearsonr, spearmanr
from skimage.filters import threshold_otsu
from tqdm import tqdm

warnings.filterwarnings("ignore")


# ═══════════════════════════════════════════════════════════════════════════════
# CONFIG
# ═══════════════════════════════════════════════════════════════════════════════


DEFAULT_THRESHOLD_METHOD = "otsu"
PERINUCLEAR_RING_UM = 5.0
N_COSTES_RANDOMIZATIONS = 100
COSTES_BLOCK_SIZE = 16


# ═══════════════════════════════════════════════════════════════════════════════
# I/O — load channels from OME-TIFF
# ═══════════════════════════════════════════════════════════════════════════════


def normalize_to_float01(img, bit_depth):
    img = img.astype(np.float32)
    if bit_depth == 8:
        return img / 255.0
    if bit_depth == 12:
        return img / 4095.0
    if bit_depth == 16:
        return img / 65535.0
    a_min, a_max = float(img.min()), float(img.max())
    if a_max - a_min < 1e-9:
        return np.zeros_like(img, dtype=np.float32)
    return ((img - a_min) / (a_max - a_min)).astype(np.float32)


def load_channels(tiff_path, ch_dapi, ch_asyn, ch_iav, bit_depth):
    with tifffile.TiffFile(str(tiff_path)) as tf:
        data = tf.asarray()
    if data.ndim == 5:
        data = data[0]  # drop T dim
    elif data.ndim != 4:
        raise ValueError(f"Unexpected shape {data.shape}")
    return {
        "dapi": normalize_to_float01(data[ch_dapi], bit_depth),
        "asyn": normalize_to_float01(data[ch_asyn], bit_depth),
        "iav": normalize_to_float01(data[ch_iav], bit_depth),
    }


# ═══════════════════════════════════════════════════════════════════════════════
# CO-LOCALIZATION METRICS
# ═══════════════════════════════════════════════════════════════════════════════


def pearson_correlation(ch_a, ch_b, mask):
    """Pearson r within mask. NaN if mask < 30 voxels or zero variance."""
    if mask.sum() < 30:
        return np.nan, np.nan
    a = ch_a[mask].astype(np.float64)
    b = ch_b[mask].astype(np.float64)
    if np.std(a) < 1e-9 or np.std(b) < 1e-9:
        return np.nan, np.nan
    try:
        r, p = pearsonr(a, b)
        return float(r), float(p)
    except Exception:
        return np.nan, np.nan


def spearman_correlation(ch_a, ch_b, mask):
    if mask.sum() < 30:
        return np.nan, np.nan
    a = ch_a[mask].astype(np.float64)
    b = ch_b[mask].astype(np.float64)
    if np.std(a) < 1e-9 or np.std(b) < 1e-9:
        return np.nan, np.nan
    try:
        rho, p = spearmanr(a, b)
        return float(rho), float(p)
    except Exception:
        return np.nan, np.nan


def compute_threshold(arr, method="otsu"):
    """Get binarization threshold for Manders."""
    if len(arr) == 0:
        return 0.0
    if isinstance(method, (int, float)):
        return float(method)
    if method == "percentile":
        return float(np.percentile(arr, 75))
    # otsu
    try:
        if np.std(arr) < 1e-9:
            return float(np.median(arr))
        return float(threshold_otsu(arr))
    except Exception:
        return float(np.median(arr))


def manders_coefficients(ch_a, ch_b, mask, thresh_method=DEFAULT_THRESHOLD_METHOD):
    """
    M1 = (sum of ch_a voxels where ch_b > threshold_b) / (total ch_a sum)
    M2 = (sum of ch_b voxels where ch_a > threshold_a) / (total ch_b sum)
    Returns dict with M1, M2, thresholds, and % overlap area.
    """
    if mask.sum() < 30:
        return {
            "M1": np.nan,
            "M2": np.nan,
            "thresh_a": np.nan,
            "thresh_b": np.nan,
            "overlap_pct": np.nan,
        }

    a_vals = ch_a[mask].astype(np.float64)
    b_vals = ch_b[mask].astype(np.float64)

    thresh_a = compute_threshold(a_vals, thresh_method)
    thresh_b = compute_threshold(b_vals, thresh_method)

    a_above = a_vals > thresh_a
    b_above = b_vals > thresh_b

    sum_a = a_vals.sum()
    sum_b = b_vals.sum()

    # M1: fraction of A signal that co-occurs with B+
    m1 = float(a_vals[b_above].sum() / sum_a) if sum_a > 0 else np.nan
    # M2: fraction of B signal that co-occurs with A+
    m2 = float(b_vals[a_above].sum() / sum_b) if sum_b > 0 else np.nan

    overlap = float((a_above & b_above).sum() / mask.sum() * 100)

    return {
        "M1": m1,
        "M2": m2,
        "thresh_a": thresh_a,
        "thresh_b": thresh_b,
        "overlap_pct": overlap,
    }


def costes_randomization(
    ch_a,
    ch_b,
    mask,
    observed_r,
    n_permutations=N_COSTES_RANDOMIZATIONS,
    block_size=COSTES_BLOCK_SIZE,
):

    if mask.sum() < 100 or np.isnan(observed_r):
        return np.nan, np.nan, np.nan

    coords = np.where(mask)
    if len(coords[0]) == 0:
        return np.nan, np.nan, np.nan

    bbox = tuple(slice(c.min(), c.max() + 1) for c in coords)
    a_crop = ch_a[bbox].copy()
    b_crop = ch_b[bbox].copy()
    mask_crop = mask[bbox]

    if a_crop.size < 100:
        return np.nan, np.nan, np.nan

    rng = np.random.default_rng(seed=42)
    random_rs = []

    for _ in range(n_permutations):
        # Block-shuffle b_crop
        b_shuf = b_crop.copy()
        if b_shuf.ndim == 3:
            # 3D: shuffle along YX in blocks
            for z in range(b_shuf.shape[0]):
                b_shuf[z] = _block_shuffle_2d(b_shuf[z], block_size, rng)
        else:
            b_shuf = _block_shuffle_2d(b_shuf, block_size, rng)

        try:
            a_in = a_crop[mask_crop]
            b_in = b_shuf[mask_crop]
            if np.std(a_in) > 1e-9 and np.std(b_in) > 1e-9:
                r, _ = pearsonr(a_in, b_in)
                random_rs.append(r)
        except Exception:
            continue

    if len(random_rs) < 10:
        return np.nan, np.nan, np.nan

    random_rs = np.array(random_rs)
    mean_r = float(random_rs.mean())
    std_r = float(random_rs.std())

    p_value = float((random_rs >= observed_r).sum() / len(random_rs))

    return p_value, mean_r, std_r


def _block_shuffle_2d(arr_2d, block_size, rng):
    """Shuffle 2D array in block units."""
    h, w = arr_2d.shape
    n_blocks_y = h // block_size
    n_blocks_x = w // block_size
    if n_blocks_y < 2 or n_blocks_x < 2:
        # too small, fall back to pixel shuffle
        flat = arr_2d.flatten()
        rng.shuffle(flat)
        return flat.reshape(arr_2d.shape)
    # Extract blocks
    blocks = []
    for i in range(n_blocks_y):
        for j in range(n_blocks_x):
            blocks.append(
                arr_2d[
                    i * block_size : (i + 1) * block_size,
                    j * block_size : (j + 1) * block_size,
                ]
            )
    rng.shuffle(blocks)
    # Rebuild
    out = arr_2d.copy()
    idx = 0
    for i in range(n_blocks_y):
        for j in range(n_blocks_x):
            out[
                i * block_size : (i + 1) * block_size,
                j * block_size : (j + 1) * block_size,
            ] = blocks[idx]
            idx += 1
    return out


# ═══════════════════════════════════════════════════════════════════════════════
# PER-CELL ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════════
def analyze_one_image(
    row, results_dir, n_costes=N_COSTES_RANDOMIZATIONS, do_costes=True
):
    """
    For 1 image: load TIFF + nuclei mask + cytoring → per-cell coloc metrics.
    Returns DataFrame (1 row per cell) or None on failure.
    """
    image_id = row["image_id"]
    tiff_path = row["tiff_path"]
    bit_depth = int(row["bit_depth"])

    labels_dir = Path(results_dir) / "labels"
    nuc_path = labels_dir / f"{image_id}_nuclei.tif"
    cyt_path = labels_dir / f"{image_id}_cytoring.tif"

    if not nuc_path.exists():
        return None, f"nuclei mask missing: {nuc_path}"
    if not Path(tiff_path).exists():
        return None, f"TIFF missing: {tiff_path}"

    try:
        nuclei = tifffile.imread(str(nuc_path))
        cytoring = tifffile.imread(str(cyt_path)) if cyt_path.exists() else None
        channels = load_channels(
            tiff_path,
            int(row["ch_dapi_idx"]),
            int(row["ch_asyn_idx"]),
            int(row["ch_iav_idx"]),
            bit_depth,
        )
    except Exception as e:
        return None, f"load failed: {e}"

    rows = []
    nuc_ids = np.unique(nuclei)
    nuc_ids = nuc_ids[nuc_ids > 0]

    for nuc_id in nuc_ids:
        nuc_mask = nuclei == nuc_id
        cyto_mask = (cytoring == nuc_id) if cytoring is not None else None

        result = {"nucleus_id": int(nuc_id), "image_id": image_id}

        # ── aSyn × IAV (nuclear) ──
        r, p = pearson_correlation(channels["asyn"], channels["iav"], nuc_mask)
        result["nuc_asyn_iav_pearson_r"] = r
        result["nuc_asyn_iav_pearson_p"] = p

        rho, p_s = spearman_correlation(channels["asyn"], channels["iav"], nuc_mask)
        result["nuc_asyn_iav_spearman_rho"] = rho
        result["nuc_asyn_iav_spearman_p"] = p_s

        m = manders_coefficients(channels["asyn"], channels["iav"], nuc_mask)
        result["nuc_asyn_iav_M1"] = m["M1"]
        result["nuc_asyn_iav_M2"] = m["M2"]
        result["nuc_asyn_iav_overlap_pct"] = m["overlap_pct"]

        if do_costes and not np.isnan(r):
            costes_p, mean_rand_r, std_rand_r = costes_randomization(
                channels["asyn"], channels["iav"], nuc_mask, r, n_costes
            )
            result["nuc_asyn_iav_costes_p"] = costes_p
            result["nuc_asyn_iav_costes_mean_random_r"] = mean_rand_r

        # ── aSyn × DAPI chromatin density (nuclear only) ──
        r_chr, p_chr = pearson_correlation(channels["asyn"], channels["dapi"], nuc_mask)
        result["nuc_asyn_dapi_pearson_r"] = r_chr
        result["nuc_asyn_dapi_pearson_p"] = p_chr

        rho_chr, p_s_chr = spearman_correlation(
            channels["asyn"], channels["dapi"], nuc_mask
        )
        result["nuc_asyn_dapi_spearman_rho"] = rho_chr

        # ── aSyn × IAV (perinuclear) — to test redistribution ──
        if cyto_mask is not None and cyto_mask.sum() > 30:
            r_p, p_p = pearson_correlation(channels["asyn"], channels["iav"], cyto_mask)
            result["cyto_asyn_iav_pearson_r"] = r_p
            result["cyto_asyn_iav_pearson_p"] = p_p

            m_p = manders_coefficients(channels["asyn"], channels["iav"], cyto_mask)
            result["cyto_asyn_iav_M1"] = m_p["M1"]
            result["cyto_asyn_iav_M2"] = m_p["M2"]

        rows.append(result)

    if not rows:
        return None, "no cells with valid measurements"

    df = pd.DataFrame(rows)
    # Add metadata
    for col in [
        "combo",
        "antibody",
        "condition",
        "replicate",
        "source_lif",
        "magnification_class",
        "scene_subtype",
    ]:
        if col in row.index:
            df[col] = row[col]

    return df, None


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--manifest",
        type=str,
        required=True,
        help="Path to manifest_final.csv (after step 01)",
    )
    parser.add_argument(
        "--results-dir",
        type=str,
        required=True,
        help="Output folder of step 02 (contains labels/ and master_cells.csv)",
    )
    parser.add_argument(
        "--tiff-dir",
        type=str,
        default=None,
        help="Path to tiff_3d folder (auto-resolved if --rewrite-paths)",
    )
    parser.add_argument(
        "--rewrite-paths",
        action="store_true",
        help="Rewrite tiff_path in manifest using --tiff-dir basename matching",
    )
    parser.add_argument(
        "--limit", type=int, default=None, help="Process only the first N images"
    )
    parser.add_argument(
        "--n-costes",
        type=int,
        default=N_COSTES_RANDOMIZATIONS,
        help=f"Costes permutations (default {N_COSTES_RANDOMIZATIONS}; 0 = skip)",
    )
    parser.add_argument("--filter-combo", type=str, default=None)
    parser.add_argument("--filter-condition", type=str, default=None)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    print("═" * 70)
    print(f"  CO-LOCALIZATION ANALYSIS — HMC3 aSyn")
    print(f"  Started: {datetime.now()}")
    print("═" * 70)

    results_dir = Path(args.results_dir).expanduser().resolve()
    coloc_dir = results_dir / "colocalization"
    coloc_dir.mkdir(parents=True, exist_ok=True)

    # Load manifest
    manifest = pd.read_csv(args.manifest)
    print(f"\nManifest: {len(manifest)} rows")

    # Optional rewrite paths
    if args.rewrite_paths and args.tiff_dir:
        tiff_dir = Path(args.tiff_dir).expanduser().resolve()
        manifest["tiff_path"] = manifest["tiff_path"].apply(
            lambda p: (
                str(tiff_dir / os.path.basename(str(p))) if pd.notna(p) and p else ""
            )
        )

    # Filter 3D + existing
    df = manifest[manifest["scene_type"] == "3D_stack"].copy()
    df = df[df["tiff_path"].apply(lambda p: bool(p) and os.path.exists(p))]

    # Exclude PFF
    if "condition" in df.columns:
        df = df[
            ~df["condition"].astype(str).str.contains("PFF", case=False, na=False)
        ].copy()

    # Exclude high magnification
    if "magnification_class" in df.columns:
        df = df[
            ~df["magnification_class"]
            .astype(str)
            .str.contains("high", case=False, na=False)
        ].copy()

    if args.filter_combo:
        df = df[df["combo"] == args.filter_combo]
    if args.filter_condition:
        df = df[df["condition"] == args.filter_condition]

    # Drop missing metadata
    required = ["combo", "antibody", "condition", "replicate"]
    df = df.dropna(subset=required)
    df = df[(df[required] != "").all(axis=1)]

    if args.limit:
        df = df.head(args.limit)

    print(f"Processing {len(df)} images")
    print(
        f"Costes randomizations: {args.n_costes} per cell "
        f"({'DISABLED' if args.n_costes == 0 else 'enabled'})"
    )

    do_costes = args.n_costes > 0

    # Process
    all_results = []
    n_done = n_skipped = n_failed = 0

    for _, row in tqdm(df.iterrows(), total=len(df), desc="Co-loc"):
        image_id = row["image_id"]
        out_path = coloc_dir / f"{image_id}_coloc.csv"

        if out_path.exists() and not args.overwrite:
            try:
                all_results.append(pd.read_csv(out_path))
                n_skipped += 1
            except Exception:
                pass
            continue

        result_df, err = analyze_one_image(
            row, results_dir, n_costes=args.n_costes, do_costes=do_costes
        )
        if result_df is None:
            print(f"\n[FAIL] {image_id}: {err}")
            n_failed += 1
            continue

        result_df.to_csv(out_path, index=False)
        all_results.append(result_df)
        n_done += 1

    print(f"\n\n--- Summary ---")
    print(f"  Processed: {n_done}")
    print(f"  Skipped:   {n_skipped} (already done)")
    print(f"  Failed:    {n_failed}")

    # Master CSV
    if all_results:
        master = pd.concat(all_results, ignore_index=True)
        master_path = results_dir / "master_coloc.csv"
        master.to_csv(master_path, index=False)
        print(f"\nmaster_coloc.csv: {len(master)} cells → {master_path}")

        # Merge into master_cells.csv → master_cells_coloc.csv
        cells_path = results_dir / "master_cells.csv"
        if cells_path.exists():
            cells = pd.read_csv(cells_path)
            merged = cells.merge(
                master.drop(
                    columns=[
                        c
                        for c in master.columns
                        if c in cells.columns and c not in ["nucleus_id", "image_id"]
                    ]
                ),
                on=["nucleus_id", "image_id"],
                how="left",
            )
            merged_path = results_dir / "master_cells_coloc.csv"
            merged.to_csv(merged_path, index=False)
            print(
                f"master_cells_coloc.csv (cells + coloc): {len(merged)} → {merged_path}"
            )

    print("\n" + "═" * 70)
    print("  DONE")
    print("═" * 70)


if __name__ == "__main__":
    main()
