#!/usr/bin/env python3
"""
═══════════════════════════════════════════════════════════════════════════════
02_Cellpose-Segmentation-Quantification.py
═══════════════════════════════════════════════════════════════════════════════
"""

import os
import sys
import time
import gc
import glob
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from tqdm import tqdm

import matplotlib

matplotlib.use("Agg")  # No display needed
import matplotlib.pyplot as plt

from scipy import ndimage as ndi
from scipy.ndimage import gaussian_filter
from scipy.stats import pearsonr as _pearsonr, spearmanr as _spearmanr

from skimage import morphology
from skimage.filters import threshold_otsu, gaussian
from skimage.morphology import remove_small_objects, ball
from skimage.segmentation import watershed, clear_border, find_boundaries
from skimage.feature import blob_log, peak_local_max

from cellpose import models as cp_models
from csbdeep.utils import normalize as cd_normalize

warnings.filterwarnings("ignore")


# ═══════════════════════════════════════════════════════════════════════════════
# CONFIG
# ═══════════════════════════════════════════════════════════════════════════════

# ─── PATHS ───
BASE_DIR = "Confocal"
MANIFEST_FILENAME = "manifest_final.csv"
TIFF_SUBFOLDER = "tiff_3d"

DATA_DIR = BASE_DIR
OUTPUT_DIR = os.path.join(BASE_DIR, "results_segmentation")
TMP_DIR = os.path.join(BASE_DIR, "temp_processing")

MANIFEST_PATH = os.path.join(DATA_DIR, MANIFEST_FILENAME)
TIFF_DIR = os.path.join(DATA_DIR, TIFF_SUBFOLDER)

# ─── PROCESSING ───
LIMIT = None
OVERWRITE_EXISTING = True
USE_WATERSHED_ONLY = False

# ─── CONDITION RENAMING ───
CONDITION_RENAME = {
    "A/PR/8/34": "V1",
    "A/NWS/33": "V2",
    "PFF+A/PR/8/34": "PFF+V1",
    "PFF+A/NWS/33": "PFF+V2",
    "PFF": "PFF",
    "Control": "Control",
}

# ─── A-PRIORI DESIGN EXCLUSIONS ───
EXCLUDE_PFF = True
EXCLUDE_HIGH_MAG = True

# ─── FILTERS (optional) ───
FILTER_MAG = None
FILTER_COMBO = None
FILTER_CONDITION = None

# ─── SEGMENTATION (Cellpose) ───
CELLPOSE_DIAMETER_UM = 15.0
CELLPOSE_STITCH_THRESHOLD = 0.25
CELLPOSE_FLOW_THRESHOLD = 0.4
CELLPOSE_MIN_SIZE = 15

# ─── PREPROCESSING ───
ROLLING_BALL_RADIUS_UM = 5.0
TOPHAT_RADIUS_UM = 2.0
CYTORING_WIDTH_UM = 5.0
DAPI_SMOOTH_SIGMA_PX = 1.5
IAV_BGSUB_SIGMA_UM = 5.0

# ─── PUNCTA DETECTION ───
PUNCTA_MIN_SIGMA = 1.0
PUNCTA_MAX_SIGMA = 3.0
PUNCTA_NUM_SIGMA = 4
PUNCTA_THRESHOLD = 0.20
PUNCTA_OVERLAP = 0.5

# ─── IAV CLASSIFICATION (dual threshold) ───
N_CONTROL_FOVS_FOR_THRESHOLD = 3
IAV_CELL_SHELL_UM = 3.0

# ─── CO-LOCALIZATION ───
COLOC_ENABLE = True
COLOC_DO_COSTES = False
COLOC_N_COSTES = 100
COLOC_BLOCK_SIZE = 16
COLOC_THRESHOLD_METHOD = "otsu"

# ─── QC OVERLAY ───
QC_DPI = 120
QC_MIP_AXIS = 0


# ─── SETUP & VALIDATE ───
if not os.path.exists(MANIFEST_PATH):
    raise FileNotFoundError(f"Manifest not found: {MANIFEST_PATH}")
if not os.path.exists(TIFF_DIR):
    raise FileNotFoundError(f"TIFF dir not found: {TIFF_DIR}")

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(TMP_DIR, exist_ok=True)
for sub in ["labels", "qc_overlays", "per_image_csv"]:
    os.makedirs(os.path.join(TMP_DIR, sub), exist_ok=True)
    os.makedirs(os.path.join(OUTPUT_DIR, sub), exist_ok=True)

OUT_DIRS = {
    "labels": os.path.join(OUTPUT_DIR, "labels"),
    "qc": os.path.join(OUTPUT_DIR, "qc_overlays"),
    "per_image_csv": os.path.join(OUTPUT_DIR, "per_image_csv"),
}


# ═══════════════════════════════════════════════════════════════════════════════
# LOAD & FILTER MANIFEST
# ═══════════════════════════════════════════════════════════════════════════════


def load_and_filter_manifest():
    print("=" * 70)
    print(f"Loading manifest: {MANIFEST_PATH}")
    manifest_raw = pd.read_csv(MANIFEST_PATH)
    print(f"  {len(manifest_raw)} rows total")

    # Point tiff_path to the local TIFF folder
    def rewrite_path_local(p):
        if pd.isna(p) or not p:
            return ""
        return os.path.join(TIFF_DIR, os.path.basename(str(p)))

    manifest_raw["tiff_path"] = manifest_raw["tiff_path"].apply(rewrite_path_local)

    # Rename conditions
    if "condition" in manifest_raw.columns:
        manifest_raw["condition"] = manifest_raw["condition"].replace(CONDITION_RENAME)

    # Keep 3D only + existing files
    df = manifest_raw[manifest_raw["scene_type"] == "3D_stack"].copy()
    df = df[df["tiff_path"].apply(lambda p: bool(p) and os.path.exists(p))]
    print(f"  3D scenes with TIFF on disk: {len(df)}")

    # A-priori design exclusions
    n_before = len(df)
    if EXCLUDE_HIGH_MAG and "magnification_class" in df.columns:
        df = df[df["magnification_class"] != "high_mag"].copy()
    if EXCLUDE_PFF and "condition" in df.columns:
        df = df[~df["condition"].astype(str).str.contains("PFF", na=False)].copy()
    print(f"  After design exclusions (PFF + high_mag): {n_before} → {len(df)}")

    # Optional user filters
    if FILTER_MAG and "magnification_class" in df.columns:
        df = df[df["magnification_class"] == FILTER_MAG]
    if FILTER_COMBO and "combo" in df.columns:
        df = df[df["combo"] == FILTER_COMBO]
    if FILTER_CONDITION and "condition" in df.columns:
        df = df[df["condition"] == FILTER_CONDITION]

    # Required metadata
    required = ["combo", "antibody", "condition", "replicate"]
    df = df.dropna(subset=required)
    df = df[(df[required] != "").all(axis=1)]

    if LIMIT:
        df = df.head(LIMIT)

    print(f"\nFinal: {len(df)} FOVs to process")
    print("\nCombo × Condition breakdown:")
    print(df.groupby(["combo", "condition"]).size().unstack(fill_value=0))
    print("=" * 70)
    return df


# ═══════════════════════════════════════════════════════════════════════════════
# CORE FUNCTIONS (channel I/O + preprocessing)
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


def load_ome_tiff(tiff_path, ch_dapi, ch_asyn, ch_iav, bit_depth):
    with tifffile.TiffFile(str(tiff_path)) as tf:
        data = tf.asarray()
    if data.ndim == 5:
        data = data[0]
    elif data.ndim != 4:
        raise ValueError(f"Unexpected shape {data.shape}")
    result = {}
    for ch_name, ch_idx in [("dapi", ch_dapi), ("asyn", ch_asyn), ("iav", ch_iav)]:
        if ch_idx < 0 or ch_idx >= data.shape[0]:
            raise ValueError(f"Channel idx {ch_idx} out of range for {ch_name}")
        result[ch_name] = normalize_to_float01(data[ch_idx], bit_depth)
    return result


def load_iav_only(row):
    return load_ome_tiff(
        row["tiff_path"],
        int(row["ch_dapi_idx"]),
        int(row["ch_asyn_idx"]),
        int(row["ch_iav_idx"]),
        int(row["bit_depth"]),
    )


def create_cytoplasm_ring(nuclei_labels, ring_width_um, voxel_size_um):
    if nuclei_labels.max() == 0:
        return np.zeros_like(nuclei_labels, dtype=np.uint16)
    bg_mask = nuclei_labels == 0
    distances, indices = ndi.distance_transform_edt(
        bg_mask, sampling=voxel_size_um, return_indices=True
    )
    nearest_labels = nuclei_labels[tuple(indices)]
    ring_mask = bg_mask & (distances <= ring_width_um) & (distances > 0)
    ring_labels = np.zeros_like(nuclei_labels, dtype=np.uint16)
    ring_labels[ring_mask] = nearest_labels[ring_mask]
    return ring_labels


def gaussian_smooth_dapi(dapi, sigma_px):
    if sigma_px <= 0:
        return dapi.astype(np.float32)
    sigmas = (
        (0.0, float(sigma_px), float(sigma_px))
        if dapi.ndim == 3
        else (float(sigma_px), float(sigma_px))
    )
    return gaussian_filter(dapi.astype(np.float32), sigma=sigmas, mode="reflect")


def rolling_ball_bg_subtract(img, radius_um, voxel_size_um):
    """Gaussian high-pass: FFT-fast, equivalent rolling-ball for sparse cells."""
    sigma_z = radius_um / voxel_size_um[0] if voxel_size_um[0] > 0 else 0
    sigma_y = radius_um / voxel_size_um[1] if voxel_size_um[1] > 0 else 0
    sigma_x = radius_um / voxel_size_um[2] if voxel_size_um[2] > 0 else 0
    sigmas = (sigma_z, sigma_y, sigma_x) if img.ndim == 3 else (sigma_y, sigma_x)
    bg = gaussian_filter(img.astype(np.float32), sigma=sigmas, mode="reflect")
    out = img - bg
    out[out < 0] = 0
    return out.astype(np.float32)


def white_tophat_3d(img, radius_um, voxel_size_um):
    """Gaussian high-pass with a smaller radius, used for puncta enhancement."""
    return rolling_ball_bg_subtract(img, radius_um, voxel_size_um)


# ═══════════════════════════════════════════════════════════════════════════════
# SEGMENTATION (Cellpose-SAM)
# ═══════════════════════════════════════════════════════════════════════════════

_cellpose_model = None


def get_cellpose_model():
    """Load the Cellpose-SAM model once and cache it."""
    global _cellpose_model
    if _cellpose_model is None:
        import torch

        device = "mps" if torch.backends.mps.is_available() else "cpu"
        print(f"Loading Cellpose-SAM on device={device}...")
        _cellpose_model = cp_models.CellposeModel(gpu=(device != "cpu"))
    return _cellpose_model


def segment_nuclei_cellpose(dapi_data, voxel_size_um):
    """Cellpose-SAM nuclear segmentation: 2D masks per plane stitched into 3D objects."""
    model = get_cellpose_model()
    diam_px = CELLPOSE_DIAMETER_UM / voxel_size_um[1]
    dapi_norm = cd_normalize(dapi_data, 1, 99.8)
    masks, _, _ = model.eval(
        dapi_norm,
        diameter=diam_px,
        do_3D=False,
        stitch_threshold=CELLPOSE_STITCH_THRESHOLD,
        flow_threshold=CELLPOSE_FLOW_THRESHOLD,
        min_size=CELLPOSE_MIN_SIZE,
        z_axis=0,  # stacks are (Z, Y, X)
    )
    if masks is None:
        return np.zeros_like(dapi_data, dtype=np.uint16)

    # Remove objects touching the lateral (X/Y) borders; objects reaching the
    # first/last Z plane are kept
    border_ids = set()
    for z in range(masks.shape[0]):
        sl = masks[z]
        edge_vals = np.concatenate([sl[0, :], sl[-1, :], sl[:, 0], sl[:, -1]])
        for nid in np.unique(edge_vals):
            if nid > 0:
                border_ids.add(int(nid))

    if border_ids:
        mask_remove = np.isin(masks, list(border_ids))
        masks[mask_remove] = 0

    masks = remove_small_objects(
        masks.astype(np.int32), min_size=CELLPOSE_MIN_SIZE * 10
    ).astype(np.uint16)
    return masks


def segment_nuclei_watershed_fallback(dapi_data, voxel_size_um, min_volume_um3=50):
    """Distance-transform watershed, used only if Cellpose segmentation fails."""
    dapi_smooth = gaussian(dapi_data, sigma=1.5)
    try:
        thresh = threshold_otsu(dapi_smooth)
    except Exception:
        thresh = float(dapi_smooth.mean() + dapi_smooth.std())
    mask = dapi_smooth > thresh
    distance = ndi.distance_transform_edt(mask, sampling=voxel_size_um)
    local_max = peak_local_max(distance, min_distance=5, labels=mask)
    markers = np.zeros_like(mask, dtype=np.int32)
    for i, (z, y, x) in enumerate(local_max, start=1):
        markers[z, y, x] = i
    labels = watershed(-distance, markers, mask=mask)
    voxel_vol = np.prod(voxel_size_um)
    min_voxels = int(min_volume_um3 / voxel_vol)
    labels = remove_small_objects(labels.astype(np.int32), min_size=min_voxels)
    labels = clear_border(labels)
    return labels.astype(np.uint16)


# ═══════════════════════════════════════════════════════════════════════════════
# DUAL THRESHOLD COMPUTATION (intensity + area, derived from Control fields)
# ═══════════════════════════════════════════════════════════════════════════════


def compute_objective_dual_thresholds(manifest_df, n_control_fovs=3):
    combo_thresholds = {}
    print("=" * 70)
    print("Computing dual IAV thresholds (intensity p99 + area mean+3σ)")
    print("=" * 70)

    for combo in manifest_df["combo"].dropna().unique():
        controls = manifest_df[
            (manifest_df["combo"] == combo) & (manifest_df["condition"] == "Control")
        ].head(n_control_fovs)

        if len(controls) == 0:
            print(f"⚠️ {combo}: no Control → fallback (intensity=0.15, area=50)")
            combo_thresholds[combo] = {"intensity": 0.15, "area": 50}
            continue

        all_control_pixels = []
        control_cells_data = []

        for _, row in controls.iterrows():
            try:
                label_path = os.path.join(
                    OUT_DIRS["labels"], f"{row['image_id']}_nuclei.tif"
                )
                ch = load_iav_only(row)
                iav_raw = ch["iav"]
                voxel_size = (
                    float(row["voxel_z_um"]),
                    float(row["voxel_y_um"]),
                    float(row["voxel_x_um"]),
                )

                if os.path.exists(label_path):
                    nuclei = tifffile.imread(label_path)
                else:
                    print(
                        f"   {row['image_id']}: segmenting for threshold computation..."
                    )
                    dapi_smooth = gaussian_smooth_dapi(ch["dapi"], DAPI_SMOOTH_SIGMA_PX)
                    nuclei = segment_nuclei_cellpose(dapi_smooth, voxel_size)
                    tifffile.imwrite(label_path, nuclei, compression="zlib")

                perinuc = create_cytoplasm_ring(nuclei, IAV_CELL_SHELL_UM, voxel_size)
                full_cell_mask = (nuclei > 0) | (perinuc > 0)
                if full_cell_mask.sum() > 0:
                    all_control_pixels.append(iav_raw[full_cell_mask])
                    control_cells_data.append(
                        {"iav": iav_raw, "nuc": nuclei, "cyto": perinuc}
                    )
            except Exception as e:
                print(f"   {row['image_id']}: load failed ({e})")
                continue

        if not all_control_pixels:
            combo_thresholds[combo] = {"intensity": 0.15, "area": 50}
            continue

        # INTENSITY threshold: p99 of Control cell pixels
        merged_pixels = np.concatenate(all_control_pixels)
        intensity_thresh = float(np.percentile(merged_pixels, 99))

        # AREA threshold: collect bright counts per Control cell
        control_noise_areas = []
        for data in control_cells_data:
            iav_raw, nuclei, cyto_ring = data["iav"], data["nuc"], data["cyto"]
            bright_mask = iav_raw > intensity_thresh
            unique_nids = np.unique(nuclei)
            for nuc_id in unique_nids[unique_nids > 0]:
                cell_mask = (nuclei == nuc_id) | (cyto_ring == nuc_id)
                if cell_mask.sum() < 10:
                    continue
                noise_area = int((bright_mask & cell_mask).sum())
                control_noise_areas.append(noise_area)

        if control_noise_areas:
            areas_arr = np.array(control_noise_areas)
            p95_cutoff = np.percentile(areas_arr, 95)
            trimmed = areas_arr[areas_arr <= p95_cutoff]
            mean_noise = float(trimmed.mean())
            std_noise = (
                float(trimmed.std()) if len(trimmed) > 1 else max(mean_noise * 0.2, 5)
            )
            area_thresh = int(mean_noise + 3 * std_noise)
            print(
                f"   {combo}: Control n={len(areas_arr)} cells, "
                f"after trimming >p95 n={len(trimmed)}, "
                f"mean={mean_noise:.1f}, std={std_noise:.1f}"
            )
        else:
            area_thresh = 50

        if combo == "Syn1_16bit":
            old_area = area_thresh
            area_thresh = max(area_thresh, 7000)
            if area_thresh != old_area:
                print(
                    f"   ⚠️ {combo}: area_thresh raised from {old_area} → {area_thresh} (min floor for 16-bit)"
                )
        elif combo == "Abcam_EPR20535_8bit":
            old_area = area_thresh
            area_thresh = max(area_thresh, 9000)
            if area_thresh != old_area:
                print(
                    f"   ⚠️ {combo}: area_thresh raised from {old_area} → {area_thresh} (min floor for 8-bit)"
                )

        combo_thresholds[combo] = {"intensity": intensity_thresh, "area": area_thresh}
        print(f"✅ {combo} → Intensity={intensity_thresh:.4f} | Area={area_thresh} px")

    return combo_thresholds


# ═══════════════════════════════════════════════════════════════════════════════
# QUANTIFICATION FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════


def classify_infected_cells(
    nuclei_labels, iav_data, voxel_size_um, method=None, global_threshold=None
):
    perinuc_ring = create_cytoplasm_ring(
        nuclei_labels, IAV_CELL_SHELL_UM, voxel_size_um
    )

    if isinstance(global_threshold, dict) and "intensity" in global_threshold:
        intensity_thresh = float(global_threshold["intensity"])
        area_thresh = int(global_threshold["area"])
        bright_mask = iav_data > intensity_thresh

        info = {}
        for nuc_id in np.unique(nuclei_labels):
            if nuc_id == 0:
                continue
            cell_mask = (nuclei_labels == nuc_id) | (perinuc_ring == nuc_id)
            cell_size = int(cell_mask.sum())
            if cell_size < 10:
                continue
            bright_area = int((bright_mask & cell_mask).sum())
            iav_mean = float(iav_data[cell_mask].mean())
            info[int(nuc_id)] = {
                "iav_mean": iav_mean,
                "iav_bright_area": bright_area,
                "cell_size_voxels": cell_size,
                "infected": bool(bright_area >= area_thresh),
            }
        return info, intensity_thresh

    # Otsu fallback if no dual threshold is available
    means, nuc_ids = [], []
    for nuc_id in np.unique(nuclei_labels):
        if nuc_id == 0:
            continue
        cell_mask = (nuclei_labels == nuc_id) | (perinuc_ring == nuc_id)
        if cell_mask.sum() == 0:
            continue
        means.append(float(iav_data[cell_mask].mean()))
        nuc_ids.append(int(nuc_id))
    if len(means) == 0:
        return {}, 0.0
    means_arr = np.array(means)
    if isinstance(global_threshold, (int, float)):
        thresh = float(global_threshold)
    else:
        try:
            thresh = float(threshold_otsu(means_arr))
        except Exception:
            thresh = float(means_arr.mean() + means_arr.std())
    return (
        {
            nid: {"iav_mean": m, "infected": m > thresh}
            for nid, m in zip(nuc_ids, means)
        },
        thresh,
    )


def measure_nc_ratio(asyn_data, nuclei_labels, cyto_ring, infection_info=None):
    rows = []
    for nuc_id in np.unique(nuclei_labels):
        if nuc_id == 0:
            continue
        nuc_mask = nuclei_labels == nuc_id
        cyto_mask = cyto_ring == nuc_id
        if nuc_mask.sum() == 0 or cyto_mask.sum() == 0:
            continue
        nuc_vals = asyn_data[nuc_mask]
        cyto_vals = asyn_data[cyto_mask]
        nuc_mean = float(nuc_vals.mean())
        nuc_std = float(nuc_vals.std())
        cyto_mean = float(cyto_vals.mean())
        cyto_std = float(cyto_vals.std())
        nc_ratio = nuc_mean / (cyto_mean + 1e-9)
        row = {
            "nucleus_id": int(nuc_id),
            "nucleus_volume_voxels": int(nuc_mask.sum()),
            "cytoplasm_volume_voxels": int(cyto_mask.sum()),
            "asyn_nuclear_mean": nuc_mean,
            "asyn_nuclear_sum": float(nuc_vals.sum()),
            "asyn_cytoplasmic_mean": cyto_mean,
            "asyn_cytoplasmic_sum": float(cyto_vals.sum()),
            "asyn_nuclear_std": nuc_std,
            "asyn_nuclear_CV": float(nuc_std / (nuc_mean + 1e-9)),
            "asyn_nuclear_p90": float(np.percentile(nuc_vals, 90)),
            "asyn_cytoplasmic_std": cyto_std,
            "asyn_cytoplasmic_CV": float(cyto_std / (cyto_mean + 1e-9)),
            "asyn_cytoplasmic_p90": float(np.percentile(cyto_vals, 90)),
            "NC_ratio": nc_ratio,
            "log2_NC_ratio": float(np.log2(nc_ratio + 1e-9)),
        }
        if infection_info and nuc_id in infection_info:
            row["iav_mean"] = infection_info[nuc_id]["iav_mean"]
            row["iav_bright_area"] = infection_info[nuc_id].get("iav_bright_area", 0)
            row["cell_size_voxels"] = infection_info[nuc_id].get("cell_size_voxels", 0)
            row["infected"] = bool(infection_info[nuc_id]["infected"])
        rows.append(row)
    return pd.DataFrame(rows)


def detect_nuclear_puncta(asyn_tophat, nuclei_labels, **kwargs):
    """Laplacian-of-Gaussian blob detection on the high-pass aSyn channel; keeps blobs inside nuclei."""
    if asyn_tophat.max() < 1e-9:
        return pd.DataFrame()
    a_min, a_max = float(asyn_tophat.min()), float(asyn_tophat.max())
    if a_max - a_min < 1e-9:
        return pd.DataFrame()
    norm = (asyn_tophat - a_min) / (a_max - a_min)
    blobs = blob_log(
        norm,
        min_sigma=PUNCTA_MIN_SIGMA,
        max_sigma=PUNCTA_MAX_SIGMA,
        num_sigma=PUNCTA_NUM_SIGMA,
        threshold=PUNCTA_THRESHOLD,
        overlap=PUNCTA_OVERLAP,
    )
    if len(blobs) == 0:
        return pd.DataFrame()
    rows = []
    n_outside = n_wrong = 0
    for blob in blobs:
        if asyn_tophat.ndim == 3:
            z, y, x, sigma = blob
            z, y, x = int(z), int(y), int(x)
            if not (
                0 <= z < nuclei_labels.shape[0]
                and 0 <= y < nuclei_labels.shape[1]
                and 0 <= x < nuclei_labels.shape[2]
            ):
                continue
            nuc_id = int(nuclei_labels[z, y, x])
        else:
            y, x, sigma = blob
            y, x = int(y), int(x)
            if not (
                0 <= y < nuclei_labels.shape[0] and 0 <= x < nuclei_labels.shape[1]
            ):
                continue
            z = 0
            nuc_id = int(nuclei_labels[y, x])
        if nuc_id == 0:
            n_outside += 1
            continue
        intensity = float(
            asyn_tophat[z, y, x] if asyn_tophat.ndim == 3 else asyn_tophat[y, x]
        )
        rows.append(
            {
                "puncta_z": z,
                "puncta_y": y,
                "puncta_x": x,
                "puncta_sigma": float(sigma),
                "puncta_intensity": intensity,
                "nucleus_id": nuc_id,
            }
        )
    print(f"    Puncta: {len(rows)} (outside: {n_outside}, wrong nuc: {n_wrong})")
    return pd.DataFrame(rows)


def aggregate_puncta_per_cell(puncta_df, nc_df):
    """Aggregate puncta count + intensity per cell."""
    if "nucleus_id" not in nc_df.columns:
        return nc_df
    nc_df["n_nuclear_puncta"] = 0
    nc_df["mean_puncta_intensity"] = 0.0
    nc_df["max_puncta_intensity"] = 0.0
    if len(puncta_df) == 0:
        return nc_df
    agg = (
        puncta_df.groupby("nucleus_id")
        .agg(
            n=("puncta_intensity", "count"),
            mean_i=("puncta_intensity", "mean"),
            max_i=("puncta_intensity", "max"),
        )
        .reset_index()
    )
    nc_df = nc_df.merge(agg, on="nucleus_id", how="left")
    nc_df["n_nuclear_puncta"] = nc_df.pop("n").fillna(0).astype(int)
    nc_df["mean_puncta_intensity"] = nc_df.pop("mean_i").fillna(0)
    nc_df["max_puncta_intensity"] = nc_df.pop("max_i").fillna(0)
    return nc_df


def measure_puncta_radial_position(puncta_df, nuclei_labels, voxel_size_um):
    """Distance (µm) of each punctum from its nucleus centroid."""
    if len(puncta_df) == 0:
        return puncta_df
    centroids = ndi.center_of_mass(
        np.ones_like(nuclei_labels),
        labels=nuclei_labels,
        index=np.unique(nuclei_labels)[1:],
    )
    centroids_dict = {
        int(nid): c for nid, c in zip(np.unique(nuclei_labels)[1:], centroids)
    }
    distances = []
    for _, row in puncta_df.iterrows():
        nid = int(row["nucleus_id"])
        if nid not in centroids_dict:
            distances.append(np.nan)
            continue
        cz, cy, cx = centroids_dict[nid]
        dz = (row["puncta_z"] - cz) * voxel_size_um[0]
        dy = (row["puncta_y"] - cy) * voxel_size_um[1]
        dx = (row["puncta_x"] - cx) * voxel_size_um[2]
        distances.append(float(np.sqrt(dz**2 + dy**2 + dx**2)))
    puncta_df = puncta_df.copy()
    puncta_df["radial_dist_um"] = distances
    return puncta_df


def correlate_asyn_with_chromatin(asyn_data, dapi_data, nuclei_labels):
    """Per-cell Pearson correlation aSyn × DAPI (chromatin density proxy)."""
    rows = []
    valid_rs = []
    for nuc_id in np.unique(nuclei_labels):
        if nuc_id == 0:
            continue
        mask = nuclei_labels == nuc_id
        if mask.sum() < 30:
            continue
        a_vals = asyn_data[mask].astype(np.float64)
        d_vals = dapi_data[mask].astype(np.float64)
        if np.std(a_vals) < 1e-9 or np.std(d_vals) < 1e-9:
            continue
        try:
            r, p = _pearsonr(a_vals, d_vals)
            rows.append(
                {
                    "nucleus_id": int(nuc_id),
                    "asyn_dapi_pearson_r": float(r),
                    "asyn_dapi_pearson_p": float(p),
                }
            )
            valid_rs.append(float(r))
        except Exception:
            continue
    df = pd.DataFrame(rows)
    stats = {
        "pearson_r": float(np.mean(valid_rs)) if valid_rs else np.nan,
        "n_cells": len(valid_rs),
    }
    return df, stats


# ─── Coloc functions ───
def coloc_pearson(ch_a, ch_b, mask):
    if mask.sum() < 30:
        return np.nan, np.nan
    a = ch_a[mask].astype(np.float64)
    b = ch_b[mask].astype(np.float64)
    if np.std(a) < 1e-9 or np.std(b) < 1e-9:
        return np.nan, np.nan
    try:
        r, p = _pearsonr(a, b)
        return float(r), float(p)
    except Exception:
        return np.nan, np.nan


def coloc_spearman(ch_a, ch_b, mask):
    if mask.sum() < 30:
        return np.nan, np.nan
    a = ch_a[mask].astype(np.float64)
    b = ch_b[mask].astype(np.float64)
    if np.std(a) < 1e-9 or np.std(b) < 1e-9:
        return np.nan, np.nan
    try:
        rho, p = _spearmanr(a, b)
        return float(rho), float(p)
    except Exception:
        return np.nan, np.nan


def coloc_threshold(arr, method):
    if len(arr) == 0:
        return 0.0
    if isinstance(method, (int, float)):
        return float(method)
    if method == "percentile":
        return float(np.percentile(arr, 75))
    try:
        if np.std(arr) < 1e-9:
            return float(np.median(arr))
        return float(threshold_otsu(arr))
    except Exception:
        return float(np.median(arr))


def coloc_manders(ch_a, ch_b, mask, thresh_method=COLOC_THRESHOLD_METHOD):
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
    thresh_a = coloc_threshold(a_vals, thresh_method)
    thresh_b = coloc_threshold(b_vals, thresh_method)
    a_above = a_vals > thresh_a
    b_above = b_vals > thresh_b
    sum_a, sum_b = a_vals.sum(), b_vals.sum()
    m1 = float(a_vals[b_above].sum() / sum_a) if sum_a > 0 else np.nan
    m2 = float(b_vals[a_above].sum() / sum_b) if sum_b > 0 else np.nan
    overlap = float((a_above & b_above).sum() / mask.sum() * 100)
    return {
        "M1": m1,
        "M2": m2,
        "thresh_a": thresh_a,
        "thresh_b": thresh_b,
        "overlap_pct": overlap,
    }


def _block_shuffle_2d(arr_2d, block_size, rng):
    h, w = arr_2d.shape
    nh, nw = h // block_size, w // block_size
    if nh == 0 or nw == 0:
        flat = arr_2d.flatten()
        rng.shuffle(flat)
        return flat.reshape(arr_2d.shape)
    blocks = []
    for i in range(nh):
        for j in range(nw):
            blocks.append(
                arr_2d[
                    i * block_size : (i + 1) * block_size,
                    j * block_size : (j + 1) * block_size,
                ].copy()
            )
    idx = list(range(len(blocks)))
    rng.shuffle(idx)
    blocks = [blocks[i] for i in idx]
    out = arr_2d.copy()
    k = 0
    for i in range(nh):
        for j in range(nw):
            out[
                i * block_size : (i + 1) * block_size,
                j * block_size : (j + 1) * block_size,
            ] = blocks[k]
            k += 1
    return out


def coloc_costes(ch_a, ch_b, mask, observed_r, n_permutations=None, block_size=None):
    if n_permutations is None:
        n_permutations = COLOC_N_COSTES
    if block_size is None:
        block_size = COLOC_BLOCK_SIZE
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
        b_shuf = b_crop.copy()
        if b_shuf.ndim == 3:
            for z in range(b_shuf.shape[0]):
                b_shuf[z] = _block_shuffle_2d(b_shuf[z], block_size, rng)
        else:
            b_shuf = _block_shuffle_2d(b_shuf, block_size, rng)
        try:
            a_in = a_crop[mask_crop]
            b_in = b_shuf[mask_crop]
            if np.std(a_in) > 1e-9 and np.std(b_in) > 1e-9:
                r, _ = _pearsonr(a_in, b_in)
                random_rs.append(r)
        except Exception:
            continue
    if len(random_rs) < 10:
        return np.nan, np.nan, np.nan
    random_rs = np.array(random_rs)
    p_value = float((random_rs >= observed_r).sum() / len(random_rs))
    return p_value, float(random_rs.mean()), float(random_rs.std())


def compute_coloc_for_image(channels, nuclei_labels, cyto_ring, do_costes=True):
    """Per-cell coloc metrics: aSyn×IAV (nuc + perinuc), aSyn×DAPI (nuc)."""
    rows = []
    nuc_ids = np.unique(nuclei_labels)
    nuc_ids = nuc_ids[nuc_ids > 0]
    for nuc_id in nuc_ids:
        nuc_mask = nuclei_labels == int(nuc_id)
        cyto_mask = cyto_ring == int(nuc_id) if cyto_ring is not None else None
        result = {"nucleus_id": int(nuc_id)}

        r, p = coloc_pearson(channels["asyn"], channels["iav"], nuc_mask)
        result["nuc_asyn_iav_pearson_r"] = r
        result["nuc_asyn_iav_pearson_p"] = p
        rho, p_s = coloc_spearman(channels["asyn"], channels["iav"], nuc_mask)
        result["nuc_asyn_iav_spearman_rho"] = rho
        result["nuc_asyn_iav_spearman_p"] = p_s
        m = coloc_manders(channels["asyn"], channels["iav"], nuc_mask)
        result["nuc_asyn_iav_M1"] = m["M1"]
        result["nuc_asyn_iav_M2"] = m["M2"]
        result["nuc_asyn_iav_overlap_pct"] = m["overlap_pct"]
        if do_costes and not np.isnan(r):
            costes_p, mean_rr, std_rr = coloc_costes(
                channels["asyn"], channels["iav"], nuc_mask, r
            )
            result["nuc_asyn_iav_costes_p"] = costes_p
            result["nuc_asyn_iav_costes_mean_random_r"] = mean_rr

        r_chr, p_chr = coloc_pearson(channels["asyn"], channels["dapi"], nuc_mask)
        result["nuc_asyn_dapi_pearson_r"] = r_chr
        result["nuc_asyn_dapi_pearson_p"] = p_chr
        rho_chr, _ = coloc_spearman(channels["asyn"], channels["dapi"], nuc_mask)
        result["nuc_asyn_dapi_spearman_rho"] = rho_chr

        if cyto_mask is not None and cyto_mask.sum() > 30:
            r_p, p_p = coloc_pearson(channels["asyn"], channels["iav"], cyto_mask)
            result["cyto_asyn_iav_pearson_r"] = r_p
            result["cyto_asyn_iav_pearson_p"] = p_p
            m_p = coloc_manders(channels["asyn"], channels["iav"], cyto_mask)
            result["cyto_asyn_iav_M1"] = m_p["M1"]
            result["cyto_asyn_iav_M2"] = m_p["M2"]
        rows.append(result)
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def make_qc_overlay(
    channels,
    nuclei_labels,
    cyto_ring,
    puncta_df,
    out_path,
    image_id,
    info_dict,
    voxel_size,
):
    """4-panel QC overlay (MIP)."""

    def mip(arr):
        return arr.max(axis=0) if arr.ndim == 3 else arr

    dapi_mip = mip(channels["dapi"])
    asyn_mip = mip(channels["asyn"])
    iav_mip = mip(channels["iav"])
    nuc_mip = mip(nuclei_labels)
    boundaries = find_boundaries(nuc_mip, mode="outer")

    fig, axes = plt.subplots(2, 2, figsize=(14, 14))
    title = (
        f"{image_id} | {info_dict.get('combo', '')} | "
        f"{info_dict.get('condition', '')} | {info_dict.get('replicate', '')} | "
        f"voxel=({voxel_size[0]:.2f},{voxel_size[1]:.3f},{voxel_size[2]:.3f})μm"
    )
    fig.suptitle(title, fontsize=10)

    # DAPI + outlines
    axes[0, 0].imshow(dapi_mip, cmap="gray")
    overlay = np.zeros((*nuc_mip.shape, 4))
    overlay[boundaries] = [0, 1, 1, 1]
    axes[0, 0].imshow(overlay)
    axes[0, 0].set_title(
        f"DAPI + nuclei outlines (n={int(nuc_mip.max())})", fontsize=11
    )
    axes[0, 0].axis("off")

    # Merged RGB
    rgb = np.stack([iav_mip, asyn_mip, dapi_mip], axis=-1)
    rgb = np.clip(rgb / (rgb.max() + 1e-9) * 2.5, 0, 1)
    axes[0, 1].imshow(rgb)
    overlay2 = np.zeros((*nuc_mip.shape, 4))
    overlay2[boundaries] = [1, 1, 1, 1]
    axes[0, 1].imshow(overlay2)
    axes[0, 1].set_title("Merged: R=IAV  G=aSyn  B=DAPI", fontsize=11)
    axes[0, 1].axis("off")

    # aSyn + puncta
    axes[1, 0].imshow(asyn_mip, cmap="Greens")
    overlay3 = np.zeros((*nuc_mip.shape, 4))
    overlay3[boundaries] = [0, 1, 1, 1]
    axes[1, 0].imshow(overlay3)
    n_punc = len(puncta_df)
    if n_punc > 0:
        axes[1, 0].scatter(
            puncta_df["puncta_x"],
            puncta_df["puncta_y"],
            s=15,
            edgecolors="magenta",
            facecolors="none",
            linewidths=0.8,
        )
    axes[1, 0].set_title(f"aSyn MIP + nuclear puncta (n={n_punc})", fontsize=11)
    axes[1, 0].axis("off")

    # IAV
    axes[1, 1].imshow(iav_mip, cmap="Reds")
    overlay4 = np.zeros((*nuc_mip.shape, 4))
    overlay4[boundaries] = [0, 1, 1, 1]
    axes[1, 1].imshow(overlay4)
    info = info_dict.get("infection_info", {})
    n_inf = sum(1 for v in info.values() if v.get("infected", False))
    if n_inf > 0:
        from scipy.ndimage import center_of_mass

        for nid, v in info.items():
            if v.get("infected", False):
                try:
                    cz, cy, cx = center_of_mass(nuclei_labels == nid)
                    axes[1, 1].plot(
                        cx,
                        cy,
                        "*",
                        color="yellow",
                        markersize=15,
                        markeredgecolor="black",
                        markeredgewidth=0.5,
                    )
                except Exception:
                    pass
    axes[1, 1].set_title(
        f"IAV MIP: {n_inf}/{int(nuc_mip.max())} cells IAV+", fontsize=11
    )
    axes[1, 1].axis("off")

    plt.tight_layout()
    plt.savefig(out_path, dpi=QC_DPI, bbox_inches="tight")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════════
# PROCESS ONE IMAGE (main per-FOV pipeline)
# ═══════════════════════════════════════════════════════════════════════════════


def process_one_image(row, out_dirs, use_cellpose, global_iav_thresholds):
    """Process one field of view. Returns ((nc_df, puncta_df, coloc_df, summary), None) or (None, err)."""
    image_id = row["image_id"]
    voxel_size = (
        float(row["voxel_z_um"]),
        float(row["voxel_y_um"]),
        float(row["voxel_x_um"]),
    )
    ch_dapi = int(row["ch_dapi_idx"])
    ch_asyn = int(row["ch_asyn_idx"])
    ch_iav = int(row["ch_iav_idx"])
    combo = row.get("combo", "")
    global_thresh = global_iav_thresholds.get(combo, None)

    if ch_dapi < 0 or ch_asyn < 0 or ch_iav < 0:
        return None, "missing channel indices"

    try:
        channels = load_ome_tiff(
            row["tiff_path"], ch_dapi, ch_asyn, ch_iav, int(row["bit_depth"])
        )
    except Exception as e:
        return None, f"load failed: {e}"

    # Preprocess: DAPI smoothing, aSyn high-pass (5 µm intensity / 2 µm puncta)
    dapi_for_seg = gaussian_smooth_dapi(channels["dapi"], DAPI_SMOOTH_SIGMA_PX)
    asyn_bgsub = rolling_ball_bg_subtract(
        channels["asyn"], ROLLING_BALL_RADIUS_UM, voxel_size
    )
    asyn_tophat = white_tophat_3d(channels["asyn"], TOPHAT_RADIUS_UM, voxel_size)
    iav_bgsub = rolling_ball_bg_subtract(
        channels["iav"], IAV_BGSUB_SIGMA_UM, voxel_size
    )

    # Segment (re-use an existing mask so that it matches the one used for the thresholds)
    label_path = os.path.join(out_dirs["labels"], f"{image_id}_nuclei.tif")
    if os.path.exists(label_path):
        nuclei_labels = tifffile.imread(label_path)
    else:
        try:
            if use_cellpose:
                nuclei_labels = segment_nuclei_cellpose(dapi_for_seg, voxel_size)
            else:
                nuclei_labels = segment_nuclei_watershed_fallback(
                    dapi_for_seg, voxel_size
                )
        except Exception as e:
            print(
                f"  WARNING: Cellpose segmentation failed for {image_id} ({e}); watershed fallback used"
            )
            nuclei_labels = segment_nuclei_watershed_fallback(dapi_for_seg, voxel_size)
        tifffile.imwrite(label_path, nuclei_labels, compression="zlib")

    n_nuclei = int(len(np.unique(nuclei_labels)) - 1)
    if n_nuclei == 0:
        return None, "no nuclei detected"

    cyto_ring = create_cytoplasm_ring(nuclei_labels, CYTORING_WIDTH_UM, voxel_size)
    tifffile.imwrite(
        os.path.join(out_dirs["labels"], f"{image_id}_cytoring.tif"),
        cyto_ring,
        compression="zlib",
    )

    # Classify IAV+ cells on the raw IAV channel with the dual threshold
    infection_info, iav_thresh = classify_infected_cells(
        nuclei_labels, channels["iav"], voxel_size, global_threshold=global_thresh
    )

    # Manual QC correction: two cells in this A/NWS/33 field (Syn1 panel) that
    # passed the automated dual threshold were reassigned to the bystander
    # compartment after visual inspection of the segmentation/QC overlay.
    # This is the only manual intervention in the pipeline (2 of 470 cells).
    if image_id == "4.S33-C__scene05":
        for nid in [31, 36]:
            if nid in infection_info:
                infection_info[nid]["infected"] = False
                print(
                    f"   ℹ️ Manual override: cell {nid} set to infected=False in {image_id}"
                )

    # Quantify
    nc_df = measure_nc_ratio(asyn_bgsub, nuclei_labels, cyto_ring, infection_info)
    puncta_df = detect_nuclear_puncta(asyn_tophat, nuclei_labels)
    nc_df = aggregate_puncta_per_cell(puncta_df, nc_df)
    if len(puncta_df) > 0:
        puncta_df = measure_puncta_radial_position(puncta_df, nuclei_labels, voxel_size)
    chromatin_df, chromatin_stats = correlate_asyn_with_chromatin(
        asyn_bgsub, channels["dapi"], nuclei_labels
    )

    # Colocalisation (high-pass aSyn; raw IAV and DAPI)
    coloc_df = pd.DataFrame()
    if COLOC_ENABLE:
        try:
            chans_for_coloc = {
                "asyn": asyn_bgsub,
                "iav": channels["iav"],
                "dapi": channels["dapi"],
            }
            coloc_df = compute_coloc_for_image(
                chans_for_coloc, nuclei_labels, cyto_ring, do_costes=COLOC_DO_COSTES
            )
        except Exception as e:
            print(f"  Coloc failed: {e}")

    # Add metadata to all DFs
    meta_cols = [
        "combo",
        "antibody",
        "condition",
        "replicate",
        "source_lif",
        "magnification_class",
        "scene_subtype",
    ]
    for df_ in [nc_df, puncta_df, chromatin_df, coloc_df]:
        if len(df_) > 0:
            df_["image_id"] = image_id
            for col in meta_cols:
                df_[col] = row.get(col, "")
            df_["voxel_y_um"] = voxel_size[1]
            df_["voxel_z_um"] = voxel_size[0]

    # Save per-image CSVs
    nc_df.to_csv(
        os.path.join(out_dirs["per_image_csv"], f"{image_id}_cells.csv"), index=False
    )
    puncta_df.to_csv(
        os.path.join(out_dirs["per_image_csv"], f"{image_id}_puncta.csv"), index=False
    )
    chromatin_df.to_csv(
        os.path.join(out_dirs["per_image_csv"], f"{image_id}_chromatin.csv"),
        index=False,
    )
    if len(coloc_df) > 0:
        coloc_df.to_csv(
            os.path.join(out_dirs["per_image_csv"], f"{image_id}_coloc.csv"),
            index=False,
        )

    # QC overlay
    info_dict = {
        "combo": row.get("combo", ""),
        "condition": row.get("condition", ""),
        "replicate": row.get("replicate", ""),
        "infection_info": infection_info,
    }
    qc_path = os.path.join(out_dirs["qc"], f"{image_id}_qc.png")
    try:
        make_qc_overlay(
            channels,
            nuclei_labels,
            cyto_ring,
            puncta_df,
            qc_path,
            image_id,
            info_dict,
            voxel_size,
        )
    except Exception as e:
        print(f"  QC overlay failed: {e}")

    summary = {
        "image_id": image_id,
        "n_nuclei": n_nuclei,
        "n_infected": int(nc_df["infected"].sum())
        if "infected" in nc_df.columns
        else 0,
        "NC_ratio_median": float(nc_df["NC_ratio"].median()),
        "asyn_nuclear_CV_median": float(nc_df["asyn_nuclear_CV"].median())
        if "asyn_nuclear_CV" in nc_df.columns
        else np.nan,
        "n_puncta_total": len(puncta_df),
        "mean_puncta_per_cell": float(nc_df["n_nuclear_puncta"].mean()),
        "iav_threshold": float(iav_thresh),
        "chromatin_pearson_r": chromatin_stats.get("pearson_r", np.nan),
        "coloc_n_cells": len(coloc_df),
    }
    return (nc_df, puncta_df, coloc_df, summary), None


# ═══════════════════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════════════════


def check_already_processed(image_id, output_dir):
    target = os.path.join(output_dir, "per_image_csv", f"{image_id}_cells.csv")
    return os.path.exists(target) and os.path.getsize(target) > 0


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN LOOP + BUILD MASTER CSVS
# ═══════════════════════════════════════════════════════════════════════════════


def main():
    # Load + filter manifest
    df = load_and_filter_manifest()
    if len(df) == 0:
        print("❌ No FOVs to process. Exit.")
        return

    # Compute dual thresholds (requires Control masks; if missing, segment quick)
    global_iav_thresholds = compute_objective_dual_thresholds(
        df, N_CONTROL_FOVS_FOR_THRESHOLD
    )
    print("\nFinal thresholds:")
    for combo, thresh in global_iav_thresholds.items():
        print(
            f"  {combo}: intensity={thresh['intensity']:.4f}, area={thresh['area']} px"
        )

    use_cellpose = not USE_WATERSHED_ONLY

    # Main loop
    print("\n" + "=" * 70)
    print(f"Starting batch at {datetime.now()}, {len(df)} images queued")
    print("=" * 70)

    all_summaries = []
    n_processed = n_skipped = n_failed = 0
    t_start = time.time()

    pbar = tqdm(df.iterrows(), total=len(df), desc="Processing")
    for _, row in pbar:
        image_id = row["image_id"]
        pbar.set_postfix_str(image_id[:40])

        if not OVERWRITE_EXISTING and check_already_processed(image_id, OUTPUT_DIR):
            n_skipped += 1
            continue

        t0 = time.time()
        try:
            result, err = process_one_image(
                row, OUT_DIRS, use_cellpose, global_iav_thresholds
            )
            if result is None:
                print(f"\n[FAIL] {image_id}: {err}")
                n_failed += 1
                continue
            nc_df, puncta_df, coloc_df, summary = result
            summary["elapsed_sec"] = time.time() - t0
            all_summaries.append(summary)
            n_processed += 1
            del nc_df, puncta_df, coloc_df, summary
            gc.collect()
        except Exception as e:
            print(f"\n[ERROR] {image_id}: {e}")
            n_failed += 1
            gc.collect()
            continue

    t_total = time.time() - t_start
    print(f"\n{'=' * 70}")
    print(f"✅ DONE")
    print(f"  Processed : {n_processed}")
    print(f"  Skipped   : {n_skipped} (already done)")
    print(f"  Failed    : {n_failed}")
    print(
        f"  Total time: {t_total / 60:.1f} min ({t_total / max(n_processed, 1):.1f}s/image avg)"
    )

    # Build master CSVs
    print("\n" + "=" * 70)
    print("Building master CSVs...")
    print("=" * 70)
    build_master_csvs()

    print("\n✅ Pipeline complete. Outputs in:")
    print(f"   {OUTPUT_DIR}")


def build_master_csvs():
    csv_dir = os.path.join(OUTPUT_DIR, "per_image_csv")

    # Master cells
    cell_files = sorted(glob.glob(os.path.join(csv_dir, "*_cells.csv")))
    all_cells = []
    for fp in tqdm(cell_files, desc="Loading cells"):
        try:
            all_cells.append(pd.read_csv(fp))
        except Exception as e:
            print(f"Skip {fp}: {e}")
    if all_cells:
        master_cells = pd.concat(all_cells, ignore_index=True)
        master_cells.to_csv(os.path.join(OUTPUT_DIR, "master_cells.csv"), index=False)
        print(f"  master_cells.csv: {len(master_cells)} cells")

    # Master puncta
    punc_files = sorted(glob.glob(os.path.join(csv_dir, "*_puncta.csv")))
    all_punc = [
        pd.read_csv(fp)
        for fp in tqdm(punc_files, desc="Loading puncta")
        if os.path.getsize(fp) > 0
    ]
    all_punc = [d for d in all_punc if len(d) > 0]
    if all_punc:
        master_puncta = pd.concat(all_punc, ignore_index=True)
        master_puncta.to_csv(os.path.join(OUTPUT_DIR, "master_puncta.csv"), index=False)
        print(f"  master_puncta.csv: {len(master_puncta)} puncta")

    # Master coloc
    coloc_files = sorted(glob.glob(os.path.join(csv_dir, "*_coloc.csv")))
    all_coloc = []
    for fp in tqdm(coloc_files, desc="Loading coloc"):
        try:
            d = pd.read_csv(fp)
            if len(d) > 0:
                all_coloc.append(d)
        except Exception:
            pass
    if all_coloc:
        master_coloc = pd.concat(all_coloc, ignore_index=True)
        master_coloc.to_csv(os.path.join(OUTPUT_DIR, "master_coloc.csv"), index=False)
        print(f"  master_coloc.csv: {len(master_coloc)} cells")

        # master_cells_coloc = cells + coloc merged
        if all_cells:
            merge_keys = ["image_id", "nucleus_id"]
            coloc_cols = [
                c
                for c in master_coloc.columns
                if c not in master_cells.columns or c in merge_keys
            ]
            master_coloc_merged = master_cells.merge(
                master_coloc[coloc_cols], on=merge_keys, how="left"
            )
            master_coloc_merged.to_csv(
                os.path.join(OUTPUT_DIR, "master_cells_coloc.csv"), index=False
            )
            n_with = master_coloc_merged["nuc_asyn_iav_pearson_r"].notna().sum()
            print(
                f"  master_cells_coloc.csv: {len(master_coloc_merged)} cells ({n_with} with coloc)"
            )


if __name__ == "__main__":
    main()
