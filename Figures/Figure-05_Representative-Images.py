#!/usr/bin/env python3

import os

os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")

import json
import logging
import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.offsetbox import AnchoredOffsetbox, HPacker, TextArea
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from skimage import io as skio
from skimage.measure import find_contours, regionprops
from scipy.ndimage import gaussian_filter, zoom, binary_fill_holes
import aicsimageio
import tifffile

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "DejaVu Sans"]

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
)
logger = logging.getLogger(__name__)


USE_N2V_ASYN = True
HIDE_IAV_IN_NAIVE = False
FIGURE_BG = "white"
TEXT_ON_BG = "black" if FIGURE_BG == "white" else "white"


CONFIG = {
    "base_dir": Path("."),
    "tiff_dir": "Confocal/tiff_3d",
    "manifest_file": "Confocal/manifest_final.csv",
    "cells_csv": "Confocal/results_segmentation/per_image_csv",
    "mask_dir": "Confocal/results_segmentation/labels",
    "n2v_models_dir": "Confocal/n2v/models",
    "output_dir": "Figures/output/Figure5_Representative",
    "dpi_png": 600,
    "scalebar_um": 20,
    "pct_vmin": 40,
    "pct_vmax": 99.5,
    "global_vmax_mult": 1.3,
    "n2v_floor_pct": 10,
}

# Displayed fields and the Z-plane shown for each.
MANUAL_Z_MAPPING = {
    "3.HMC3-No-Virus__scene00": 2,
    "3.HMC3-PR8__scene02": 3,
    "3.HMC3-S33__scene01": 6,
    "4-C__scene03": 2,
    "4.PR8-C__scene06": 4,
    "4.S33-C__scene05": 9,
}

CLONE_DISPLAY = {
    "Abcam_EPR20535_8bit": "EPR20535",
    "Syn1_16bit": "Syn1",
}
CONDITION_DISPLAY = {
    "Control": "Naive",
    "A/PR/8/34": "A/PR/8/34",
    "A/NWS/33": "A/NWS/33",
}
CONDITION_ORDER = ["Control", "A/PR/8/34", "A/NWS/33"]
COMBO_ORDER = ["Syn1_16bit", "Abcam_EPR20535_8bit"]

COMBO_OVERRIDE = {"Abcam_EPR20535_8bit": {"iav_vmax_mult": 1.4, "dapi_vmax_mult": 0.6}}

# Cells shown in panel b (nucleus IDs). ROI label 'i' = bystander,
# 'ii' = infected; a Naive field has a single 'i'.
ROI_OVERRIDE = {
    "4-C__scene03": {"naive": 8},
    "4.PR8-C__scene06": {"bystander": 17, "infected": 19},
    "4.S33-C__scene05": {"bystander": 33, "infected": 7},
    "3.HMC3-No-Virus__scene00": {"naive": 17},
    "3.HMC3-PR8__scene02": {"bystander": 24, "infected": 15},
    "3.HMC3-S33__scene01": {"bystander": 49, "infected": 35},
}

IAV_GAUSSIAN_SIGMA = 0.8
DAPI_GAUSSIAN_SIGMA = 0.5
ASYN_FALLBACK_SIGMA = 0.6
DAPI_PCT_VMIN = 90
DAPI_PCT_VMAX = 99.0

TILE_SIZE = (256, 256)
TILE_OVERLAP = (48, 48)

PANEL_XY_IN = 6.0
PANEL_STRIP_IN = 0.6
PANEL_THUMB_IN = 1.9
THUMB_INTER_GAP = 0.15
GAP_IN = 0.10
PANEL_A_CELL_GAP = 0.55

GROUP_UNDERLINE_LW = 5.0
GROUP_UNDERLINE_COLOR = "#1a1a1a"

CROSSHAIR_LW = 2.5
STRIP_BORDER_LW = 2.0
STRIP_BORDER_COLOR = "#FFD500"
XY_CROSSHAIR_LW = 2.0
XY_CROSSHAIR_COLOR = "#FFD500"

ROI_BOX_COLOR = "white"
ROI_BOX_LW = 2.6
ROI_LABEL_COLOR = "white"

SCALEBAR_COLOR = "white"
SCALEBAR_LW = 3.4

ROI_CROP_SIZE = 280
VALID_CONDITIONS = ["Control", "A/PR/8/34", "A/NWS/33"]

FS_GROUP_TITLE = 30
FS_TARGET_LABEL = 18
FS_SCALEBAR = 14
FS_THUMB_LABEL = 13
FS_ROI_LABEL = 14
FS_ROI_KIND_LABEL = 14

FW_BOLD = "black"


def load_ometiff(path: Path):
    try:
        img = aicsimageio.AICSImage(path)
        return img.get_image_data("CZYX")
    except Exception as e:
        logger.error(f"Load error {path.name}: {e}")
        return None


def load_nuclei_mask(image_id, mask_dir: Path):
    mask_file = mask_dir / f"{image_id}_nuclei.tif"
    if not mask_file.exists():
        return None
    try:
        m = skio.imread(mask_file)
        if m.ndim == 2:
            m = np.expand_dims(m, 0)
        return m.astype(np.uint16)
    except Exception as e:
        logger.warning(f"  ⚠ mask load: {e}")
        return None


def extract_cell_bbox(mask, nucleus_id):
    voxels = np.where(mask == nucleus_id)
    if len(voxels[0]) == 0:
        return None
    return {
        "centroid_y": float(voxels[1].mean()),
        "centroid_x": float(voxels[2].mean()),
    }


def is_control_image(image_id: str) -> bool:
    return "No-Virus" in image_id or image_id.startswith("4-C__")


def find_n2v_model(combo: str, base: Path):
    for c in [
        base / f"{combo}_asyn" / "best_model.ckpt",
        base / f"{combo}_asyn" / "last_model.ckpt",
        base / f"{combo}_asyn" / "checkpoints" / f"n2v_{combo}_asyn.ckpt",
    ]:
        if c.exists():
            return c
    return None


def n2v_predict_slice(careamist, slice_2d):
    for arr, axes, tile in [
        (slice_2d[np.newaxis, ...], "SYX", True),
        (slice_2d, "YX", True),
        (slice_2d[np.newaxis, ...], "SYX", False),
    ]:
        try:
            kw = {
                "source": arr.astype(np.float32),
                "axes": axes,
                "data_type": "array",
                "batch_size": 1,
            }
            if tile:
                kw["tile_size"] = TILE_SIZE
                kw["tile_overlap"] = TILE_OVERLAP
            out = careamist.predict(**kw)
            if isinstance(out, list):
                out = out[0] if out else None
            if out is None:
                continue
            out = np.squeeze(np.asarray(out))
            if out.ndim != 2:
                if out.ndim == 3 and 1 in out.shape:
                    out = out.squeeze()
                else:
                    continue
            if out.shape != slice_2d.shape:
                if (
                    out.shape[0] >= slice_2d.shape[0]
                    and out.shape[1] >= slice_2d.shape[1]
                ):
                    dh = (out.shape[0] - slice_2d.shape[0]) // 2
                    dw = (out.shape[1] - slice_2d.shape[1]) // 2
                    out = out[dh : dh + slice_2d.shape[0], dw : dw + slice_2d.shape[1]]
                else:
                    continue
            return out.astype(np.float32)
        except Exception:
            continue
    return None


def n2v_predict_volume(careamist, vol):
    out = np.zeros_like(vol, dtype=np.float32)
    for z in range(vol.shape[0]):
        s = n2v_predict_slice(careamist, vol[z])
        if s is None:
            return None
        out[z] = s
    return out


def process_iav(channel_3d):
    return gaussian_filter(
        channel_3d.astype(np.float32),
        sigma=(0.0, IAV_GAUSSIAN_SIGMA, IAV_GAUSSIAN_SIGMA),
    ).astype(np.float32)


def process_dapi(channel_3d):
    return gaussian_filter(
        channel_3d.astype(np.float32),
        sigma=(0.3, DAPI_GAUSSIAN_SIGMA, DAPI_GAUSSIAN_SIGMA),
    ).astype(np.float32)


def process_asyn(asyn_raw_3d, careamist):
    """Noise2Void + floor subtraction (USE_N2V_ASYN = True) or Gaussian fallback."""
    if USE_N2V_ASYN and careamist is not None:
        n2v_vol = n2v_predict_volume(careamist, asyn_raw_3d.astype(np.float32))
        if n2v_vol is None:
            return None, None, "failed"
        floor = float(np.percentile(n2v_vol, CONFIG["n2v_floor_pct"]))
        cleaned = np.maximum(n2v_vol - floor, 0.0).astype(np.float32)
        return cleaned, floor, "n2v_floor"
    # Fallback: Gaussian σ
    smoothed = gaussian_filter(
        asyn_raw_3d.astype(np.float32),
        sigma=(0.0, ASYN_FALLBACK_SIGMA, ASYN_FALLBACK_SIGMA),
    ).astype(np.float32)
    return smoothed, 0.0, "gaussian"


def compute_per_combo_ranges(per_image):
    ranges = {}
    for combo in COMBO_ORDER:
        imgs = [d for d in per_image if d["combo"] == combo]
        if not imgs:
            continue
        asyn_pool = np.concatenate([d["asyn_proc_3d"].flatten()[::8] for d in imgs])
        iav_pool = np.concatenate([d["iav_3d"].flatten()[::8] for d in imgs])
        a_vmin = float(np.percentile(asyn_pool, CONFIG["pct_vmin"]))
        a_vmax = float(np.percentile(asyn_pool, CONFIG["pct_vmax"]))
        i_vmin = float(np.percentile(iav_pool, CONFIG["pct_vmin"]))
        i_vmax = float(np.percentile(iav_pool, CONFIG["pct_vmax"]))
        ov = COMBO_OVERRIDE.get(combo, {})
        iav_m = ov.get("iav_vmax_mult", 1.0)
        dapi_m = ov.get("dapi_vmax_mult", 1.0)
        ranges[combo] = {
            "asyn_vmin": a_vmin,
            "asyn_vmax": a_vmax * CONFIG["global_vmax_mult"],
            "iav_vmin": i_vmin,
            "iav_vmax": i_vmax * iav_m * CONFIG["global_vmax_mult"],
            "dapi_mult": dapi_m,
        }
        logger.info(
            f"  '{combo}':  aSyn=({a_vmin:.1f}, {ranges[combo]['asyn_vmax']:.1f})  "
            f"IAV=({i_vmin:.1f}, {ranges[combo]['iav_vmax']:.1f})"
        )
    return ranges


def normalize_channel(ch, vmin, vmax):
    return np.clip((ch - vmin) / (vmax - vmin + 1e-8), 0.0, 1.0)


def create_rgb(asyn_n, iav_n, dapi_n):
    rgb = np.zeros((*asyn_n.shape, 3), dtype=np.float32)
    rgb[:, :, 0] = iav_n
    rgb[:, :, 1] = asyn_n
    rgb[:, :, 2] = dapi_n * 0.45
    return np.clip(rgb, 0.0, 1.0)


def build_rgb_volume(
    asyn_3d, iav_3d, dapi_3d, ranges, dapi_vmin, dapi_vmax_final, is_ctrl
):
    z, y, x = asyn_3d.shape
    rgb_vol = np.zeros((z, y, x, 3), dtype=np.float32)
    for k in range(z):
        a_n = normalize_channel(asyn_3d[k], ranges["asyn_vmin"], ranges["asyn_vmax"])
        d_n = normalize_channel(dapi_3d[k], dapi_vmin, dapi_vmax_final)
        i_n = (
            np.zeros_like(iav_3d[k], dtype=np.float32)
            if is_ctrl
            else normalize_channel(iav_3d[k], ranges["iav_vmin"], ranges["iav_vmax"])
        )
        rgb_vol[k] = create_rgb(a_n, i_n, d_n)
    return rgb_vol


def make_orthogonal_slices(rgb_vol, y_idx, x_idx, vz, vx):
    z_scale = vz / vx
    xz = rgb_vol[:, y_idx, :]
    yz = rgb_vol[:, :, x_idx]
    xz = zoom(xz, (z_scale, 1, 1), order=1)
    yz = zoom(yz, (z_scale, 1, 1), order=1)
    yz = np.transpose(yz, (1, 0, 2))
    return xz, yz


def make_channel_thumb(ch_norm, color_rgb):
    rgb = np.zeros((*ch_norm.shape, 3), dtype=np.float32)
    rgb[:, :, 0] = ch_norm * color_rgb[0]
    rgb[:, :, 1] = ch_norm * color_rgb[1]
    rgb[:, :, 2] = ch_norm * color_rgb[2]
    return np.clip(rgb, 0.0, 1.0)


# ── ROI / cell selection ──────────────────────────────────────────────────────


def find_best_coloc_position(image_id, cells_csv_dir, mask):
    if mask is None:
        return None, None
    cf = cells_csv_dir / f"{image_id}_cells.csv"
    if not cf.exists():
        return None, None
    try:
        cells = pd.read_csv(cf)
        coloc_data = None
        cocf = cells_csv_dir / f"{image_id}_coloc.csv"
        if cocf.exists():
            try:
                coloc_data = pd.read_csv(cocf)
            except:
                pass
        if cells["infected"].sum() > 0:
            inf = cells[cells["infected"] == True].copy()
            if coloc_data is not None:
                inf = inf.merge(
                    coloc_data[["nucleus_id", "nuc_asyn_iav_M1"]],
                    on="nucleus_id",
                    how="left",
                )
                best = inf.nlargest(1, "nuc_asyn_iav_M1").iloc[0]
            else:
                best = inf.nlargest(1, "asyn_nuclear_mean").iloc[0]
        else:
            best = cells.nlargest(1, "asyn_nuclear_mean").iloc[0]
        props = regionprops(mask)
        centroids = {p.label: p.centroid for p in props}
        nid = int(best["nucleus_id"])
        if nid in centroids:
            _, y, x = centroids[nid]
            return int(y), int(x)
    except Exception as e:
        logger.warning(f"  ⚠ coloc: {e}")
    return None, None


def select_roi_cells(image_id, cells_csv_dir, mask, condition):
    """Returns a list of (roi_label, bbox, nucleus_id, kind).

    Convention: roi_label 'i' = Bystander, 'ii' = Infected
    for infected groups. Naive groups get one ROI with label 'i'.

    Honors ROI_OVERRIDE first (manual nucleus_id picks); falls back to
    auto-selection (max NC_ratio for infected, min for bystander, median
    for naive) if no override provided.
    """
    if mask is None:
        return []

    # ── Manual override ──────────────────────────────────────────
    override = ROI_OVERRIDE.get(image_id)
    if override:
        rois = []
        if "bystander" in override:
            bb = extract_cell_bbox(mask, override["bystander"])
            if bb:
                rois.append(("i", bb, int(override["bystander"]), "bystander"))
        if "infected" in override:
            bb = extract_cell_bbox(mask, override["infected"])
            if bb:
                rois.append(("ii", bb, int(override["infected"]), "infected"))
        if "naive" in override:
            bb = extract_cell_bbox(mask, override["naive"])
            if bb:
                rois.append(("i", bb, int(override["naive"]), "naive"))
        if rois:
            return rois

    # ── Auto-selection fallback ──────────────────────────────────
    cf = cells_csv_dir / f"{image_id}_cells.csv"
    if not cf.exists():
        return []
    try:
        cells = pd.read_csv(cf)
        props = regionprops(mask)
        centroids = {p.label: p.centroid for p in props}
        cells["centroid_y"] = cells["nucleus_id"].map(
            lambda nid: centroids[nid][1] if nid in centroids else np.nan
        )
        cells["centroid_x"] = cells["nucleus_id"].map(
            lambda nid: centroids[nid][2] if nid in centroids else np.nan
        )
        cells = cells.dropna(subset=["centroid_y", "centroid_x"])

        rois = []
        if condition in ("A/PR/8/34", "A/NWS/33"):
            byst = cells[cells["infected"] == False]
            if not byst.empty:
                r = byst.nsmallest(1, "NC_ratio").iloc[0]
                nid = int(r["nucleus_id"])
                bb = extract_cell_bbox(mask, nid)
                if bb:
                    rois.append(("i", bb, nid, "bystander"))
            inf = cells[cells["infected"] == True]
            if not inf.empty:
                r = inf.nlargest(1, "NC_ratio").iloc[0]
                nid = int(r["nucleus_id"])
                bb = extract_cell_bbox(mask, nid)
                if bb:
                    rois.append(("ii", bb, nid, "infected"))
        else:
            med = cells["NC_ratio"].median()
            cells["d2m"] = (cells["NC_ratio"] - med).abs()
            picks = cells.nsmallest(1, "d2m")
            for _, r in picks.iterrows():
                nid = int(r["nucleus_id"])
                bb = extract_cell_bbox(mask, nid)
                if bb:
                    rois.append(("i", bb, nid, "naive"))
        return rois
    except Exception as e:
        logger.warning(f"  ⚠ ROI select: {e}")
        return []


# ── Label helpers ───────────────────────────────────────────────────────────


def add_combined_label(ax, clone_name, y_pos=0.04):
    """Bottom-left label 'α-Syn (Clone) / IAV / DAPI' (no background patch)."""
    # Channel-tinted label colours
    asyn_color = "#00FF66"
    iav_color = "#FF3030"
    dapi_color = "#3399FF"
    sep_color = "white"

    parts = [
        TextArea(
            f"α-Syn ({clone_name})",
            textprops=dict(
                color=asyn_color,
                fontweight=FW_BOLD,
                fontsize=FS_TARGET_LABEL,
                family="sans-serif",
                fontname="Helvetica",
            ),
        ),
        TextArea(
            " / ",
            textprops=dict(
                color=sep_color,
                fontweight=FW_BOLD,
                fontsize=FS_TARGET_LABEL,
                family="sans-serif",
                fontname="Helvetica",
            ),
        ),
        TextArea(
            "IAV",
            textprops=dict(
                color=iav_color,
                fontweight=FW_BOLD,
                fontsize=FS_TARGET_LABEL,
                family="sans-serif",
                fontname="Helvetica",
            ),
        ),
        TextArea(
            " / ",
            textprops=dict(
                color=sep_color,
                fontweight=FW_BOLD,
                fontsize=FS_TARGET_LABEL,
                family="sans-serif",
                fontname="Helvetica",
            ),
        ),
        TextArea(
            "DAPI",
            textprops=dict(
                color=dapi_color,
                fontweight=FW_BOLD,
                fontsize=FS_TARGET_LABEL,
                family="sans-serif",
                fontname="Helvetica",
            ),
        ),
    ]
    hbox = HPacker(children=parts, align="center", pad=0, sep=0)
    anchored = AnchoredOffsetbox(
        loc="lower left",
        child=hbox,
        pad=0.2,
        borderpad=0.2,
        bbox_to_anchor=(0.02, y_pos),
        bbox_transform=ax.transAxes,
        frameon=False,
    )
    ax.add_artist(anchored)


def add_scalebar(ax, img_h, img_w, voxel_um, scalebar_um=20, y_pos=0.035):
    """White scale bar (line + text with a black outline; no background patch)."""
    if voxel_um <= 0 or img_w <= 0:
        return
    sb_px = scalebar_um / voxel_um
    w_frac = sb_px / img_w
    actual = scalebar_um
    while w_frac > 0.5 and actual > 5:
        actual = {20: 10, 10: 5}.get(actual, 5)
        sb_px = actual / voxel_um
        w_frac = sb_px / img_w
    x1 = 0.96
    x0 = x1 - w_frac
    ax.plot(
        [x0, x1],
        [y_pos, y_pos],
        color=SCALEBAR_COLOR,
        transform=ax.transAxes,
        linewidth=SCALEBAR_LW,
        solid_capstyle="butt",
    )
    # White text with subtle black outline (path_effects) for legibility on
    # any background — no opaque bbox needed.
    import matplotlib.patheffects as pe

    txt = ax.text(
        (x0 + x1) / 2,
        y_pos + 0.022,
        f"{int(actual)} μm",
        color=SCALEBAR_COLOR,
        fontsize=FS_SCALEBAR,
        fontweight=FW_BOLD,
        ha="center",
        va="bottom",
        family="sans-serif",
        fontname="Helvetica",
        transform=ax.transAxes,
    )
    txt.set_path_effects(
        [
            pe.Stroke(linewidth=2.5, foreground="black"),
            pe.Normal(),
        ]
    )


def draw_roi_box(ax, bbox, crop_size=ROI_CROP_SIZE, label="", color=None, lw=None):
    """Dashed white ROI box with its label."""
    if color is None:
        color = ROI_BOX_COLOR
    if lw is None:
        lw = ROI_BOX_LW
    half = crop_size / 2
    y0 = bbox["centroid_y"] - half
    x0 = bbox["centroid_x"] - half
    rect = Rectangle(
        (x0, y0),
        crop_size,
        crop_size,
        fill=False,
        edgecolor=color,
        linewidth=lw,
        linestyle="--",
    )
    ax.add_patch(rect)
    if label:
        import matplotlib.patheffects as pe

        txt = ax.text(
            x0 + 6,
            y0 - 12,
            label,
            color=ROI_LABEL_COLOR,
            fontsize=FS_ROI_LABEL,
            fontweight=FW_BOLD,
            family="sans-serif",
            fontname="Helvetica",
        )
        txt.set_path_effects(
            [
                pe.Stroke(linewidth=2.5, foreground="black"),
                pe.Normal(),
            ]
        )


def style_strip_axes(ax, with_border=True):
    """Yellow border, no ticks (orthogonal strips)."""
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        if with_border:
            s.set_visible(True)
            s.set_color(STRIP_BORDER_COLOR)
            s.set_linewidth(STRIP_BORDER_LW)
        else:
            s.set_visible(False)


def style_image_axes(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def lock_data_extent(ax, img_2d_or_3d):
    """Lock the axes to the image extent so that an ROI box reaching past the
    image edge cannot rescale the axes."""
    h, w = img_2d_or_3d.shape[:2]
    ax.set_xlim(-0.5, w - 0.5)
    ax.set_ylim(h - 0.5, -0.5)  # y inverted (image convention)


def save_multichannel_tiff(
    out_path, asyn_3d, iav_3d, dapi_3d, vz, vy, vx, label_asyn="alpha-Syn"
):
    stack = np.stack(
        [
            asyn_3d.astype(np.float32),
            iav_3d.astype(np.float32),
            dapi_3d.astype(np.float32),
        ],
        axis=0,
    )
    stack_imagej = stack.transpose(1, 0, 2, 3)[np.newaxis, ...]
    metadata = {
        "axes": "TZCYX",
        "spacing": float(vz),
        "unit": "um",
        "Labels": [label_asyn, "IAV", "DAPI"],
        "channels": 3,
        "mode": "composite",
    }
    res = (1.0 / float(vx), 1.0 / float(vy))
    tifffile.imwrite(
        out_path,
        stack_imagej.astype(np.float32),
        imagej=True,
        resolution=res,
        metadata=metadata,
    )


def render_panel1_for_image(d, combo_ranges, paths, dpi=600, fig_facecolor=None):
    """Layout:

    [thumb α-Syn]  [   XY MAIN   ]  [YZ strip]
    [thumb IAV  ]  [             ]
    [thumb DAPI ]  [             ]
    [           ]  [   XZ strip  ]
    """
    if fig_facecolor is None:
        fig_facecolor = FIGURE_BG

    combo, cond, clone = d["combo"], d["condition"], d["clone"]
    is_ctrl = d["is_control"] and HIDE_IAV_IN_NAIVE
    vx, vz = d["voxel_x_um"], d["voxel_z_um"]
    iid = d["image_id"]
    rng = combo_ranges[combo]

    # DAPI per-image range
    dapi_vmin = float(np.percentile(d["dapi_3d"], DAPI_PCT_VMIN))
    dapi_vmax_base = float(np.percentile(d["dapi_3d"], DAPI_PCT_VMAX))
    dapi_vmax_final = dapi_vmax_base * rng["dapi_mult"]

    # XY slice
    best_z = d["manual_z"]
    a_xy = d["asyn_proc_3d"][best_z]
    i_xy = d["iav_3d"][best_z]
    p_xy = d["dapi_3d"][best_z]
    a_n_xy = normalize_channel(a_xy, rng["asyn_vmin"], rng["asyn_vmax"])
    p_n_xy = normalize_channel(p_xy, dapi_vmin, dapi_vmax_final)
    if is_ctrl:
        i_n_xy = np.zeros_like(i_xy, dtype=np.float32)
    else:
        i_n_xy = normalize_channel(i_xy, rng["iav_vmin"], rng["iav_vmax"])
    rgb_xy = create_rgb(a_n_xy, i_n_xy, p_n_xy)

    # Orthogonal RGB volume
    rgb_vol = build_rgb_volume(
        d["asyn_proc_3d"],
        d["iav_3d"],
        d["dapi_3d"],
        rng,
        dapi_vmin,
        dapi_vmax_final,
        is_ctrl,
    )

    mask = load_nuclei_mask(iid, paths["mask_dir"])
    y_max, x_max = find_best_coloc_position(iid, paths["cells_csv"], mask)
    if y_max is None:
        y_max, x_max = rgb_xy.shape[0] // 2, rgb_xy.shape[1] // 2

    xz, yz = make_orthogonal_slices(rgb_vol, y_max, x_max, vz, vx)

    # Single-channel thumbnails (same hue as in the merge)
    thumb_asyn = make_channel_thumb(a_n_xy, (0.0, 1.0, 0.0))
    if is_ctrl:
        thumb_iav = np.zeros((*i_n_xy.shape, 3), dtype=np.float32)
    else:
        thumb_iav = make_channel_thumb(i_n_xy, (1.0, 0.0, 0.0))
    thumb_dapi = make_channel_thumb(p_n_xy, (0.0, 0.0, 1.0))  # pure blue

    # ── Layout: thumbnails on LEFT, then XY, then YZ on right ─────────
    margin = 0.05
    fig_w = (
        margin
        + PANEL_THUMB_IN
        + GAP_IN
        + PANEL_XY_IN
        + GAP_IN
        + PANEL_STRIP_IN
        + margin
    )
    fig_h = margin + PANEL_STRIP_IN + GAP_IN + PANEL_XY_IN + margin

    fig = plt.figure(figsize=(fig_w, fig_h), dpi=dpi)
    fig.patch.set_facecolor(fig_facecolor)

    def axn(x_in, y_in, w_in, h_in):
        a = fig.add_axes([x_in / fig_w, y_in / fig_h, w_in / fig_w, h_in / fig_h])
        a.set_facecolor(fig_facecolor)
        return a

    y_xy_bot = margin + PANEL_STRIP_IN + GAP_IN
    y_xz_bot = margin

    # Thumbnails on the LEFT side (top: aSyn, mid: IAV, bottom: DAPI)
    # 3 thumbnails + 2 gaps = XY height
    thumb_x = margin
    top_y = y_xy_bot + PANEL_XY_IN - PANEL_THUMB_IN
    mid_y = top_y - PANEL_THUMB_IN - THUMB_INTER_GAP
    bot_y = mid_y - PANEL_THUMB_IN - THUMB_INTER_GAP
    import matplotlib.patheffects as pe

    thumb_specs = [
        (thumb_asyn, top_y, f"α-Syn ({clone})", "#00FF66"),
        (thumb_iav, mid_y, "IAV", "#FF3030"),
        (thumb_dapi, bot_y, "DAPI", "#3399FF"),
    ]
    for thumb_img, y_in, lbl, lc in thumb_specs:
        axt = axn(thumb_x, y_in, PANEL_THUMB_IN, PANEL_THUMB_IN)
        axt.imshow(thumb_img, interpolation="none", aspect="auto")
        lock_data_extent(axt, thumb_img)
        style_image_axes(axt)
        txt = axt.text(
            0.04,
            0.96,
            lbl,
            color=lc,
            fontsize=FS_THUMB_LABEL,
            fontweight=FW_BOLD,
            va="top",
            ha="left",
            transform=axt.transAxes,
            family="sans-serif",
            fontname="Helvetica",
        )
        txt.set_path_effects(
            [
                pe.Stroke(linewidth=2.5, foreground="black"),
                pe.Normal(),
            ]
        )

    # XY main
    xy_x = margin + PANEL_THUMB_IN + GAP_IN
    ax_xy = axn(xy_x, y_xy_bot, PANEL_XY_IN, PANEL_XY_IN)
    ax_xy.imshow(rgb_xy, interpolation="none", aspect="auto")
    lock_data_extent(ax_xy, rgb_xy)
    style_image_axes(ax_xy)
    ax_xy.axhline(
        y_max,
        color=XY_CROSSHAIR_COLOR,
        linestyle=":",
        linewidth=XY_CROSSHAIR_LW,
        alpha=0.75,
    )
    ax_xy.axvline(
        x_max,
        color=XY_CROSSHAIR_COLOR,
        linestyle=":",
        linewidth=XY_CROSSHAIR_LW,
        alpha=0.75,
    )
    add_combined_label(ax_xy, clone)
    add_scalebar(ax_xy, rgb_xy.shape[0], rgb_xy.shape[1], vx, CONFIG["scalebar_um"])

    # ROI boxes
    rois = select_roi_cells(iid, paths["cells_csv"], mask, cond)
    roi_info = []
    for lbl, bb, nid, kind in rois:
        draw_roi_box(ax_xy, bb, label=lbl, lw=1.4)
        roi_info.append({"label": lbl, "bbox": bb, "nucleus_id": nid, "kind": kind})

    # YZ strip (right of XY)
    yz_x = xy_x + PANEL_XY_IN + GAP_IN
    ax_yz = axn(yz_x, y_xy_bot, PANEL_STRIP_IN, PANEL_XY_IN)
    ax_yz.imshow(yz, interpolation="none", aspect="auto")
    lock_data_extent(ax_yz, yz)
    ax_yz.axhline(
        y_max,
        color=XY_CROSSHAIR_COLOR,
        linestyle=":",
        linewidth=CROSSHAIR_LW,
        alpha=0.75,
    )
    style_strip_axes(ax_yz, with_border=True)

    # XZ strip (below XY)
    ax_xz = axn(xy_x, y_xz_bot, PANEL_XY_IN, PANEL_STRIP_IN)
    ax_xz.imshow(xz, interpolation="none", aspect="auto")
    lock_data_extent(ax_xz, xz)
    ax_xz.axvline(
        x_max,
        color=XY_CROSSHAIR_COLOR,
        linestyle=":",
        linewidth=CROSSHAIR_LW,
        alpha=0.75,
    )
    style_strip_axes(ax_xz, with_border=True)

    return fig, roi_info, mask, dapi_vmin, dapi_vmax_final, rgb_xy, fig_w, fig_h


def compute_zoom_crop(asyn_raw_3d, bbox, crop_size=ROI_CROP_SIZE):
    asyn_mip = asyn_raw_3d.max(axis=0).astype(np.float32)
    half = crop_size // 2
    cy = int(round(bbox["centroid_y"]))
    cx = int(round(bbox["centroid_x"]))
    y0, y1 = cy - half, cy + half
    x0, x1 = cx - half, cx + half
    if y0 < 0:
        y1 -= y0
        y0 = 0
    if y1 > asyn_mip.shape[0]:
        y0 -= y1 - asyn_mip.shape[0]
        y1 = asyn_mip.shape[0]
    if x0 < 0:
        x1 -= x0
        x0 = 0
    if x1 > asyn_mip.shape[1]:
        x0 -= x1 - asyn_mip.shape[1]
        x1 = asyn_mip.shape[1]
    y0 = max(0, y0)
    y1 = min(asyn_mip.shape[0], y1)
    x0 = max(0, x0)
    x1 = min(asyn_mip.shape[1], x1)
    return asyn_mip[y0:y1, x0:x1], (y0, y1, x0, x1)


def render_panel2_zoom(
    asyn_raw_3d,
    mask_3d,
    nucleus_id,
    bbox,
    vx,
    vmin,
    vmax,
    dpi=600,
    show_colorbar=True,
    clone_name="Syn1",
    kind_label="",
):
    """Single-cell zoom of the RAW α-synuclein MIP with the nucleus contour."""
    asyn_crop, (y0, y1, x0, x1) = compute_zoom_crop(asyn_raw_3d, bbox)

    fig_w = 5.0
    fig_h = 5.0
    fig = plt.figure(figsize=(fig_w, fig_h), dpi=dpi)
    fig.patch.set_facecolor(FIGURE_BG)

    ax = fig.add_axes([0.05, 0.05, 0.82 if show_colorbar else 0.9, 0.9])
    ax.set_facecolor(FIGURE_BG)
    im = ax.imshow(
        asyn_crop,
        cmap="inferno",
        vmin=vmin,
        vmax=vmax,
        interpolation="none",
        aspect="equal",
    )
    lock_data_extent(ax, asyn_crop)

    # Nucleus contour
    if mask_3d is not None and nucleus_id is not None:
        try:
            nuc_2d = (mask_3d == nucleus_id).max(axis=0)
            nuc_crop = nuc_2d[y0:y1, x0:x1]
            filled = binary_fill_holes(nuc_crop)
            cs = find_contours(filled.astype(float), 0.5)
            big = sorted([c for c in cs if len(c) >= 20], key=len, reverse=True)
            for c in big[:1]:
                ax.plot(c[:, 1], c[:, 0], color="white", linewidth=1.4, alpha=0.95)
        except Exception:
            pass

    style_image_axes(ax)
    add_scalebar(ax, asyn_crop.shape[0], asyn_crop.shape[1], vx, CONFIG["scalebar_um"])

    # Bottom-left label: "α-Syn (Clone)" + cell type
    import matplotlib.patheffects as pe

    parts = [f"α-Syn ({clone_name})"]
    if kind_label:
        parts.append(kind_label)
    label_text = "\n".join(parts)
    txt = ax.text(
        0.04,
        0.04,
        label_text,
        color="white",
        fontsize=FS_TARGET_LABEL,
        fontweight=FW_BOLD,
        va="bottom",
        ha="left",
        family="sans-serif",
        fontname="Helvetica",
        transform=ax.transAxes,
    )
    txt.set_path_effects(
        [
            pe.Stroke(linewidth=2.5, foreground="black"),
            pe.Normal(),
        ]
    )

    if show_colorbar:
        cax = fig.add_axes([0.90, 0.18, 0.025, 0.64])
        cb = fig.colorbar(im, cax=cax)
        cb.ax.tick_params(colors=TEXT_ON_BG, labelsize=10)
        cb.outline.set_edgecolor(TEXT_ON_BG)
        cb.set_label(
            "α-Syn fluorescence intensity (a.u.)",
            color=TEXT_ON_BG,
            fontsize=12,
            fontweight=FW_BOLD,
            family="sans-serif",
            fontname="Helvetica",
            labelpad=6,
        )

    return fig


def _draw_panel1_cell(
    fig, d, combo_ranges, paths, fig_w_total, fig_h_total, cell_x_in, cell_y_in
):
    """Render one image's orthogonal layout into the combined figure at offset."""
    combo, cond, clone = d["combo"], d["condition"], d["clone"]
    is_ctrl = d["is_control"] and HIDE_IAV_IN_NAIVE
    vx, vz = d["voxel_x_um"], d["voxel_z_um"]
    iid = d["image_id"]
    rng = combo_ranges[combo]

    dapi_vmin = float(np.percentile(d["dapi_3d"], DAPI_PCT_VMIN))
    dapi_vmax_base = float(np.percentile(d["dapi_3d"], DAPI_PCT_VMAX))
    dapi_vmax_final = dapi_vmax_base * rng["dapi_mult"]

    best_z = d["manual_z"]
    a_xy = d["asyn_proc_3d"][best_z]
    i_xy = d["iav_3d"][best_z]
    p_xy = d["dapi_3d"][best_z]
    a_n_xy = normalize_channel(a_xy, rng["asyn_vmin"], rng["asyn_vmax"])
    p_n_xy = normalize_channel(p_xy, dapi_vmin, dapi_vmax_final)
    if is_ctrl:
        i_n_xy = np.zeros_like(i_xy, dtype=np.float32)
    else:
        i_n_xy = normalize_channel(i_xy, rng["iav_vmin"], rng["iav_vmax"])
    rgb_xy = create_rgb(a_n_xy, i_n_xy, p_n_xy)

    rgb_vol = build_rgb_volume(
        d["asyn_proc_3d"],
        d["iav_3d"],
        d["dapi_3d"],
        rng,
        dapi_vmin,
        dapi_vmax_final,
        is_ctrl,
    )

    mask = load_nuclei_mask(iid, paths["mask_dir"])
    y_max, x_max = find_best_coloc_position(iid, paths["cells_csv"], mask)
    if y_max is None:
        y_max, x_max = rgb_xy.shape[0] // 2, rgb_xy.shape[1] // 2
    xz, yz = make_orthogonal_slices(rgb_vol, y_max, x_max, vz, vx)

    thumb_asyn = make_channel_thumb(a_n_xy, (0.0, 1.0, 0.0))
    if is_ctrl:
        thumb_iav = np.zeros((*i_n_xy.shape, 3), dtype=np.float32)
    else:
        thumb_iav = make_channel_thumb(i_n_xy, (1.0, 0.0, 0.0))
    thumb_dapi = make_channel_thumb(p_n_xy, (0.0, 0.0, 1.0))

    margin = 0.05

    def axn(x_in_local, y_in_local, w_in, h_in):
        x_abs = cell_x_in + x_in_local
        y_abs = cell_y_in + y_in_local
        a = fig.add_axes(
            [
                x_abs / fig_w_total,
                y_abs / fig_h_total,
                w_in / fig_w_total,
                h_in / fig_h_total,
            ]
        )
        a.set_facecolor(FIGURE_BG)
        return a

    y_xy_bot = margin + PANEL_STRIP_IN + GAP_IN
    y_xz_bot = margin

    # Thumbs on LEFT
    thumb_x = margin
    top_y = y_xy_bot + PANEL_XY_IN - PANEL_THUMB_IN
    mid_y = top_y - PANEL_THUMB_IN - THUMB_INTER_GAP
    bot_y = mid_y - PANEL_THUMB_IN - THUMB_INTER_GAP
    import matplotlib.patheffects as pe

    for thumb_img, y_in, lbl, lc in [
        (thumb_asyn, top_y, f"α-Syn ({clone})", "#00FF66"),
        (thumb_iav, mid_y, "IAV", "#FF3030"),
        (thumb_dapi, bot_y, "DAPI", "#3399FF"),
    ]:
        axt = axn(thumb_x, y_in, PANEL_THUMB_IN, PANEL_THUMB_IN)
        axt.imshow(thumb_img, interpolation="none", aspect="auto")
        lock_data_extent(axt, thumb_img)
        style_image_axes(axt)
        txt = axt.text(
            0.04,
            0.96,
            lbl,
            color=lc,
            fontsize=FS_THUMB_LABEL,
            fontweight=FW_BOLD,
            va="top",
            ha="left",
            transform=axt.transAxes,
            family="sans-serif",
            fontname="Helvetica",
        )
        txt.set_path_effects(
            [
                pe.Stroke(linewidth=2.5, foreground="black"),
                pe.Normal(),
            ]
        )

    # XY main
    xy_x = margin + PANEL_THUMB_IN + GAP_IN
    ax_xy = axn(xy_x, y_xy_bot, PANEL_XY_IN, PANEL_XY_IN)
    ax_xy.imshow(rgb_xy, interpolation="none", aspect="auto")
    lock_data_extent(ax_xy, rgb_xy)
    style_image_axes(ax_xy)
    ax_xy.axhline(
        y_max,
        color=XY_CROSSHAIR_COLOR,
        linestyle=":",
        linewidth=XY_CROSSHAIR_LW,
        alpha=0.75,
    )
    ax_xy.axvline(
        x_max,
        color=XY_CROSSHAIR_COLOR,
        linestyle=":",
        linewidth=XY_CROSSHAIR_LW,
        alpha=0.75,
    )
    add_combined_label(ax_xy, clone)
    add_scalebar(ax_xy, rgb_xy.shape[0], rgb_xy.shape[1], vx, CONFIG["scalebar_um"])
    rois = select_roi_cells(iid, paths["cells_csv"], mask, cond)
    for lbl, bb, nid, kind in rois:
        draw_roi_box(ax_xy, bb, label=lbl, lw=1.2)

    # YZ strip
    yz_x = xy_x + PANEL_XY_IN + GAP_IN
    ax_yz = axn(yz_x, y_xy_bot, PANEL_STRIP_IN, PANEL_XY_IN)
    ax_yz.imshow(yz, interpolation="none", aspect="auto")
    lock_data_extent(ax_yz, yz)
    ax_yz.axhline(
        y_max,
        color=XY_CROSSHAIR_COLOR,
        linestyle=":",
        linewidth=CROSSHAIR_LW,
        alpha=0.75,
    )
    style_strip_axes(ax_yz, with_border=True)

    # XZ strip
    ax_xz = axn(xy_x, y_xz_bot, PANEL_XY_IN, PANEL_STRIP_IN)
    ax_xz.imshow(xz, interpolation="none", aspect="auto")
    lock_data_extent(ax_xz, xz)
    ax_xz.axvline(
        x_max,
        color=XY_CROSSHAIR_COLOR,
        linestyle=":",
        linewidth=CROSSHAIR_LW,
        alpha=0.75,
    )
    style_strip_axes(ax_xz, with_border=True)


def build_combined_panel_A(per_image, combo_ranges, paths, output_path, dpi=600):
    n_rows, n_cols = 2, 3
    margin = 0.05
    sub_fig_w = (
        margin
        + PANEL_THUMB_IN
        + GAP_IN
        + PANEL_XY_IN
        + GAP_IN
        + PANEL_STRIP_IN
        + margin
    )
    sub_fig_h = margin + PANEL_STRIP_IN + GAP_IN + PANEL_XY_IN + margin

    header_h = 1.10
    cell_gap = PANEL_A_CELL_GAP

    fig_w_total = n_cols * sub_fig_w + (n_cols - 1) * cell_gap
    fig_h_total = n_rows * sub_fig_h + (n_rows - 1) * cell_gap + header_h

    fig = plt.figure(figsize=(fig_w_total, fig_h_total), dpi=dpi)
    fig.patch.set_facecolor(FIGURE_BG)

    by_key = {(d["combo"], d["condition"]): d for d in per_image}

    # ── Column headers + underline bars ───────────────────────
    # Title text vertical center, bar lower
    title_y_frac = 1.0 - (header_h * 0.40) / fig_h_total
    bar_y_frac = 1.0 - (header_h * 0.78) / fig_h_total

    from matplotlib.lines import Line2D

    for c, cond in enumerate(CONDITION_ORDER):
        col_x_in = c * (sub_fig_w + cell_gap)
        x_center_in = col_x_in + sub_fig_w / 2
        # Bar spans from leftmost content (margin) to rightmost (sub_fig_w − margin)
        bar_x0_in = col_x_in + margin
        bar_x1_in = col_x_in + sub_fig_w - margin

        # Column header
        fig.text(
            x_center_in / fig_w_total,
            title_y_frac,
            CONDITION_DISPLAY[cond],
            ha="center",
            va="center",
            fontsize=FS_GROUP_TITLE,
            fontweight=FW_BOLD,
            color=TEXT_ON_BG,
            family="sans-serif",
            fontname="Helvetica",
        )

        # Underline bar
        bar = Line2D(
            [bar_x0_in / fig_w_total, bar_x1_in / fig_w_total],
            [bar_y_frac, bar_y_frac],
            color=GROUP_UNDERLINE_COLOR,
            linewidth=GROUP_UNDERLINE_LW,
            transform=fig.transFigure,
            solid_capstyle="butt",
            clip_on=False,
        )
        fig.add_artist(bar)

    for r, combo in enumerate(COMBO_ORDER):
        for c, cond in enumerate(CONDITION_ORDER):
            d = by_key.get((combo, cond))
            if d is None:
                continue
            cell_x_in = c * (sub_fig_w + cell_gap)
            cell_y_in = (n_rows - 1 - r) * (sub_fig_h + cell_gap)
            _draw_panel1_cell(
                fig,
                d,
                combo_ranges,
                paths,
                fig_w_total,
                fig_h_total,
                cell_x_in,
                cell_y_in,
            )

    fig.savefig(output_path, dpi=dpi, facecolor=FIGURE_BG, edgecolor="none")
    plt.close(fig)


def build_combined_panel_B(zoom_records, output_path, dpi=600):
    """Panel b: wider gaps between groups than between the two cells of a group.

    Column layout (left → right):
      0: Naive
      1: (i) Bystander  PR8       — intra-gap from infected at col 2
      2: (ii) Infected  PR8
      3: (i) Bystander  S33       — intra-gap from infected at col 4
      4: (ii) Infected  S33

    Gaps: inter-group (Naive↔PR8, PR8↔S33) is ~3× wider than intra-group.
    """
    n_rows, n_cols = 2, 5
    panel_in = 4.5
    intra_gap = 0.30  # within combo (Bystander ↔ Infected)
    inter_gap = 0.95  # between combos (Naive ↔ PR8, PR8 ↔ S33)
    header_h = 1.20
    cbar_section_w = 2.2

    # Column descriptors (cond, roi_label, kind_label, group_id)
    col_descriptors = [
        ("Control", "i", "Naive", "g_naive"),
        ("A/PR/8/34", "i", "(i) Bystander", "g_pr8"),
        ("A/PR/8/34", "ii", "(ii) Infected", "g_pr8"),
        ("A/NWS/33", "i", "(i) Bystander", "g_s33"),
        ("A/NWS/33", "ii", "(ii) Infected", "g_s33"),
    ]

    # Compute x position for each column based on group transitions
    col_x_in = []
    x = 0.0
    prev_group = None
    for _, _, _, group in col_descriptors:
        if prev_group is None:
            pass
        elif group == prev_group:
            x += intra_gap
        else:
            x += inter_gap
        col_x_in.append(x)
        x += panel_in
        prev_group = group
    grid_right_in = col_x_in[-1] + panel_in

    fig_w_total = grid_right_in + cbar_section_w
    fig_h_total = n_rows * panel_in + (n_rows - 1) * intra_gap + header_h

    fig = plt.figure(figsize=(fig_w_total, fig_h_total), dpi=dpi)
    fig.patch.set_facecolor(FIGURE_BG)

    # ── Group titles + underline bars ──────────────────────────
    title_y_frac = 1.0 - (header_h * 0.40) / fig_h_total
    bar_y_frac = 1.0 - (header_h * 0.78) / fig_h_total

    # Compute group spans from col_descriptors
    groups_layout = []
    cur_group = None
    cur_start_col = 0
    cur_cond = None
    for c, (cond, _, _, group) in enumerate(col_descriptors):
        if group != cur_group:
            if cur_group is not None:
                groups_layout.append((cur_cond, cur_start_col, c - cur_start_col))
            cur_group = group
            cur_start_col = c
            cur_cond = cond
    groups_layout.append(
        (cur_cond, cur_start_col, len(col_descriptors) - cur_start_col)
    )

    from matplotlib.lines import Line2D

    for cond, start_col, n_span in groups_layout:
        # Bar spans from start of first col to end of last col in the group
        x_start_in = col_x_in[start_col]
        x_end_in = col_x_in[start_col + n_span - 1] + panel_in
        x_center_in = (x_start_in + x_end_in) / 2

        fig.text(
            x_center_in / fig_w_total,
            title_y_frac,
            CONDITION_DISPLAY[cond],
            ha="center",
            va="center",
            fontsize=FS_GROUP_TITLE,
            fontweight=FW_BOLD,
            color=TEXT_ON_BG,
            family="sans-serif",
            fontname="Helvetica",
        )
        bar = Line2D(
            [x_start_in / fig_w_total, x_end_in / fig_w_total],
            [bar_y_frac, bar_y_frac],
            color=GROUP_UNDERLINE_COLOR,
            linewidth=GROUP_UNDERLINE_LW,
            transform=fig.transFigure,
            solid_capstyle="butt",
            clip_on=False,
        )
        fig.add_artist(bar)

    by_key = {(r["combo"], r["condition"], r["roi_label"]): r for r in zoom_records}

    import matplotlib.patheffects as pe

    last_im = None
    for r_idx, combo in enumerate(COMBO_ORDER):
        for c_idx, (cond, roi, kind_tag, _) in enumerate(col_descriptors):
            rec = by_key.get((combo, cond, roi))
            if rec is None:
                continue
            cell_x_in = col_x_in[c_idx]
            cell_y_in = (n_rows - 1 - r_idx) * (panel_in + intra_gap)
            ax = fig.add_axes(
                [
                    cell_x_in / fig_w_total,
                    cell_y_in / fig_h_total,
                    panel_in / fig_w_total,
                    panel_in / fig_h_total,
                ]
            )
            ax.set_facecolor(FIGURE_BG)
            im = ax.imshow(
                rec["asyn_crop"],
                cmap="inferno",
                vmin=rec["vmin"],
                vmax=rec["vmax"],
                interpolation="none",
                aspect="equal",
            )
            lock_data_extent(ax, rec["asyn_crop"])
            last_im = im

            if rec.get("nucleus_contour"):
                for contour in rec["nucleus_contour"]:
                    ax.plot(
                        contour[:, 1],
                        contour[:, 0],
                        color="white",
                        linewidth=1.4,
                        alpha=0.95,
                    )

            style_image_axes(ax)
            add_scalebar(
                ax,
                rec["asyn_crop"].shape[0],
                rec["asyn_crop"].shape[1],
                rec["voxel_x_um"],
                CONFIG["scalebar_um"],
            )

            label_text = f"α-Syn ({rec['clone']})\n{kind_tag}"
            txt = ax.text(
                0.04,
                0.04,
                label_text,
                color="white",
                fontsize=FS_TARGET_LABEL,
                fontweight=FW_BOLD,
                va="bottom",
                ha="left",
                family="sans-serif",
                fontname="Helvetica",
                transform=ax.transAxes,
            )
            txt.set_path_effects(
                [
                    pe.Stroke(linewidth=2.5, foreground="black"),
                    pe.Normal(),
                ]
            )

    # ── Colorbar fully inside figure ───────────────────────────
    cbar_x_left_in = grid_right_in + 0.45
    cbar_width_in = 0.18
    cax = fig.add_axes(
        [cbar_x_left_in / fig_w_total, 0.18, cbar_width_in / fig_w_total, 0.62]
    )
    sm = ScalarMappable(norm=Normalize(vmin=0, vmax=1), cmap="inferno")
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax, ticks=[0, 0.5, 1])
    cb.ax.set_yticklabels(
        ["Low", "Mid", "High"],
        color=TEXT_ON_BG,
        fontsize=12,
        fontweight=FW_BOLD,
        family="sans-serif",
        fontname="Helvetica",
    )
    cb.outline.set_edgecolor(TEXT_ON_BG)
    cb.set_label(
        "α-Syn fluorescence intensity (a.u.)",
        color=TEXT_ON_BG,
        fontsize=14,
        fontweight=FW_BOLD,
        family="sans-serif",
        fontname="Helvetica",
        labelpad=14,
    )

    fig.savefig(output_path, dpi=dpi, facecolor=FIGURE_BG, edgecolor="none")
    plt.close(fig)


def compute_panel2_records_all(per_image, paths):
    """Zoom records: one per Naive field, two (bystander, infected) per infected field."""
    records = []
    for d in per_image:
        iid = d["image_id"]
        combo = d["combo"]
        cond = d["condition"]
        clone = d["clone"]
        vx = d["voxel_x_um"]
        mask = load_nuclei_mask(iid, paths["mask_dir"])
        rois = select_roi_cells(iid, paths["cells_csv"], mask, cond)

        for lbl, bbox, nid, kind in rois:
            asyn_crop, (y0, y1, x0, x1) = compute_zoom_crop(d["asyn_raw_3d"], bbox)
            vmin = float(np.percentile(asyn_crop, 30))
            vmax = float(np.percentile(asyn_crop, 99.5))
            nuc_contour = None
            if mask is not None:
                try:
                    nuc_2d = (mask == nid).max(axis=0)
                    nuc_crop = nuc_2d[y0:y1, x0:x1]
                    filled = binary_fill_holes(nuc_crop)
                    cs = find_contours(filled.astype(float), 0.5)
                    big = sorted([c for c in cs if len(c) >= 20], key=len, reverse=True)
                    nuc_contour = big[:1]
                except Exception:
                    pass
            records.append(
                {
                    "image_id": iid,
                    "combo": combo,
                    "condition": cond,
                    "clone": clone,
                    "roi_label": lbl,
                    "kind": kind,
                    "asyn_crop": asyn_crop,
                    "vmin": vmin,
                    "vmax": vmax,
                    "voxel_x_um": vx,
                    "nucleus_contour": nuc_contour,
                    "nucleus_id": nid,
                    "bbox": bbox,
                    "mask": mask,
                }
            )
    return records


# ─── Loader ───────────────────────────────────────────────────────────────────


def load_and_process_all(paths, manifest_df):
    from careamics import CAREamist

    per_combo_model = {}
    if USE_N2V_ASYN:
        for combo in COMBO_ORDER:
            p = find_n2v_model(combo, paths["n2v_models"])
            if p is None:
                logger.error(f"  Missing N2V model: {combo}")
                return None
            per_combo_model[combo] = p
            logger.info(f"  ✓ {combo}_asyn model: {p.name}")
    else:
        logger.info(f"  USE_N2V_ASYN=False → Gaussian σ={ASYN_FALLBACK_SIGMA} fallback")

    per_image = []
    for combo in COMBO_ORDER:
        rows = manifest_df[manifest_df["combo"] == combo]
        if rows.empty:
            continue
        careamist = None
        if USE_N2V_ASYN:
            logger.info(f"\nLoad N2V model for {combo}...")
            careamist = CAREamist(source=str(per_combo_model[combo]))

        rows = rows.assign(
            _cr=rows["condition"].map({c: i for i, c in enumerate(CONDITION_ORDER)})
        ).sort_values("_cr")

        for _, row in rows.iterrows():
            iid = row["image_id"]
            print(f"  [{iid}] ", end="", flush=True)
            t_img = time.time()
            tiff_path = paths["tiff_dir"] / Path(row["tiff_path"]).name
            if not tiff_path.exists():
                print("TIFF MISSING")
                continue
            data = load_ometiff(tiff_path)
            if data is None:
                print("LOAD FAIL")
                continue
            ai = int(row["ch_asyn_idx"])
            ii = int(row["ch_iav_idx"])
            di = int(row["ch_dapi_idx"])
            asyn_raw_3d = data[ai].astype(np.float32)
            iav_raw_3d = data[ii].astype(np.float32)
            dapi_raw_3d = data[di].astype(np.float32)
            del data
            print(f"load {time.time() - t_img:.1f}s ", end="", flush=True)

            iav_3d = process_iav(iav_raw_3d)
            dapi_3d = process_dapi(dapi_raw_3d)

            t_a = time.time()
            asyn_proc_3d, floor, mode = process_asyn(asyn_raw_3d, careamist)
            if asyn_proc_3d is None:
                print("\n    aSyn FAIL — skip")
                continue
            print(
                f"asyn-{mode} {time.time() - t_a:.1f}s floor={floor:.1f}  "
                f"total {time.time() - t_img:.1f}s"
            )

            per_image.append(
                {
                    "image_id": iid,
                    "combo": combo,
                    "condition": row["condition"],
                    "antibody": row["antibody"],
                    "clone": CLONE_DISPLAY.get(combo, row["antibody"]),
                    "manual_z": int(row["manual_z"]),
                    "voxel_z_um": float(row["voxel_z_um"]),
                    "voxel_y_um": float(row["voxel_y_um"]),
                    "voxel_x_um": float(row["voxel_x_um"]),
                    "is_control": is_control_image(iid),
                    "asyn_raw_3d": asyn_raw_3d,
                    "asyn_proc_3d": asyn_proc_3d,
                    "iav_3d": iav_3d,
                    "dapi_3d": dapi_3d,
                    "asyn_floor": floor,
                    "asyn_mode": mode,
                }
            )

        if careamist is not None:
            del careamist
        import gc

        gc.collect()

    return per_image


def main():
    base = CONFIG["base_dir"]
    paths = {
        "tiff_dir": (base / CONFIG["tiff_dir"]).resolve(),
        "manifest": (base / CONFIG["manifest_file"]).resolve(),
        "cells_csv": (base / CONFIG["cells_csv"]).resolve(),
        "mask_dir": (base / CONFIG["mask_dir"]).resolve(),
        "n2v_models": (base / CONFIG["n2v_models_dir"]).resolve(),
        "output_dir": (base / CONFIG["output_dir"]).resolve(),
    }
    output_dir = paths["output_dir"]
    panel1_dir = output_dir / "Panel1_individual"
    panel2_dir = output_dir / "Panel2_individual"
    tiff_dir = output_dir / "tiff-last"
    panel1_dir.mkdir(parents=True, exist_ok=True)
    panel2_dir.mkdir(parents=True, exist_ok=True)
    tiff_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 70)
    logger.info(
        f"Representative figures (Fig. 5a–b)  ·  USE_N2V_ASYN={USE_N2V_ASYN}  "
        f"·  bg={FIGURE_BG}"
    )
    logger.info("=" * 70)

    for name, p in paths.items():
        if name == "output_dir":
            continue
        if not p.exists():
            logger.error(f"  ✗ Missing: {name} = {p}")
            return

    manifest_df = pd.read_csv(paths["manifest"])
    df = manifest_df[
        manifest_df["image_id"].isin(MANUAL_Z_MAPPING.keys())
        & (manifest_df["scene_type"] == "3D_stack")
        & (manifest_df["ch_dapi_idx"] >= 0)
        & (manifest_df["ch_asyn_idx"] >= 0)
        & (manifest_df["ch_iav_idx"] >= 0)
    ].copy()
    df["manual_z"] = df["image_id"].map(MANUAL_Z_MAPPING)
    logger.info(f"  Filtered to {len(df)} images")

    logger.info("\nPHASE 1 — Load + process")
    t0 = time.time()
    per_image = load_and_process_all(paths, df)
    if not per_image:
        logger.error("No images processed")
        return
    logger.info(f"  Phase 1 done in {time.time() - t0:.1f}s")

    logger.info("\nPHASE 2 — Per-combo display ranges")
    combo_ranges = compute_per_combo_ranges(per_image)

    logger.info("\nPHASE 3 — Individual Panel 1 + TIFF + Panel 2 zooms")
    panel2_records = []
    for d in per_image:
        iid = d["image_id"]
        safe_cond = d["condition"].replace("/", "_")
        clone = d["clone"]
        scene_tag = iid.split("__")[-1]

        # Panel 1
        t = time.time()
        fig1, roi_info, mask, dvmin, dvmax, rgb_xy, fw, fh = render_panel1_for_image(
            d, combo_ranges, paths, dpi=CONFIG["dpi_png"]
        )
        p1_name = f"{clone}__{safe_cond}__{scene_tag}__panel1.png"
        # Exact figure size (no bbox_inches="tight")
        fig1.savefig(
            panel1_dir / p1_name,
            dpi=CONFIG["dpi_png"],
            facecolor=FIGURE_BG,
            edgecolor="none",
        )
        plt.close(fig1)
        logger.info(f"  ✓ Panel1 {p1_name}  ({time.time() - t:.1f}s)")

        # TIFF
        tiff_name = f"{clone}__{safe_cond}__{scene_tag}.tif"
        save_multichannel_tiff(
            tiff_dir / tiff_name,
            d["asyn_proc_3d"],
            d["iav_3d"],
            d["dapi_3d"],
            d["voxel_z_um"],
            d["voxel_y_um"],
            d["voxel_x_um"],
            label_asyn=f"alpha-Syn ({clone})",
        )
        logger.info(f"  ✓ TIFF   {tiff_name}")

    # Compute all Panel-2 zoom records (1 per naive, 2 per infected)
    logger.info("\nPHASE 3b — Panel 2 zoom records")
    panel2_records = compute_panel2_records_all(per_image, paths)
    logger.info(f"  {len(panel2_records)} zoom records ready")

    # Save individual Panel 2 PNGs
    for rec in panel2_records:
        iid = rec["image_id"]
        safe_cond = rec["condition"].replace("/", "_")
        scene_tag = iid.split("__")[-1]
        clone = rec["clone"]
        roi = rec["roi_label"]
        kind = rec["kind"]
        # Find raw aSyn from per_image
        d = next(x for x in per_image if x["image_id"] == iid)
        # Display label: Naive has no (i)/(ii) prefix
        if kind == "naive":
            kind_display = "Naive"
        elif kind == "bystander":
            kind_display = f"({roi}) Bystander"
        elif kind == "infected":
            kind_display = f"({roi}) Infected"
        else:
            kind_display = kind.capitalize()
        fig2 = render_panel2_zoom(
            d["asyn_raw_3d"],
            rec["mask"],
            rec["nucleus_id"],
            rec["bbox"],
            d["voxel_x_um"],
            rec["vmin"],
            rec["vmax"],
            dpi=CONFIG["dpi_png"],
            show_colorbar=True,
            clone_name=clone,
            kind_label=kind_display,
        )
        p2_name = f"{clone}__{safe_cond}__{scene_tag}__zoom_{roi}_{kind}.png"
        fig2.savefig(
            panel2_dir / p2_name,
            dpi=CONFIG["dpi_png"],
            facecolor=FIGURE_BG,
            edgecolor="none",
        )
        plt.close(fig2)
        logger.info(f"  ✓ Panel2 {p2_name}")

    logger.info("\nPHASE 4 — Combined Panel A")
    panelA = output_dir / "PanelA_combined_orthogonal.png"
    build_combined_panel_A(
        per_image, combo_ranges, paths, panelA, dpi=CONFIG["dpi_png"]
    )
    logger.info(f"  ✓ {panelA.name}")

    logger.info("\nPHASE 5 — Combined Panel B (2×5 grid)")
    panelB = output_dir / "PanelB_combined_zoom.png"
    build_combined_panel_B(panel2_records, panelB, dpi=CONFIG["dpi_png"])
    logger.info(f"  ✓ {panelB.name}")

    meta = {
        "timestamp": pd.Timestamp.now().isoformat(),
        "USE_N2V_ASYN": USE_N2V_ASYN,
        "FIGURE_BG": FIGURE_BG,
        "n_images": len(per_image),
        "n_zoom_records": len(panel2_records),
        "combo_ranges": {
            c: {k: float(v) for k, v in vals.items() if isinstance(v, (int, float))}
            for c, vals in combo_ranges.items()
        },
    }
    with open(output_dir / "figure_metadata.json", "w") as f:
        json.dump(meta, f, indent=2)

    logger.info("\n" + "=" * 70)
    logger.info(f"Done in {time.time() - t0:.1f}s ({(time.time() - t0) / 60:.1f} min)")
    logger.info(f"   Output: {output_dir}")
    logger.info("=" * 70)


if __name__ == "__main__":
    main()
