#!/usr/bin/env python3

import logging
import time
from pathlib import Path
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import aicsimageio

plt.rcParams["font.family"] = "sans-serif"

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
)
logger = logging.getLogger(__name__)


# ─── Config ───────────────────────────────────────────────────────────────────

TIFF_DIR = Path("Confocal/tiff_3d").resolve()
MANIFEST_FILE = Path("Confocal/manifest_final.csv").resolve()
OUTPUT_DIR = Path("Confocal/n2v").resolve()
PATCH_DIR = OUTPUT_DIR / "patches"

# 8 candidate representative fields + 3 EPR20535 fields used for training only.
# N2V uses every Z-plane; manual_z (4 for the extra fields) is not used here.
MANUAL_Z_MAPPING = {
    "3.HMC3-No-Virus__scene00": 2,
    "3.HMC3-PR8__scene02": 3,
    "3.HMC3-S33__scene01": 6,
    "3.HMC3-S33__scene00": 6,
    "4-C__scene03": 2,
    "4.PR8-C__scene04": 7,
    "4.PR8-C__scene06": 4,
    "4.S33-C__scene05": 9,
    "3.HMC3-No-Virus__scene01": 4,  # N2V training only — placeholder z
    "3.HMC3-PR8__scene04": 4,  # N2V training only — placeholder z
    "3.HMC3-PR8__scene05": 4,  # N2V training only — placeholder z
}

PATCH_SIZE = 64
PATCH_STRIDE = 48  # ~25 % overlap → ~2× patches vs stride=64
TRAIN_VAL_SPLIT = 0.9
RNG_SEED = 42

# Quality filters
ZERO_PIXEL_FRAC_MAX = 0.5  # Reject patches with >50 % zero pixels
PREVIEW_SAMPLE_SIZE = 8

# N2V works on aSyn + IAV only (DAPI excluded — no N2V for nuclei reference)
CHANNELS = [
    {"name": "asyn", "col": "ch_asyn_idx"},
    {"name": "iav", "col": "ch_iav_idx"},
]

# Sanity warning thresholds
MIN_PATCHES_WARN = 1_000
MAX_PATCHES_INFO = 100_000


# ─── Helpers ──────────────────────────────────────────────────────────────────


def load_ometiff(path: Path):
    """Load OME-TIFF as (C, Z, Y, X) array."""
    try:
        img = aicsimageio.AICSImage(path)
        return img.get_image_data("CZYX")
    except Exception as e:
        logger.error(f"Load error {path.name}: {e}")
        return None


def extract_patches_from_slice(slice_2d: np.ndarray, image_id: str, z_idx: int):
    """Extract 64x64 patches from one 2D slice with quality filtering.

    Filter 1: patch p95 < slice p5  → patch is darker than slice's background.
    Filter 2: >50 % pixels exactly 0 → likely outside cell mask.
    Filter 3: NaN/Inf detected      → numerical corruption, log + skip.

    Returns: (patches_list, metadata_list, counts_dict)
    """
    h, w = slice_2d.shape
    slice_bg_thresh = float(np.percentile(slice_2d, 5))

    patches = []
    metadata = []
    counts = {
        "total": 0,
        "rejected_dark": 0,
        "rejected_zero": 0,
        "rejected_nan": 0,
        "kept": 0,
    }

    for y in range(0, h - PATCH_SIZE + 1, PATCH_STRIDE):
        for x in range(0, w - PATCH_SIZE + 1, PATCH_STRIDE):
            counts["total"] += 1
            patch = slice_2d[y : y + PATCH_SIZE, x : x + PATCH_SIZE]

            # Filter 1: skip near-empty patches (darker than slice background)
            if float(np.percentile(patch, 95)) < slice_bg_thresh:
                counts["rejected_dark"] += 1
                continue

            # Filter 2: skip patches with too many zero pixels
            zero_frac = float((patch == 0).sum()) / patch.size
            if zero_frac > ZERO_PIXEL_FRAC_MAX:
                counts["rejected_zero"] += 1
                continue

            # Filter 3: reject NaN/Inf
            if not np.isfinite(patch).all():
                counts["rejected_nan"] += 1
                logger.error(
                    f"NaN/Inf in patch from {image_id} z={z_idx} y={y} x={x} — skipping"
                )
                continue

            patches.append(patch.astype(np.float32, copy=True))
            metadata.append({"image_id": image_id, "z": z_idx, "y": y, "x": x})
            counts["kept"] += 1

    return patches, metadata, counts


def process_one_image(row, channel_col: str, channel_name: str):
    """Load one TIFF, extract patches from every Z-slice."""
    image_id = row["image_id"]
    tiff_path = TIFF_DIR / Path(row["tiff_path"]).name

    if not tiff_path.exists():
        logger.warning(f"    TIFF not found: {image_id} — skip")
        return [], [], None

    ch_idx = int(row[channel_col])
    if ch_idx < 0:
        logger.warning(f"    {image_id}: {channel_name} channel idx={ch_idx} — skip")
        return [], [], None

    print(f"    {image_id}... ", end="", flush=True)
    t_img = time.time()

    data = load_ometiff(tiff_path)
    if data is None:
        return [], [], None

    z_stack = data[ch_idx, :, :, :].astype(np.float32)
    n_z = z_stack.shape[0]
    del data

    all_patches = []
    all_metadata = []
    img_counts = {
        "total": 0,
        "rejected_dark": 0,
        "rejected_zero": 0,
        "rejected_nan": 0,
        "kept": 0,
        "n_z": n_z,
    }

    for z in range(n_z):
        patches, meta, counts = extract_patches_from_slice(z_stack[z], image_id, z)
        all_patches.extend(patches)
        all_metadata.extend(meta)
        for k in ("total", "rejected_dark", "rejected_zero", "rejected_nan", "kept"):
            img_counts[k] += counts[k]

    del z_stack
    print(
        f"Z={n_z}  kept {img_counts['kept']}/{img_counts['total']}  "
        f"(dark={img_counts['rejected_dark']}, "
        f"zero={img_counts['rejected_zero']}, "
        f"nan={img_counts['rejected_nan']})  "
        f"{time.time() - t_img:.1f}s"
    )
    return all_patches, all_metadata, img_counts


def process_group(group_df: pd.DataFrame, channel_col: str, channel_name: str):
    """Extract patches across all images in one (combo, channel) group."""
    all_patches = []
    all_metadata = []
    group_counts = {
        "total": 0,
        "rejected_dark": 0,
        "rejected_zero": 0,
        "rejected_nan": 0,
        "kept": 0,
        "n_images": 0,
    }

    for _, row in group_df.iterrows():
        patches, meta, counts = process_one_image(row, channel_col, channel_name)
        if counts is None:
            continue
        all_patches.extend(patches)
        all_metadata.extend(meta)
        group_counts["n_images"] += 1
        for k in ("total", "rejected_dark", "rejected_zero", "rejected_nan", "kept"):
            group_counts[k] += counts[k]

    return all_patches, all_metadata, group_counts


def split_train_val(patches_arr: np.ndarray, train_frac: float, seed: int):
    """Reproducible 90/10 split via permutation."""
    rng = np.random.default_rng(seed)
    indices = rng.permutation(len(patches_arr))
    split_idx = int(train_frac * len(patches_arr))
    return patches_arr[indices[:split_idx]], patches_arr[indices[split_idx:]]


def render_preview(preview_records, output_path: Path):
    """4 rows × 8 cols grid of random patches (one row per group).

    preview_records: list of (group_name, sample_patches[N,64,64], total_count)
    """
    n_rows = len(preview_records)
    n_cols = PREVIEW_SAMPLE_SIZE

    fig = plt.figure(figsize=(18, 2.3 * n_rows + 1))
    gs = gridspec.GridSpec(
        n_rows,
        n_cols + 1,
        figure=fig,
        width_ratios=[1.3] + [1] * n_cols,
        wspace=0.08,
        hspace=0.25,
        left=0.02,
        right=0.98,
        top=0.92,
        bottom=0.02,
    )

    for row_idx, (group_name, sample_patches, total_count) in enumerate(
        preview_records
    ):
        # Label cell
        lbl_ax = fig.add_subplot(gs[row_idx, 0])
        lbl_ax.text(
            0.5,
            0.5,
            f"{group_name}\n\n{total_count:,} patches",
            fontsize=10,
            fontweight="bold",
            ha="center",
            va="center",
            transform=lbl_ax.transAxes,
        )
        lbl_ax.axis("off")

        for col_idx in range(n_cols):
            ax = fig.add_subplot(gs[row_idx, col_idx + 1])
            if col_idx < len(sample_patches):
                patch = sample_patches[col_idx]
                p1, p99 = np.percentile(patch, [1, 99])
                if p99 - p1 < 1e-8:
                    p99 = p1 + 1.0
                ax.imshow(
                    patch, cmap="gray", vmin=p1, vmax=p99, interpolation="nearest"
                )
            ax.axis("off")

    fig.suptitle(
        f"N2V training patches preview ({PREVIEW_SAMPLE_SIZE} random per group)",
        fontsize=12,
        fontweight="bold",
        y=0.995,
    )
    fig.savefig(output_path, dpi=120, bbox_inches="tight")
    plt.close(fig)


# ─── Main ─────────────────────────────────────────────────────────────────────


def main():
    t_total = time.time()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PATCH_DIR.mkdir(parents=True, exist_ok=True)

    logger.info(f"📁 Output : {OUTPUT_DIR}")
    logger.info(f"📊 TIFF   : {TIFF_DIR}")
    logger.info(
        f"patch={PATCH_SIZE}×{PATCH_SIZE} "
        f"stride={PATCH_STRIDE}  "
        f"train/val={TRAIN_VAL_SPLIT:.0%}/{1 - TRAIN_VAL_SPLIT:.0%}  "
        f"seed={RNG_SEED}"
    )
    logger.info(
        f"quality filters: dark (p95<slice_p5) | "
        f"zero (>{ZERO_PIXEL_FRAC_MAX:.0%} pixels=0) | NaN/Inf reject"
    )
    logger.info("")

    # Load manifest, restrict to the mapped 3D-stack images with all channels
    manifest_df = pd.read_csv(MANIFEST_FILE)
    df_all = manifest_df[
        manifest_df["image_id"].isin(MANUAL_Z_MAPPING.keys())
        & (manifest_df["scene_type"] == "3D_stack")
        & (manifest_df["ch_dapi_idx"] >= 0)
        & (manifest_df["ch_asyn_idx"] >= 0)
        & (manifest_df["ch_iav_idx"] >= 0)
    ].copy()
    df_all["manual_z"] = df_all["image_id"].map(MANUAL_Z_MAPPING)

    logger.info(f"{len(df_all)} images selected (of {len(MANUAL_Z_MAPPING)} mapped):")
    combos = sorted(df_all["combo"].unique())
    for c in combos:
        sub = df_all[df_all["combo"] == c]
        logger.info(f"  combo '{c}': {len(sub)} images")
        for _, r in sub.iterrows():
            logger.info(f"    {r['image_id']}  z={int(r['manual_z'])}")
    logger.info("")

    # ────────────────────────────────────────────────────────────────────
    # Process each (combo, channel) group → save .npz → record preview sample
    # ────────────────────────────────────────────────────────────────────
    preview_records = []
    preview_rng = np.random.default_rng(RNG_SEED + 1)  # independent of split rng

    for combo in combos:
        group_df = df_all[df_all["combo"] == combo]

        for channel in CHANNELS:
            ch_name = channel["name"]
            ch_col = channel["col"]

            logger.info("=" * 60)
            logger.info(f"Group: combo='{combo}'  channel='{ch_name}'")
            logger.info("=" * 60)

            t_group = time.time()
            patches, metadata, gcounts = process_group(group_df, ch_col, ch_name)

            if len(patches) == 0:
                logger.warning(f"  No patches extracted — skip group {combo}_{ch_name}")
                continue

            patches_arr = np.stack(patches).astype(np.float32)
            del patches  # free list-of-arrays

            # Final sanity checks
            if not np.isfinite(patches_arr).all():
                logger.error(f"  NaN/Inf in stacked array — aborting group")
                continue

            n_all_zero = int(np.sum(np.all(patches_arr == 0, axis=(1, 2))))
            if n_all_zero > 0:
                logger.warning(
                    f"  {n_all_zero} all-zero patches slipped through "
                    f"(filter only catches >{ZERO_PIXEL_FRAC_MAX:.0%}=0)"
                )

            n_total = len(patches_arr)
            if n_total < MIN_PATCHES_WARN:
                logger.warning(f"  Only {n_total} patches — risk of overfitting")
            elif n_total > MAX_PATCHES_INFO:
                logger.info(
                    f"  Large dataset: {n_total:,} patches — training will be slow"
                )

            # Train / val split
            train_patches, val_patches = split_train_val(
                patches_arr, TRAIN_VAL_SPLIT, RNG_SEED
            )

            # Data statistics (stored with the patches)
            data_min = float(patches_arr.min())
            data_max = float(patches_arr.max())
            data_mean = float(patches_arr.mean())
            data_std = float(patches_arr.std())

            source_images = sorted({m["image_id"] for m in metadata})
            metadata_dict = {
                "combo": combo,
                "channel": ch_name,
                "source_images": source_images,
                "patch_size": PATCH_SIZE,
                "patch_stride": PATCH_STRIDE,
                "num_patches": n_total,
                "num_train": int(len(train_patches)),
                "num_val": int(len(val_patches)),
                "data_min": data_min,
                "data_max": data_max,
                "data_mean": data_mean,
                "data_std": data_std,
                "counts": gcounts,
                "rng_seed": RNG_SEED,
                "train_val_split": TRAIN_VAL_SPLIT,
                "zero_pixel_frac_max": ZERO_PIXEL_FRAC_MAX,
            }

            # Save .npz
            output_file = PATCH_DIR / f"patches_{combo}_{ch_name}.npz"
            np.savez_compressed(
                output_file,
                train_patches=train_patches,
                val_patches=val_patches,
                metadata=np.array(metadata_dict, dtype=object),
            )

            # Summary log
            logger.info(f"")
            logger.info(f"  ✅ {combo}_{ch_name}")
            logger.info(f"     source images : {gcounts['n_images']}")
            logger.info(f"     attempted     : {gcounts['total']:,} candidate patches")
            logger.info(f"     rejected dark : {gcounts['rejected_dark']:,}")
            logger.info(f"     rejected zero : {gcounts['rejected_zero']:,}")
            logger.info(f"     rejected nan  : {gcounts['rejected_nan']:,}")
            logger.info(f"     kept (total)  : {n_total:,}")
            logger.info(
                f"     train / val   : {len(train_patches):,} / {len(val_patches):,}"
            )
            logger.info(
                f"     data stats    : min={data_min:.2f}  max={data_max:.2f}  "
                f"mean={data_mean:.2f}  std={data_std:.2f}"
            )
            logger.info(
                f"     saved         : {output_file.name}  "
                f"({output_file.stat().st_size / 1e6:.1f} MB)"
            )
            logger.info(f"     elapsed       : {time.time() - t_group:.1f}s")
            logger.info(f"")

            # Sample preview patches (independent rng — does not disturb split)
            sample_n = min(PREVIEW_SAMPLE_SIZE, n_total)
            sample_idx = preview_rng.choice(n_total, sample_n, replace=False)
            preview_records.append(
                (
                    f"{combo}_{ch_name}",
                    patches_arr[sample_idx].copy(),
                    n_total,
                )
            )

            # Free large arrays before next group
            del patches_arr, train_patches, val_patches, metadata

    # ────────────────────────────────────────────────────────────────────
    # Render preview
    # ────────────────────────────────────────────────────────────────────
    logger.info("=" * 60)
    logger.info("Rendering preview...")
    logger.info("=" * 60)
    preview_path = OUTPUT_DIR / "preview_patches.png"
    render_preview(preview_records, preview_path)
    logger.info(f"  Preview saved: {preview_path}")

    elapsed = time.time() - t_total
    logger.info(f"\n{'=' * 60}")
    logger.info(f"✅ ALL DONE in {elapsed:.1f}s ({elapsed / 60:.1f} min)")
    logger.info(f"   Patches dir : {PATCH_DIR}")
    logger.info(f"   Preview     : {preview_path}")
    logger.info(f"{'=' * 60}")


if __name__ == "__main__":
    main()
