"""
═══════════════════════════════════════════════════════════════════════════════
01_LIF-to-TIFF-Conversion.py
═══════════════════════════════════════════════════════════════════════════════
"""

import argparse
import os
import sys
import re
import logging
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import tifffile

try:
    from aicsimageio import AICSImage
    from aicsimageio.writers import OmeTiffWriter
except ImportError as e:
    print(f"ERROR: Missing dependency: {e}")
    print("  pip install 'aicsimageio[bioformats]' tifffile pillow")
    print(
        "  Also need Java: brew install openjdk (macOS) or apt-get install default-jre"
    )
    sys.exit(1)


# ═══════════════════════════════════════════════════════════════════════════════
# CONFIG
# ═══════════════════════════════════════════════════════════════════════════════

CHANNEL_PATTERNS = {
    "DAPI": [r"dapi", r"hoechst", r"405", r"blue"],
    "aSyn": [r"(a-?syn)", r"synuclein", r"(488|alexa.?488|af488)", r"green"],
    "IAV": [
        r"iav",
        r"flu",
        r"influenza",
        r"np-?iav",
        r"(568|594|af568|af594)",
        r"yellow",
        r"red",
        r"goat",
    ],
    "CD68": [r"cd-?68", r"(647|af647|alexa.?647)", r"magenta", r"far.?red"],
}

FALLBACK_CHANNEL_ORDER = ["DAPI", "aSyn", "IAV", "CD68"]

SINGLE_PLANE_Z_MAX = 1

# Output OME-TIFF compression
OME_COMPRESSION = "zlib"


# ═══════════════════════════════════════════════════════════════════════════════
# UTILITIES & RE-INDEX READER
# ═══════════════════════════════════════════════════════════════════════════════


def setup_logging(log_path):
    """Setup logging to both console and file."""
    logger = logging.getLogger("preprocess")
    logger.setLevel(logging.INFO)
    # Clear any existing handlers
    logger.handlers = []
    fmt = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S"
    )
    fh = logging.FileHandler(log_path, mode="w")
    fh.setFormatter(fmt)
    ch = logging.StreamHandler()
    ch.setFormatter(fmt)
    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger


def match_channel(name, patterns):
    """
    Match channel name (case-insensitive) against list of regex patterns.
    Returns True if any pattern matches.
    """
    if name is None:
        return False
    name_lower = str(name).lower()
    return any(re.search(p, name_lower) for p in patterns)


def map_channels(channel_names):
    """
    Map channel names → semantic labels ('DAPI', 'aSyn', 'IAV', 'CD68').
    Returns: dict {channel_index: semantic_label}
             channels that cannot be assigned are labelled 'UNKNOWN_<idx>'
    """
    result = {}
    used_labels = set()

    # Pass 1: pattern match
    for idx, name in enumerate(channel_names):
        for label, patterns in CHANNEL_PATTERNS.items():
            if label in used_labels:
                continue
            if match_channel(name, patterns):
                result[idx] = label
                used_labels.add(label)
                break

    # Pass 2: assign unmatched channels in fallback order
    fallback_iter = iter([l for l in FALLBACK_CHANNEL_ORDER if l not in used_labels])
    for idx, name in enumerate(channel_names):
        if idx not in result:
            try:
                label = next(fallback_iter)
                result[idx] = label
                used_labels.add(label)
            except StopIteration:
                result[idx] = f"UNKNOWN_{idx}"

    return result


def detect_bit_depth(img_array):
    """Detect bit depth from numpy dtype + actual max value."""
    dtype = str(img_array.dtype)
    if dtype == "uint8":
        return 8
    if dtype == "uint16":
        # Check if 12-bit data stored as uint16 (max ~4095)
        max_val = int(img_array.max())
        if max_val < 4096:
            return 12  # likely 12-bit, common for Leica
        return 16
    if dtype == "uint32":
        return 32
    if "float" in dtype:
        return 32  # treat as float32
    return -1  # unknown


def read_lif_scene_corrected(lif_path, scene_idx):
    """
    Read one LIF scene and restore the true channel/Z order of stacks acquired
    in frame-sequential mode, where the reader returns the planes in linear
    acquisition order p = c * Z + z with the true indices c = p % C, z = p // C.

    Returns:
        arr: numpy array shape (C, Z, Y, X), dtype uint16/uint8
    """
    img = AICSImage(str(lif_path))
    img.set_scene(scene_idx)

    # Raw planes as returned by the reader (interleaved for frame-sequential scans)
    data = img.get_image_data("CZYX")  # (C, Z, Y, X)

    # Scenes 0 and 1 of 3.HMC3-S33.lif were acquired in standard (non-interleaved)
    # order and are returned unchanged.
    if "3.HMC3-S33" in Path(lif_path).name and scene_idx in (0, 1):
        return data

    C, Z, Y, X = data.shape
    arr = np.zeros((C, Z, Y, X), dtype=data.dtype)

    for c_lin in range(C):
        for z_lin in range(Z):
            p = c_lin * Z + z_lin  # linear position
            c_real = p % C  # true channel
            z_real = p // C  # true Z
            arr[c_real, z_real] = data[c_lin, z_lin]

    return arr


def verify_exported_tiff(tiff_path, expected_arr_czyx, logger):
    """Check that the written OME-TIFF matches the source array (per-plane means)."""
    try:
        img = AICSImage(str(tiff_path))
        actual = img.get_image_data("CZYX")
        C, Z = actual.shape[:2]
        max_diff = 0
        for c in range(C):
            for z in range(Z):
                e = float(expected_arr_czyx[c, z].mean())
                a = float(actual[c, z].mean())
                d = abs(a - e)
                if d > max_diff:
                    max_diff = d
        if max_diff > 1.0:
            logger.warning(
                f"  VERIFY FAILED for {tiff_path.name}: max mean diff = {max_diff:.3f} > 1.0"
            )
            return False
        logger.info(f"  Verify OK: max mean diff = {max_diff:.3f}")
        return True
    except Exception as e:
        logger.error(f"  Verification check raised an exception: {e}")
        return False


# ═══════════════════════════════════════════════════════════════════════════════
# LIF SCANNING
# ═══════════════════════════════════════════════════════════════════════════════


def scan_lif_file(lif_path, logger):
    """
    Read the metadata of one LIF file; returns a list of dicts (one per scene).
    """
    logger.info(f"--- Scanning: {lif_path.name}")
    scenes_info = []

    try:
        img = AICSImage(str(lif_path))
    except Exception as e:
        logger.error(f"Could not open {lif_path.name}: {e}")
        return scenes_info

    n_scenes = len(img.scenes)
    logger.info(f"  {n_scenes} scenes detected")

    for scene_idx, scene_name in enumerate(img.scenes):
        try:
            img.set_scene(scene_idx)
            # AICSImage dimensions: T, C, Z, Y, X
            dim_t = img.dims.T
            dim_c = img.dims.C
            dim_z = img.dims.Z
            dim_y = img.dims.Y
            dim_x = img.dims.X

            # Physical pixel sizes (µm); None if missing from the metadata
            ps = img.physical_pixel_sizes  # (Z, Y, X) in μm or None
            voxel_z = float(ps.Z) if ps.Z is not None else np.nan
            voxel_y = float(ps.Y) if ps.Y is not None else np.nan
            voxel_x = float(ps.X) if ps.X is not None else np.nan

            # Channel names
            ch_names = (
                list(img.channel_names)
                if img.channel_names
                else [f"Ch{i}" for i in range(dim_c)]
            )

            # Bit depth, determined from the first channel
            try:
                small = img.get_image_data("ZYX", T=0, C=0)
                bit_depth = detect_bit_depth(small)
            except Exception:
                bit_depth = -1

            # Map channels → semantic labels
            ch_map = map_channels(ch_names)

            # Find indices for DAPI/aSyn/IAV
            ch_dapi_idx = next((i for i, l in ch_map.items() if l == "DAPI"), -1)
            ch_asyn_idx = next((i for i, l in ch_map.items() if l == "aSyn"), -1)
            ch_iav_idx = next((i for i, l in ch_map.items() if l == "IAV"), -1)
            ch_cd68_idx = next((i for i, l in ch_map.items() if l == "CD68"), -1)

            # Classify scene type
            scene_type = "3D_stack" if dim_z > SINGLE_PLANE_Z_MAX else "single_plane"

            info = {
                "source_lif": lif_path.name,
                "scene_idx": scene_idx,
                "scene_name": scene_name,
                "dim_t": dim_t,
                "dim_c": dim_c,
                "dim_z": dim_z,
                "dim_y": dim_y,
                "dim_x": dim_x,
                "voxel_z_um": voxel_z,
                "voxel_y_um": voxel_y,
                "voxel_x_um": voxel_x,
                "bit_depth": bit_depth,
                "channel_names_raw": "|".join(ch_names),
                "channel_mapping": "|".join(
                    f"{i}:{ch_map[i]}" for i in sorted(ch_map.keys())
                ),
                "ch_dapi_idx": ch_dapi_idx,
                "ch_asyn_idx": ch_asyn_idx,
                "ch_iav_idx": ch_iav_idx,
                "ch_cd68_idx": ch_cd68_idx,
                "scene_type": scene_type,
                "image_id": f"{lif_path.stem}__scene{scene_idx:02d}",
            }
            scenes_info.append(info)

            logger.info(
                f"  scene {scene_idx:2d}: {scene_type:12s} | "
                f"Z={dim_z:3d} YX={dim_y}x{dim_x} | C={dim_c} {ch_names} | "
                f"voxel=({voxel_z:.3f},{voxel_y:.3f},{voxel_x:.3f}) | "
                f"{bit_depth}-bit"
            )

        except Exception as e:
            logger.error(f"  scene {scene_idx} failed: {e}")
            continue

    return scenes_info


# ═══════════════════════════════════════════════════════════════════════════════
# EXPORT
# ═══════════════════════════════════════════════════════════════════════════════


def export_scene_to_ome_tiff(lif_path, scene_idx, out_path, scene_info, logger):
    """
    Export one LIF scene to OME-TIFF with metadata preserved.
    Saved as TCZYX (T = 1; fixed-cell imaging).
    """
    try:
        # Order-corrected array (C, Z, Y, X)
        arr_czyx = read_lif_scene_corrected(lif_path, scene_idx)
        data = arr_czyx[np.newaxis, ...]  # (1, C, Z, Y, X)

        # Channel names and physical pixel sizes from the source metadata
        img = AICSImage(str(lif_path))
        img.set_scene(scene_idx)
        ps = img.physical_pixel_sizes

        # Channel names for the OME metadata
        ch_names = (
            list(img.channel_names)
            if img.channel_names
            else [f"Ch{i}" for i in range(data.shape[1])]
        )

        # Write OME-TIFF
        OmeTiffWriter.save(
            data=data,
            uri=str(out_path),
            dim_order="TCZYX",
            channel_names=ch_names,
            physical_pixel_sizes=ps,
        )

        # Verify
        verify_exported_tiff(out_path, arr_czyx, logger)

        file_size_mb = out_path.stat().st_size / 1024 / 1024
        logger.info(f"  → exported: {out_path.name} ({file_size_mb:.1f} MB)")
        return True

    except Exception as e:
        logger.error(f"  EXPORT FAILED for scene {scene_idx}: {e}")
        return False


def export_single_plane_preview(lif_path, scene_idx, out_path, scene_info, logger):
    """
    Single-plane scenes (Z = 1; typically 4K representative scans): export a PNG
    preview (maximum projection if Z > 1, otherwise the raw plane).
    Channels are combined into an RGB display image (DAPI = blue, aSyn = green, IAV = red).
    """
    try:
        from PIL import Image

        data = read_lif_scene_corrected(lif_path, scene_idx)  # (C, Z, Y, X)

        # MIP if Z>1, else squeeze
        if data.shape[1] > 1:
            data_2d = data.max(axis=1)  # (C, Y, X)
        else:
            data_2d = data[:, 0, :, :]  # (C, Y, X)

        # RGB array (Y, X, 3): R = IAV, G = aSyn, B = DAPI
        rgb = np.zeros((data_2d.shape[1], data_2d.shape[2], 3), dtype=np.uint8)

        def norm_to_uint8(arr):
            """Min-max scale to 0-255 for display only (not used for quantification)."""
            a_min, a_max = float(arr.min()), float(arr.max())
            if a_max - a_min < 1:
                return np.zeros_like(arr, dtype=np.uint8)
            return np.clip((arr - a_min) / (a_max - a_min) * 255, 0, 255).astype(
                np.uint8
            )

        ch_dapi_idx = scene_info["ch_dapi_idx"]
        ch_asyn_idx = scene_info["ch_asyn_idx"]
        ch_iav_idx = scene_info["ch_iav_idx"]

        if 0 <= ch_iav_idx < data_2d.shape[0]:
            rgb[..., 0] = norm_to_uint8(data_2d[ch_iav_idx])
        if 0 <= ch_asyn_idx < data_2d.shape[0]:
            rgb[..., 1] = norm_to_uint8(data_2d[ch_asyn_idx])
        if 0 <= ch_dapi_idx < data_2d.shape[0]:
            rgb[..., 2] = norm_to_uint8(data_2d[ch_dapi_idx])

        # Downsample previews larger than 2048 px
        from PIL import Image as PILImage

        pil_img = PILImage.fromarray(rgb, mode="RGB")
        max_dim = 2048
        if max(pil_img.size) > max_dim:
            scale = max_dim / max(pil_img.size)
            new_size = (int(pil_img.size[0] * scale), int(pil_img.size[1] * scale))
            pil_img = pil_img.resize(new_size, PILImage.LANCZOS)

        pil_img.save(str(out_path), "PNG", optimize=True)
        logger.info(f"  → preview PNG: {out_path.name}")
        return True

    except Exception as e:
        logger.error(f"  PREVIEW FAILED for scene {scene_idx}: {e}")
        return False


# ═══════════════════════════════════════════════════════════════════════════════
# METADATA TEMPLATE GENERATION
# ═══════════════════════════════════════════════════════════════════════════════


def create_metadata_template(manifest_df, out_path, logger):
    template = manifest_df[["image_id", "source_lif"]].copy()
    template["combo"] = ""
    template["antibody"] = ""
    template["condition"] = ""
    template["replicate"] = ""
    template["notes"] = ""
    template["magnification_class"] = ""
    template["scene_subtype"] = ""
    template.to_csv(out_path, index=False)
    logger.info(f"Metadata template written (Scene-level): {out_path}")

    logger.info(
        f"  → FILL the empty columns, then rerun with --metadata-csv "
        f"{out_path.name} --merge-only"
    )


def merge_metadata_into_manifest(manifest_path, metadata_csv_path, logger):
    """
    Merge manifest.csv with the completed metadata CSV → manifest_final.csv.
    """
    manifest = pd.read_csv(manifest_path)
    meta = pd.read_csv(metadata_csv_path)

    # Validate metadata file
    required_cols = [
        "image_id",
        "combo",
        "antibody",
        "condition",
        "replicate",
        "magnification_class",
        "scene_subtype",
        "notes",
    ]
    missing = [c for c in required_cols if c not in meta.columns]
    if missing:
        logger.error(f"Metadata CSV is missing columns: {missing}")
        return None

    # Check empty
    empty_rows = meta[
        meta[
            [
                "combo",
                "antibody",
                "condition",
                "replicate",
                "magnification_class",
                "scene_subtype",
            ]
        ]
        .isnull()
        .any(axis=1)
        | (
            meta[
                [
                    "combo",
                    "antibody",
                    "condition",
                    "replicate",
                    "magnification_class",
                    "scene_subtype",
                ]
            ]
            == ""
        ).any(axis=1)
    ]
    if len(empty_rows) > 0:
        logger.warning(
            f"{len(empty_rows)} metadata rows are incomplete and will be dropped:"
        )
        for _, row in empty_rows.iterrows():
            logger.warning(f"  {row.get('image_id', row.get('source_lif', '?'))}")

    # Validate combo and condition values
    valid_combos = {"Abcam_EPR20535_8bit", "Syn1_16bit"}
    valid_conditions = {"Control", "V1", "V2", "PFF", "PFF+V1", "PFF+V2"}
    bad_combo = meta[~meta["combo"].isin(valid_combos | {""})]
    bad_cond = meta[~meta["condition"].isin(valid_conditions | {""})]
    id_col = "image_id" if "image_id" in meta.columns else "source_lif"
    if len(bad_combo) > 0:
        logger.warning(f"Invalid combo values (expected {valid_combos}):")
        logger.warning(bad_combo[[id_col, "combo"]].to_string())
    if len(bad_cond) > 0:
        logger.warning(f"Invalid condition values (expected {valid_conditions}):")
        logger.warning(bad_cond[[id_col, "condition"]].to_string())

    # Merge
    merged = manifest.merge(
        meta[
            [
                "image_id",
                "combo",
                "antibody",
                "condition",
                "replicate",
                "notes",
                "magnification_class",
                "scene_subtype",
            ]
        ],
        on="image_id",
        how="left",
    )

    out_path = manifest_path.parent / "manifest_final.csv"
    merged.to_csv(out_path, index=False)
    logger.info(f"Final manifest written: {out_path}")
    logger.info(f"  Total scenes: {len(merged)}")
    logger.info(
        f"  3D stacks (for quantification): "
        f"{(merged['scene_type'] == '3D_stack').sum()}"
    )
    logger.info(
        f"  Single-plane (representative only): "
        f"{(merged['scene_type'] == 'single_plane').sum()}"
    )

    if "condition" in merged.columns:
        analyzable = merged[
            (merged["scene_type"] == "3D_stack")
            & merged["condition"].notna()
            & (merged["condition"] != "")
        ]
        if len(analyzable) > 0:
            logger.info("\n=== Breakdown 3D scenes per group ===")
            breakdown = (
                analyzable.groupby(["combo", "condition", "replicate"])
                .size()
                .reset_index(name="n_scenes")
            )
            logger.info("\n" + breakdown.to_string(index=False))

    return out_path


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════


def main():
    parser = argparse.ArgumentParser(
        description="Preprocess LIF files → OME-TIFF + manifest for the HMC3 aSyn confocal analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--lif-dir", type=str, required=True, help="Folder containing the .lif files"
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        required=True,
        help="Output folder (subfolders are created)",
    )
    parser.add_argument(
        "--metadata-csv",
        type=str,
        default=None,
        help="(optional) completed metadata CSV to merge into the manifest",
    )
    parser.add_argument(
        "--merge-only",
        action="store_true",
        help="Only merge metadata; do not re-scan or re-export",
    )
    parser.add_argument(
        "--skip-export",
        action="store_true",
        help="Scan and write the manifest only; do not export TIFFs",
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Overwrite existing TIFFs"
    )
    args = parser.parse_args()

    lif_dir = Path(args.lif_dir).expanduser().resolve()
    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    tiff_dir = out_dir / "tiff_3d"
    png_dir = out_dir / "png_singleplane"
    tiff_dir.mkdir(exist_ok=True)
    png_dir.mkdir(exist_ok=True)

    log_path = out_dir / "preprocess_log.txt"
    logger = setup_logging(log_path)

    logger.info("=" * 70)
    logger.info(f"HMC3 aSyn — LIF preprocessing | {datetime.now()}")
    logger.info("=" * 70)
    logger.info(f"LIF dir:   {lif_dir}")
    logger.info(f"Out dir:   {out_dir}")
    logger.info(
        f"Mode:      {'merge-only' if args.merge_only else 'full scan + export'}"
    )

    manifest_path = out_dir / "manifest.csv"

    # ─── MERGE-ONLY MODE ───
    if args.merge_only:
        if not manifest_path.exists():
            logger.error(
                f"manifest.csv not found at {manifest_path}. "
                "Run the full mode first (without --merge-only)."
            )
            sys.exit(1)
        if args.metadata_csv is None:
            logger.error("--merge-only requires --metadata-csv")
            sys.exit(1)
        metadata_path = Path(args.metadata_csv).expanduser().resolve()
        if not metadata_path.exists():
            logger.error(f"Metadata CSV not found: {metadata_path}")
            sys.exit(1)
        merge_metadata_into_manifest(manifest_path, metadata_path, logger)
        return

    # ─── FULL MODE: scan + export ───
    if not lif_dir.exists():
        logger.error(f"LIF folder not found: {lif_dir}")
        sys.exit(1)

    lif_files = sorted(
        [f for f in lif_dir.glob("*.lif") if not f.name.startswith("._")]
    )
    if len(lif_files) == 0:
        logger.error(f"No .lif files found in {lif_dir}")
        sys.exit(1)
    logger.info(f"Found {len(lif_files)} LIF file(s)")

    # Scan all files
    all_scenes = []
    for lif_path in lif_files:
        scenes = scan_lif_file(lif_path, logger)
        all_scenes.extend(scenes)

    if len(all_scenes) == 0:
        logger.error("No scenes could be read — exiting.")
        sys.exit(1)

    manifest_df = pd.DataFrame(all_scenes)

    # ─── EXPORT ───
    if not args.skip_export:
        logger.info("\n" + "=" * 70)
        logger.info("EXPORTING 3D stacks → OME-TIFF; single-planes → PNG preview")
        logger.info("=" * 70)

        export_status = []
        for _, row in manifest_df.iterrows():
            lif_path = lif_dir / row["source_lif"]
            image_id = row["image_id"]

            if row["scene_type"] == "3D_stack":
                out_tiff = tiff_dir / f"{image_id}.ome.tif"
                if out_tiff.exists() and not args.overwrite:
                    logger.info(f"  SKIP existing: {out_tiff.name}")
                    export_status.append("skipped_existing")
                    continue
                ok = export_scene_to_ome_tiff(
                    lif_path, row["scene_idx"], out_tiff, row.to_dict(), logger
                )
                export_status.append("exported" if ok else "failed")
            else:
                # single-plane → PNG preview
                out_png = png_dir / f"{image_id}.png"
                if out_png.exists() and not args.overwrite:
                    export_status.append("skipped_existing")
                    continue
                ok = export_single_plane_preview(
                    lif_path, row["scene_idx"], out_png, row.to_dict(), logger
                )
                export_status.append("preview" if ok else "failed")

        manifest_df["export_status"] = export_status

        # Add tiff_path column for 3D scenes
        manifest_df["tiff_path"] = manifest_df.apply(
            lambda r: (
                str(tiff_dir / f"{r['image_id']}.ome.tif")
                if r["scene_type"] == "3D_stack"
                else ""
            ),
            axis=1,
        )

    # ─── WRITE MANIFESTS ───
    manifest_df.to_csv(manifest_path, index=False)
    logger.info(f"\nManifest written: {manifest_path}")

    # Template
    template_path = out_dir / "manifest_template_metadata.csv"
    create_metadata_template(manifest_df, template_path, logger)

    # Merge immediately if a completed metadata CSV was supplied
    if args.metadata_csv:
        metadata_path = Path(args.metadata_csv).expanduser().resolve()
        if metadata_path.exists():
            merge_metadata_into_manifest(manifest_path, metadata_path, logger)
        else:
            logger.warning(f"Metadata CSV not found: {metadata_path}")

    # ─── SUMMARY ───
    logger.info("\n" + "=" * 70)
    logger.info("SUMMARY")
    logger.info("=" * 70)
    logger.info(f"Total scenes scanned : {len(manifest_df)}")
    logger.info(
        f"  3D stacks (for quantification): "
        f"{(manifest_df['scene_type'] == '3D_stack').sum()}"
    )
    logger.info(
        f"  Single-plane (representative)  : "
        f"{(manifest_df['scene_type'] == 'single_plane').sum()}"
    )
    logger.info(f"\nBit depth distribution:")
    logger.info(manifest_df["bit_depth"].value_counts().to_string())
    logger.info(f"\nVoxel sizes (Z, Y, X) seen:")
    voxel_summary = (
        manifest_df[["voxel_z_um", "voxel_y_um", "voxel_x_um"]]
        .round(4)
        .drop_duplicates()
    )
    logger.info(voxel_summary.to_string())

    logger.info("\n" + "=" * 70)
    logger.info("NEXT STEPS:")
    logger.info("=" * 70)
    logger.info(f"1. Complete the metadata template: {template_path}")
    logger.info(f"   Required columns: combo, antibody, condition, replicate")
    logger.info(f"   - combo:     'Abcam_EPR20535_8bit' or 'Syn1_16bit'")
    logger.info(f"   - antibody:  'EPR20535' or 'Syn1'")
    logger.info(f"   - condition: 'Control'/'V1'/'V2'/'PFF'/'PFF+V1'/'PFF+V2'")
    logger.info(f"   - replicate: 'R1', 'R2', 'R3', ...")
    logger.info(f"2. Then run:")
    logger.info(f"   python {Path(__file__).name} \\")
    logger.info(f"       --lif-dir {lif_dir} --out-dir {out_dir} \\")
    logger.info(f"       --metadata-csv {template_path} --merge-only")
    logger.info(
        f"3. → manifest_final.csv is used by 02_Cellpose-Segmentation-Quantification.py"
    )


if __name__ == "__main__":
    main()
