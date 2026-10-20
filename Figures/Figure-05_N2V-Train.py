#!/usr/bin/env python3
import os
os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")

# Raise the soft open-file limit (macOS default 256) for sequential training runs.
import resource
try:
    _soft, _hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    _target = min(_hard, 8192) if _hard != resource.RLIM_INFINITY else 8192
    if _soft < _target:
        resource.setrlimit(resource.RLIMIT_NOFILE, (_target, _hard))
except Exception:
    pass

import gc
import json
import logging
import shutil
import sys
import time
import traceback
import warnings
from datetime import datetime
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch

plt.rcParams["font.family"] = "sans-serif"

# ─── Config ───────────────────────────────────────────────────────────────────

PATCH_DIR  = Path("Confocal/n2v/patches").resolve()
MODELS_DIR = Path("Confocal/n2v/models").resolve()

GROUPS_TO_TRAIN = [
    "Abcam_EPR20535_8bit_asyn",
    "Abcam_EPR20535_8bit_iav",
    "Syn1_16bit_asyn",
    "Syn1_16bit_iav",
]

# Training hyperparameters
PATCH_SIZE   = 64
BATCH_SIZE   = 64
MAX_EPOCHS   = 50
PATIENCE     = 10     
MIN_EPOCHS   = 30      
LR_INITIAL   = 4e-4    
ENABLE_EARLY_STOPPING = False

# Cap on training patches per model (random subsample, reproducible seed).
MAX_PATCHES_PER_GROUP = 20_000
SUBSAMPLE_SEED        = 42

# Disk-space sanity check
MIN_FREE_GB = 5

# Progress checkpoint interval
PROGRESS_EVERY = 10    # Log "elapsed + ETA" every N epochs


# ─── Master logger setup ──────────────────────────────────────────────────────

def setup_master_logger(log_path: Path) -> logging.Logger:
    """Logger that writes to both stdout and a master log file."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("n2v_train")
    logger.setLevel(logging.DEBUG)
    # Clear any existing handlers (re-runs)
    for h in list(logger.handlers):
        logger.removeHandler(h)

    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")

    fh = logging.FileHandler(log_path, mode="a")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging.INFO)
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    return logger

.

def make_csv_metrics_callback(output_path: Path, progress_logger=None):
    """Returns a PL Callback that streams metrics to CSV on each
    validation_epoch_end. Lightweight, no external deps beyond lightning."""
    from pytorch_lightning.callbacks import Callback

    class CSVMetricsCallback(Callback):
        def __init__(self):
            super().__init__()
            self.rows = []
            self.start_time = time.time()
            self.last_progress_epoch = -1

        def _record(self, trainer):
            metrics = {}
            for k, v in trainer.callback_metrics.items():
                try:
                    metrics[k] = float(v.item() if hasattr(v, "item") else v)
                except Exception:
                    pass
            row = {
                "epoch": int(trainer.current_epoch),
                "global_step": int(trainer.global_step),
                "wall_time_s": time.time() - self.start_time,
                **metrics,
            }
            self.rows.append(row)
            # Stream to disk every epoch so a crash doesn't lose history
            pd.DataFrame(self.rows).to_csv(output_path, index=False)

            # Periodic progress log
            if (progress_logger is not None
                and row["epoch"] >= self.last_progress_epoch + PROGRESS_EVERY):
                self.last_progress_epoch = row["epoch"]
                val_loss = metrics.get("val_loss", float("nan"))
                train_loss = metrics.get("train_loss", float("nan"))
                elapsed = row["wall_time_s"]
                # Rough ETA assuming linear remaining epochs (max_epochs ceiling)
                eps_done = row["epoch"] + 1
                eta_s = (MAX_EPOCHS - eps_done) * (elapsed / max(eps_done, 1))
                progress_logger.info(
                    f"    epoch {row['epoch']:3d}  "
                    f"train_loss={train_loss:.5f}  val_loss={val_loss:.5f}  "
                    f"elapsed={elapsed/60:.1f}min  "
                    f"ETA≤{eta_s/60:.1f}min"
                )

        def on_validation_epoch_end(self, trainer, pl_module):
            # Skip the sanity-check validation that PL runs before epoch 0
            if trainer.sanity_checking:
                return
            self._record(trainer)

    return CSVMetricsCallback()


# ─── Plot training curves ─────────────────────────────────────────

def plot_training_curves(history_df: pd.DataFrame, model_dir: Path,
                          group_name: str, best_epoch: int | None,
                          stopped_epoch: int | None):
    """1×2 subplot: full log-scale curves + last-50 zoom."""
    if history_df is None or len(history_df) == 0:
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    has_train = "train_loss" in history_df.columns
    has_val   = "val_loss" in history_df.columns

    # ── Subplot 1: full curves on log Y ────────────────────────────────
    ax = axes[0]
    if has_train:
        ax.plot(history_df["epoch"], history_df["train_loss"],
                label="train_loss", color="#3366aa", lw=1.2)
    if has_val:
        ax.plot(history_df["epoch"], history_df["val_loss"],
                label="val_loss", color="#aa3333", lw=1.2)

    if best_epoch is not None:
        ax.axvline(best_epoch, color="#22aa55", ls="--", lw=1.2,
                   label=f"best (ep {best_epoch})")
    if stopped_epoch is not None and stopped_epoch != best_epoch:
        ax.axvline(stopped_epoch, color="#999999", ls=":", lw=1.2,
                   label=f"stopped (ep {stopped_epoch})")

    ax.set_yscale("log")
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Loss (log scale)")
    ax.set_title(f"{group_name} — full training history")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.3, which="both")

    # ── Subplot 2: zoom on last 50 epochs (linear scale) ───────────────
    ax = axes[1]
    last_n = 50
    sub = history_df.tail(last_n)
    if has_train:
        ax.plot(sub["epoch"], sub["train_loss"],
                label="train_loss", color="#3366aa", lw=1.2)
    if has_val:
        ax.plot(sub["epoch"], sub["val_loss"],
                label="val_loss", color="#aa3333", lw=1.2)

    if best_epoch is not None and best_epoch >= sub["epoch"].min():
        ax.axvline(best_epoch, color="#22aa55", ls="--", lw=1.2,
                   label=f"best (ep {best_epoch})")

    # Final values in legend
    final_train = history_df["train_loss"].iloc[-1] if has_train else float("nan")
    final_val   = history_df["val_loss"].iloc[-1]   if has_val   else float("nan")
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Loss")
    ax.set_title(f"Last {min(last_n, len(history_df))} epochs  ·  "
                 f"final train={final_train:.4f}  val={final_val:.4f}")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    output_path = model_dir / "training_curves.png"
    fig.savefig(output_path, dpi=120, bbox_inches="tight")
    plt.close(fig)


# ─── Locate + canonicalize checkpoints after training ─────────────

def finalize_checkpoints(model_dir: Path, logger: logging.Logger):
    """careamics writes checkpoints under model_dir/checkpoints/ with various
    names. Copy the best + last into stable names at model_dir root.

    Returns: dict with keys 'best_path', 'last_path', 'best_epoch', 'best_val_loss'
    (any may be None if the file wasn't found).
    """
    ckpts = list(model_dir.rglob("*.ckpt"))
    if not ckpts:
        logger.warning(f"  No checkpoints found under {model_dir}")
        return {"best_path": None, "last_path": None,
                "best_epoch": None, "best_val_loss": None}

    # Identify "last" checkpoint (Lightning names it last.ckpt by default)
    last_ckpts = [p for p in ckpts if p.stem == "last"
                  or p.name == "last.ckpt"
                  or "last" in p.stem.lower()]
    last_path = last_ckpts[0] if last_ckpts else None

    # Identify "best" — heuristics:
    #   1. checkpoint with min val_loss embedded in filename
    #   2. otherwise: pick the most recently modified non-last checkpoint
    best_path = None
    best_val = float("inf")
    best_epoch = None
    for p in ckpts:
        if last_path is not None and p == last_path:
            continue
        # try to parse filenames like "epoch=42-val_loss=0.0123.ckpt"
        # or pull val_loss from inside the file if needed
        name = p.stem
        val = None
        ep = None
        for tok in name.replace("=", "-").split("-"):
            tok = tok.strip()
            if not tok:
                continue
            try:
                if tok.startswith("epoch"):
                    ep_str = tok.replace("epoch", "").replace(":", "")
                    ep = int(ep_str) if ep_str else None
                elif tok.replace(".", "", 1).isdigit() and "." in tok:
                    val = float(tok)
            except ValueError:
                pass
        # Try to read actual val_loss from checkpoint if not in name
        if val is None:
            try:
                ckpt = torch.load(p, map_location="cpu", weights_only=False)
                logged_metrics = ckpt.get("hyper_parameters", {}).get("val_loss")
                if logged_metrics is None:
                    # Lightning saves callback metrics in 'callbacks' or directly
                    pass
                if isinstance(logged_metrics, (int, float)):
                    val = float(logged_metrics)
            except Exception:
                pass

        if val is not None and val < best_val:
            best_val = val
            best_path = p
            best_epoch = ep

    # Fallback: if we couldn't identify best by name, pick the one with the
    # smallest file mtime span (oldest = first written = usually the best
    # if ModelCheckpoint(save_top_k=1) was used).
    if best_path is None:
        non_last = [p for p in ckpts if p != last_path]
        if non_last:
            best_path = non_last[0]

    # Copy to canonical paths
    canon_best = model_dir / "best_model.ckpt"
    canon_last = model_dir / "last_model.ckpt"
    if best_path and best_path.exists():
        if best_path.resolve() != canon_best.resolve():
            shutil.copy2(best_path, canon_best)
    if last_path and last_path.exists():
        if last_path.resolve() != canon_last.resolve():
            shutil.copy2(last_path, canon_last)

    return {
        "best_path": str(canon_best) if canon_best.exists() else None,
        "last_path": str(canon_last) if canon_last.exists() else None,
        "best_epoch": best_epoch,
        "best_val_loss": best_val if best_val < float("inf") else None,
        "raw_ckpts": [str(p.relative_to(model_dir)) for p in ckpts],
    }


# ─── Per-group training ───────────────────────────────────────────

def train_one_group(group_name: str, npz_path: Path, model_dir: Path,
                     master_logger: logging.Logger) -> dict:
    """Run end-to-end training for a single (combo, channel) group.

    Returns a dict with keys: status, epochs, best_val_loss, elapsed_sec, error.
    Never raises — caller handles result.
    """
    result = {
        "group": group_name, "status": "FAILED", "epochs": 0,
        "best_val_loss": None, "best_epoch": None,
        "elapsed_sec": 0.0, "elapsed_min": 0.0, "error": None,
    }
    t_start = time.time()

    try:
        # ── 1. Load patches ──────────────────────────────────────────
        master_logger.info(f"  Loading {npz_path.name}...")
        data = np.load(npz_path, allow_pickle=True)
        train_patches = data["train_patches"].astype(np.float32)
        val_patches   = data["val_patches"].astype(np.float32)
        metadata      = data["metadata"].item()
        data.close()

        master_logger.info(f"    train: shape={train_patches.shape} dtype={train_patches.dtype}")
        master_logger.info(f"    val:   shape={val_patches.shape}")
        master_logger.info(f"    data range: [{metadata['data_min']:.2f}, "
                            f"{metadata['data_max']:.2f}]")
        master_logger.info(f"    data stats: mean={metadata['data_mean']:.2f} "
                            f"std={metadata['data_std']:.2f}")

        # ── Subsample if the training set exceeds the cap ──────────
        # The validation set is capped at 1/9 of the training cap (90/10).
        if len(train_patches) > MAX_PATCHES_PER_GROUP:
            master_logger.info(f"  Subsampling train patches: "
                                f"{len(train_patches):,} → {MAX_PATCHES_PER_GROUP:,}")
            rng = np.random.default_rng(SUBSAMPLE_SEED)
            train_idx = rng.choice(len(train_patches),
                                    MAX_PATCHES_PER_GROUP, replace=False)
            train_patches = train_patches[train_idx]

            # Preserve 90/10 → val cap = 1/9 of train cap (~2222 for 20K)
            val_cap = int(MAX_PATCHES_PER_GROUP * (1.0 / 9.0))
            if len(val_patches) > val_cap:
                val_idx = rng.choice(len(val_patches), val_cap, replace=False)
                val_patches = val_patches[val_idx]
            master_logger.info(f"    after subsample: "
                                f"train={train_patches.shape}  val={val_patches.shape}")

        # ── 2. Setup output dir ─────────────────────────────────────
        model_dir.mkdir(parents=True, exist_ok=True)

        # ── 3. Build careamics config ───────────────────────────────
        # careamics configures ModelCheckpoint / EarlyStopping through the
        # training config; only the CSV metrics recorder is passed as a callback.
        from careamics import CAREamist
        from careamics.config import create_n2v_configuration
        from careamics.config.callback_model import EarlyStoppingModel

        config = create_n2v_configuration(
            experiment_name=f"n2v_{group_name}",
            data_type="array",
            axes="SYX",
            patch_size=[PATCH_SIZE, PATCH_SIZE],
            batch_size=BATCH_SIZE,
            num_epochs=MAX_EPOCHS,
            # augmentations=None → default XYFlip + XYRandomRotate90
            logger="none",   # we record metrics via custom callback
            # Save best (top_k=1) + last via ModelCheckpoint defaults
            checkpoint_params={
                "monitor": "val_loss",
                "save_top_k": 1,
                "save_last": True,
                "mode": "min",
            },
            # Single-process data loading (avoids "Too many open files" on
            # macOS when several models are trained in one session).
            train_dataloader_params={
                "num_workers": 0,
                "persistent_workers": False,
                "pin_memory": False,
            },
            val_dataloader_params={
                "num_workers": 0,
                "persistent_workers": False,
                "pin_memory": False,
            },
        )

        if ENABLE_EARLY_STOPPING:
            config.training_config.early_stopping_callback = EarlyStoppingModel(
                monitor="val_loss",
                patience=PATIENCE,    # max 10 in careamics 0.0.15
                mode="min",
                min_delta=0.0,
                verbose=True,
            )
            master_logger.info("  Early stopping enabled (fails in careamics 0.0.15)")
        else:
            config.training_config.early_stopping_callback = None
            master_logger.info(f"  Early stopping disabled — fixed {MAX_EPOCHS} epochs")

        # Log the effective configuration
        master_logger.info("  Config inspection:")
        master_logger.info(f"    config type: {type(config).__name__}")
        if hasattr(config, "training_config"):
            tc = config.training_config
            master_logger.info(f"    training_config: {type(tc).__name__}")
            master_logger.info(f"    num_epochs:      {tc.num_epochs}")
            master_logger.info(f"    early_stopping:  {tc.early_stopping_callback}")
            master_logger.info(f"    checkpoint:      {tc.checkpoint_callback}")

        # ── 4. Callbacks — only our custom CSV recorder ─────────────
        csv_cb = make_csv_metrics_callback(
            output_path=model_dir / "training_log.csv",
            progress_logger=master_logger,
        )

        # ── 5. Create CAREamist with the CSV callback ───────────────
        careamist = CAREamist(
            source=config,
            work_dir=str(model_dir),
            callbacks=[csv_cb],
        )

        # ── 6. Train ────────────────────────────────────────────────
        master_logger.info(f"  Starting training (max {MAX_EPOCHS} epochs, "
                            f"patience {PATIENCE})...")
        t_train = time.time()

        careamist.train(
            train_source=train_patches,
            val_source=val_patches,
        )

        train_elapsed = time.time() - t_train
        master_logger.info(f"  Training loop finished in {train_elapsed/60:.1f}min")

        # ── 7. Locate & canonicalize checkpoints ────────────────────
        ckpt_info = finalize_checkpoints(model_dir, master_logger)

        # ── 8. Read final training history ──────────────────────────
        log_csv = model_dir / "training_log.csv"
        history_df = None
        if log_csv.exists():
            history_df = pd.read_csv(log_csv)
            result["epochs"] = int(history_df["epoch"].max()) + 1
            if "val_loss" in history_df.columns:
                vl = history_df["val_loss"].dropna()
                if len(vl):
                    best_idx = int(vl.idxmin())
                    result["best_val_loss"] = float(vl.iloc[vl.values.argmin()])
                    result["best_epoch"] = int(history_df.loc[best_idx, "epoch"])

        # Override best info with checkpoint-derived value if available
        if ckpt_info["best_val_loss"] is not None:
            if (result["best_val_loss"] is None
                or ckpt_info["best_val_loss"] < result["best_val_loss"]):
                result["best_val_loss"] = ckpt_info["best_val_loss"]
                if ckpt_info["best_epoch"] is not None:
                    result["best_epoch"] = ckpt_info["best_epoch"]

        # ── 9. Plot training curves ─────────────────────────────────
        try:
            plot_training_curves(
                history_df, model_dir, group_name,
                best_epoch=result["best_epoch"],
                stopped_epoch=(result["epochs"] - 1) if result["epochs"] else None,
            )
            master_logger.info(f"  Curves saved: training_curves.png")
        except Exception as e:
            master_logger.warning(f"  Plot failed: {e}")

        # ── 10. Save training_metadata.json ─────────────────────────
        # Record both original and effective (post-subsample) counts
        orig_n_train = int(metadata.get("num_train", len(train_patches)))
        orig_n_val   = int(metadata.get("num_val",   len(val_patches)))
        was_subsampled = orig_n_train > MAX_PATCHES_PER_GROUP

        training_metadata = {
            "group_name": group_name,
            "data": {
                "combo": metadata["combo"],
                "channel": metadata["channel"],
                "n_train": int(len(train_patches)),
                "n_val":   int(len(val_patches)),
                "n_train_original": orig_n_train,
                "n_val_original":   orig_n_val,
                "was_subsampled":   was_subsampled,
                "data_min":  metadata["data_min"],
                "data_max":  metadata["data_max"],
                "data_mean": metadata["data_mean"],
                "data_std":  metadata["data_std"],
                "source_images": metadata.get("source_images", []),
            },
            "config": {
                "patch_size": PATCH_SIZE,
                "batch_size": BATCH_SIZE,
                "max_epochs": MAX_EPOCHS,
                "patience": PATIENCE,
                "min_epochs": MIN_EPOCHS,
                "lr_initial": LR_INITIAL,
                "max_patches_per_group": MAX_PATCHES_PER_GROUP,
                "subsample_seed": SUBSAMPLE_SEED,
                "augmentations": "default (XYFlip + XYRandomRotate90)",
                "framework": "careamics 0.0.15",
            },
            "results": {
                "epochs_run": result["epochs"],
                "best_epoch": result["best_epoch"],
                "best_val_loss": result["best_val_loss"],
                "best_ckpt": ckpt_info["best_path"],
                "last_ckpt": ckpt_info["last_path"],
                "raw_checkpoints": ckpt_info["raw_ckpts"],
                "elapsed_minutes": (time.time() - t_start) / 60.0,
                "timestamp_finished": datetime.now().isoformat(),
            },
        }
        with open(model_dir / "training_metadata.json", "w") as f:
            json.dump(training_metadata, f, indent=2, default=str)

        result["status"] = "OK"

        # ── 11. Cleanup before next group ───────────────────────────
        del careamist, train_patches, val_patches, history_df

    except Exception as e:
        result["error"] = f"{type(e).__name__}: {e}"
        master_logger.error(f"  Training FAILED for {group_name}: {result['error']}")
        master_logger.error(traceback.format_exc())

    finally:
        result["elapsed_sec"] = time.time() - t_start
        result["elapsed_min"] = result["elapsed_sec"] / 60.0
        # Aggressive cleanup
        gc.collect()
        if torch.backends.mps.is_available():
            try:
                torch.mps.empty_cache()
            except Exception:
                pass

    return result


# ─── Main ─────────────────────────────────────────────────────────

def format_summary_table(results: list[dict]) -> str:
    """ASCII summary table for the master log."""
    col_w = [30, 8, 7, 10, 9]
    sep   = "┌" + "┬".join("─" * (w + 2) for w in col_w) + "┐"
    mid   = "├" + "┼".join("─" * (w + 2) for w in col_w) + "┤"
    end   = "└" + "┴".join("─" * (w + 2) for w in col_w) + "┘"
    header_cells = ["Group", "Status", "Epochs", "Best Val", "Time"]
    header = "│" + "│".join(f" {h:<{w}} " for h, w in zip(header_cells, col_w)) + "│"

    lines = [sep, header, mid]
    n_ok = 0
    total_epochs = 0
    total_min = 0.0
    for r in results:
        if r["status"] == "OK":
            status_glyph = "OK"
        elif r["status"] == "SKIPPED":
            status_glyph = "SKIP"
        else:
            status_glyph = "FAIL"
        epochs_s = str(r["epochs"]) if r["epochs"] else "-"
        bv = (f"{r['best_val_loss']:.4f}" if r["best_val_loss"] is not None
              else "-")
        time_s = f"{r['elapsed_min']:.1f} m"
        row_cells = [r["group"], status_glyph, epochs_s, bv, time_s]
        lines.append("│" + "│".join(f" {c:<{w}} "
                                     for c, w in zip(row_cells, col_w)) + "│")
        # SKIPPED = model already trained in an earlier session
        if r["status"] in ("OK", "SKIPPED"):
            n_ok += 1
        total_epochs += r["epochs"] if r["epochs"] else 0
        total_min += r["elapsed_min"]

    lines.append(mid)
    total_cells = [
        "TOTAL", f"{n_ok}/{len(results)}",
        str(total_epochs), "-", f"{total_min/60:.1f} h",
    ]
    lines.append("│" + "│".join(f" {c:<{w}} "
                                 for c, w in zip(total_cells, col_w)) + "│")
    lines.append(end)
    return "\n".join(lines)


def main():
    MODELS_DIR.mkdir(parents=True, exist_ok=True)
    master_log_path = MODELS_DIR / "training_master_log.txt"
    logger = setup_master_logger(master_log_path)

    t_total = time.time()
    timestamp_start = datetime.now().isoformat()

    logger.info("=" * 70)
    logger.info(f"Session start: {timestamp_start}")
    logger.info("=" * 70)

    # ── Sanity checks ────────────────────────────────────────────────
    # 1. MPS available?
    if not torch.backends.mps.is_available():
        logger.error("MPS NOT available — training will use CPU (very slow!)")
        device_str = "cpu"
    else:
        device_str = "mps"
        logger.info("✓ MPS available — using Apple Silicon GPU")
    logger.info(f"  torch={torch.__version__}  device={device_str}")
    logger.info(f"  PYTORCH_ENABLE_MPS_FALLBACK="
                f"{os.environ.get('PYTORCH_ENABLE_MPS_FALLBACK')}")

    # Log the effective open-file limit
    try:
        _soft, _hard = resource.getrlimit(resource.RLIMIT_NOFILE)
        logger.info(f"  fd limit: soft={_soft} hard={_hard}  "
                    f"(dataloader num_workers=0, no worker pipes)")
        if _soft < 1024:
            logger.warning(f"  ⚠ Low fd limit ({_soft}) — if training fails "
                            f"with 'Too many open files', run `ulimit -n 8192` "
                            f"before launching the script")
    except Exception as e:
        logger.warning(f"  Could not read fd limit: {e}")

    # 2. Disk space
    free_gb = shutil.disk_usage(MODELS_DIR).free / 1e9
    if free_gb < MIN_FREE_GB:
        logger.error(f"Free disk space {free_gb:.1f} GB < {MIN_FREE_GB} GB "
                      f"required — aborting")
        return
    logger.info(f"  free disk: {free_gb:.1f} GB OK")

    # 3. careamics version
    try:
        import careamics
        logger.info(f"  careamics={careamics.__version__}")
    except ImportError as e:
        logger.error(f"careamics not importable: {e}")
        return

    # 4. All input .npz files present
    missing = []
    for g in GROUPS_TO_TRAIN:
        p = PATCH_DIR / f"patches_{g}.npz"
        if not p.exists():
            missing.append(p)
    if missing:
        logger.error("Missing patch files — aborting:")
        for p in missing:
            logger.error(f"   {p}")
        return

    # ── Print pre-train summary ──────────────────────────────────────
    logger.info("")
    logger.info(f"Batch {BATCH_SIZE}, {MAX_EPOCHS} epochs (no early stopping), "
                f"max {MAX_PATCHES_PER_GROUP:,} training patches per model")
    logger.info("")

    logger.info(f"Training plan ({len(GROUPS_TO_TRAIN)} groups):")
    plan_total_h = 0.0
    for g in GROUPS_TO_TRAIN:
        p = PATCH_DIR / f"patches_{g}.npz"
        data = np.load(p, allow_pickle=True)
        n_train = len(data["train_patches"])
        n_val = len(data["val_patches"])
        data.close()
        n_effective = min(n_train, MAX_PATCHES_PER_GROUP)
        # Rough runtime estimate: ~0.39 min/epoch per 1000 patches (batch 64, MPS)
        time_per_epoch_min = 0.39 * (n_effective / 1000.0)
        time_full_h = time_per_epoch_min * MAX_EPOCHS / 60.0
        plan_total_h += time_full_h
        sub_note = f" (eff {n_effective:,})" if n_train > n_effective else ""
        logger.info(f"   {g}: train={n_train:,}{sub_note}  val={n_val:,}  "
                    f"→ ~{time_full_h:.1f}h ({MAX_EPOCHS} epochs)")
    logger.info(f"   estimated total ({MAX_EPOCHS}-epoch fixed): ~{plan_total_h:.1f}h")
    logger.info(f"   hyperparams: patch={PATCH_SIZE} batch={BATCH_SIZE} "
                f"max_epochs={MAX_EPOCHS} early_stop={ENABLE_EARLY_STOPPING}")
    logger.info("")

    # ── Per-group training loop ─────────────────────────────────────
    results = []
    for idx, group_name in enumerate(GROUPS_TO_TRAIN, start=1):
        logger.info("")
        logger.info("=" * 70)
        logger.info(f"MODEL {idx}/{len(GROUPS_TO_TRAIN)}: {group_name}")
        logger.info(f"           start: {datetime.now().isoformat()}")
        logger.info(f"           elapsed total: "
                    f"{(time.time()-t_total)/60:.1f} min so far")
        logger.info("=" * 70)

        npz_path  = PATCH_DIR / f"patches_{group_name}.npz"
        model_dir = MODELS_DIR / group_name

        # Skip models that already have a best_model.ckpt (safe re-launch).
        ckpt_path = model_dir / "best_model.ckpt"
        if ckpt_path.exists() and ckpt_path.stat().st_size > 1000:
            logger.info(f"  ⏭️  Skipping {group_name} — best_model.ckpt already exists "
                        f"({ckpt_path.stat().st_size / 1e6:.1f} MB)")
            result = {
                "group": group_name,
                "status": "SKIPPED",
                "epochs": 0,
                "best_val_loss": None,
                "best_epoch": None,
                "elapsed_sec": 0.0,
                "elapsed_min": 0.0,
                "error": None,
            }
            results.append(result)
            continue

        result = train_one_group(group_name, npz_path, model_dir, logger)
        results.append(result)

        if result["status"] == "OK":
            logger.info(f"  ✅ {group_name} done in {result['elapsed_min']:.1f} min "
                        f"({result['epochs']} epochs, "
                        f"best val_loss={result['best_val_loss']})")
        else:
            logger.error(f"  ❌ {group_name} FAILED after {result['elapsed_min']:.1f} min")

        # Write a partial summary after each model
        partial_table = format_summary_table(results)
        partial_path = MODELS_DIR / "summary_partial.txt"
        partial_path.write_text(
            f"Session start: {timestamp_start}\n"
            f"As of:         {datetime.now().isoformat()}\n"
            f"Progress:      {idx}/{len(GROUPS_TO_TRAIN)} groups complete\n\n"
            + partial_table + "\n"
        )

    # ── Final summary ────────────────────────────────────────────────
    elapsed_total_min = (time.time() - t_total) / 60.0
    logger.info("")
    logger.info("=" * 70)
    logger.info("FINAL SUMMARY")
    logger.info("=" * 70)
    summary_text = format_summary_table(results)
    for line in summary_text.split("\n"):
        logger.info(line)
    logger.info("")
    logger.info(f"  total elapsed:   {elapsed_total_min:.1f} min "
                f"({elapsed_total_min/60:.2f} h)")
    logger.info(f"  models dir:      {MODELS_DIR}")
    logger.info(f"  master log:      {master_log_path}")

    # Save final summary text alongside the partial
    (MODELS_DIR / "summary_final.txt").write_text(
        f"Session start:  {timestamp_start}\n"
        f"Session end:    {datetime.now().isoformat()}\n"
        f"Total elapsed:  {elapsed_total_min:.1f} min ({elapsed_total_min/60:.2f} h)\n\n"
        + summary_text + "\n"
    )

    # Mirror results dict to JSON for downstream scripts
    with open(MODELS_DIR / "results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)


if __name__ == "__main__":
    main()
