from __future__ import annotations

from collections import defaultdict
import json
import os
from pathlib import Path
import sys

import h5py
import numpy as np
import pandas as pd


SLEEPFM_COMMIT = "2bcbae04c3592f61352addb7ac3d4193f0a3ca25"
PROBABILITY_COLUMNS = ["Wake", "N1", "N2", "N3", "REM"]


def _default_sleepfm_root() -> Path:
    env = os.environ.get("SLEEPFM_ROOT")
    if env:
        return Path(env)
    raise EnvironmentError(
        "SLEEPFM_ROOT is not set. Point it at your sleepfm-clinical checkout."
    )


def _load_official_modules(sleepfm_root: Path):
    package_root = sleepfm_root / "sleepfm"
    if not (package_root / "models" / "models.py").exists():
        raise FileNotFoundError(
            f"SleepFM-Clinical was not found at {sleepfm_root}. Run ./scripts/setup.sh first."
        )
    package_text = str(package_root)
    if package_text not in sys.path:
        sys.path.insert(0, package_text)
    from models.dataset import SetTransformerDataset, collate_fn
    from models.models import SetTransformer, SleepEventLSTMClassifier
    from preprocessing.preprocessing import EDFToHDF5Converter

    return EDFToHDF5Converter, SetTransformerDataset, collate_fn, SetTransformer, SleepEventLSTMClassifier


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text())


def _state_dict_for_model(checkpoint, model):
    state = checkpoint.get("state_dict", checkpoint) if isinstance(checkpoint, dict) else checkpoint
    model_keys = model.state_dict().keys()
    expects_module = next(iter(model_keys)).startswith("module.")
    has_module = next(iter(state)).startswith("module.")
    if has_module and not expects_module:
        state = {key.removeprefix("module."): value for key, value in state.items()}
    elif expects_module and not has_module:
        state = {f"module.{key}": value for key, value in state.items()}
    return state


def _resolve_device(torch, requested: str):
    if requested == "auto":
        requested = "cuda" if torch.cuda.is_available() else "cpu"
    if requested.startswith("cuda") and not torch.cuda.is_available():
        raise RuntimeError("CUDA was requested but is not available")
    return torch.device(requested)


def preprocess_edf(
    edf_path: str | Path,
    hdf5_path: str | Path,
    sleepfm_root: str | Path | None = None,
    force: bool = False,
) -> Path:
    """Convert EDF with the official SleepFM converter at its required 128 Hz."""
    sleepfm_root = Path(sleepfm_root or _default_sleepfm_root())
    converter_class, *_ = _load_official_modules(sleepfm_root)
    edf_path, hdf5_path = Path(edf_path), Path(hdf5_path)
    hdf5_path.parent.mkdir(parents=True, exist_ok=True)
    if hdf5_path.exists() and not force:
        return hdf5_path
    converter = converter_class(
        root_dir=str(edf_path.parent),
        target_dir=str(hdf5_path.parent),
        resample_rate=128,
        num_threads=1,
    )
    converter.convert(str(edf_path), str(hdf5_path))
    return hdf5_path


def _validate_modalities(hdf5_path: Path, channel_groups: dict, modality_types: list[str]) -> dict[str, list[str]]:
    with h5py.File(hdf5_path, "r") as handle:
        available = list(handle.keys())
    selected = {
        modality: [channel for channel in available if channel in channel_groups[modality]]
        for modality in modality_types
    }
    missing = [modality for modality, channels in selected.items() if not channels]
    if missing:
        raise ValueError(
            f"EDF lacks channels recognized by SleepFM for modalities {missing}. "
            f"Available channels: {available}"
        )
    return selected


def infer_hypnodensity(
    edf_path: str | Path,
    output_dir: str | Path = "outputs",
    sleepfm_root: str | Path | None = None,
    device: str = "auto",
    batch_size: int = 1,
    force: bool = False,
    prepared_hdf5_path: str | Path | None = None,
    cache_tag: str | None = None,
) -> pd.DataFrame:
    """Return official SleepFM softmax probabilities at its native 5-s resolution.

    The official pretrained encoder and sleep-staging classifier are used without
    altering their weights. A cached CSV is returned on later calls.
    """
    import torch
    from torch.utils.data import DataLoader

    sleepfm_root = Path(sleepfm_root or _default_sleepfm_root())
    edf_path, output_dir = Path(edf_path), Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = edf_path.stem if cache_tag is None else f"{edf_path.stem}_{cache_tag}"
    probability_path = output_dir / f"{stem}_sleepfm_hypnodensity_5s.csv"
    hdf5_path = (
        Path(prepared_hdf5_path)
        if prepared_hdf5_path is not None
        else output_dir / f"{stem}_sleepfm_128hz.hdf5"
    )
    if probability_path.exists() and not force:
        return pd.read_csv(probability_path)

    (
        _,
        dataset_class,
        collate_fn,
        base_model_class,
        staging_model_class,
    ) = _load_official_modules(sleepfm_root)
    if prepared_hdf5_path is None:
        hdf5_path = preprocess_edf(edf_path, hdf5_path, sleepfm_root, force=force)
    elif not hdf5_path.exists():
        raise FileNotFoundError(f"Prepared SleepFM HDF5 not found: {hdf5_path}")
    checkpoints = sleepfm_root / "sleepfm" / "checkpoints"
    base_path = checkpoints / "model_base"
    staging_path = checkpoints / "model_sleep_staging"
    base_config = _load_json(base_path / "config.json")
    staging_config = _load_json(staging_path / "config.json")
    channel_groups = _load_json(sleepfm_root / "sleepfm" / "configs" / "channel_groups.json")
    base_config["num_workers"] = 1
    selected_channels = _validate_modalities(
        hdf5_path, channel_groups, base_config["modality_types"]
    )

    runtime_device = _resolve_device(torch, device)
    base_model = base_model_class(
        base_config["in_channels"],
        base_config["patch_size"],
        base_config["embed_dim"],
        base_config["num_heads"],
        base_config["num_layers"],
        pooling_head=base_config["pooling_head"],
        dropout=0.0,
    ).to(runtime_device)
    checkpoint = torch.load(base_path / "best.pt", map_location=runtime_device)
    base_model.load_state_dict(_state_dict_for_model(checkpoint, base_model))
    base_model.eval()

    dataset = dataset_class(base_config, channel_groups, hdf5_paths=[str(hdf5_path)], split="test")
    if len(dataset) == 0:
        raise ValueError("EDF is shorter than SleepFM's 5-minute encoder window")
    loader = DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=0,
        collate_fn=collate_fn,
    )
    embedding_chunks: dict[str, list[tuple[int, np.ndarray]]] = defaultdict(list)
    with torch.inference_mode():
        for batch_data, mask_list, _, _, chunk_starts in loader:
            for modality, values, mask in zip(base_config["modality_types"], batch_data, mask_list):
                values = values.to(runtime_device, dtype=torch.float32)
                mask = mask.to(runtime_device, dtype=torch.bool)
                embeddings = base_model(values, mask)[1].cpu().numpy()
                for start, embedding in zip(chunk_starts, embeddings):
                    embedding_chunks[modality].append((int(start), embedding))

    modality_embeddings = []
    for modality in base_config["modality_types"]:
        ordered = [array for _, array in sorted(embedding_chunks[modality], key=lambda item: item[0])]
        modality_embeddings.append(np.concatenate(ordered, axis=0))
    common_steps = min(map(len, modality_embeddings))
    x_data = np.stack([values[:common_steps] for values in modality_embeddings], axis=0)
    x_data = torch.from_numpy(x_data).unsqueeze(0).to(runtime_device, dtype=torch.float32)
    mask = torch.zeros((1, len(modality_embeddings), common_steps), dtype=torch.bool, device=runtime_device)

    params = staging_config["model_params"]
    staging_model = staging_model_class(**params).to(runtime_device)
    checkpoint = torch.load(staging_path / "best.pth", map_location=runtime_device)
    staging_model.load_state_dict(_state_dict_for_model(checkpoint, staging_model))
    staging_model.eval()
    with torch.inference_mode():
        logits, _ = staging_model(x_data, mask)
        probabilities = torch.softmax(logits, dim=-1)[0].cpu().numpy()

    result = pd.DataFrame(probabilities, columns=PROBABILITY_COLUMNS)
    result.insert(0, "stop_seconds", (np.arange(len(result)) + 1) * 5.0)
    result.insert(0, "start_seconds", np.arange(len(result)) * 5.0)
    result["predicted_stage_number"] = probabilities.argmax(axis=1)
    result["predicted_stage"] = result["predicted_stage_number"].map(
        dict(enumerate(PROBABILITY_COLUMNS))
    )
    result.to_csv(probability_path, index=False)

    metadata = {
        "sleepfm_commit": SLEEPFM_COMMIT,
        "device": str(runtime_device),
        "source_edf": str(edf_path),
        "prepared_hdf5": str(hdf5_path),
        "native_resolution_seconds": 5,
        "selected_channels": selected_channels,
    }
    (output_dir / f"{stem}_sleepfm_metadata.json").write_text(json.dumps(metadata, indent=2))
    return result


def head_features(prepared_hdf5_path: str | Path, sleepfm_root: str | Path | None = None) -> dict[str, np.ndarray]:
    """Frozen SleepFM up to the staging head's final linear layer, for one prepared HDF5.

    Runs the same encoder and staging-head pass as `infer_hypnodensity` (CPU,
    official checkpoints) and returns the BiLSTM output that feeds `fc`
    (`features`, steps x 128, one row per 5-s step) together with the head's
    softmax probabilities (`probabilities`) for verification.
    """
    import torch
    from torch.utils.data import DataLoader

    sleepfm_root = Path(sleepfm_root or _default_sleepfm_root())
    _, dataset_class, collate_fn, base_model_class, staging_model_class = _load_official_modules(sleepfm_root)
    checkpoints = sleepfm_root / "sleepfm" / "checkpoints"
    base_config = _load_json(checkpoints / "model_base" / "config.json")
    base_config["num_workers"] = 1
    channel_groups = _load_json(sleepfm_root / "sleepfm" / "configs" / "channel_groups.json")
    base_model = base_model_class(
        base_config["in_channels"], base_config["patch_size"], base_config["embed_dim"],
        base_config["num_heads"], base_config["num_layers"], pooling_head=base_config["pooling_head"], dropout=0.0,
    )
    base_model.load_state_dict(_state_dict_for_model(torch.load(checkpoints / "model_base" / "best.pt", map_location="cpu"), base_model))
    base_model.eval()
    dataset = dataset_class(base_config, channel_groups, hdf5_paths=[str(prepared_hdf5_path)], split="test")
    loader = DataLoader(dataset, batch_size=1, shuffle=False, num_workers=0, collate_fn=collate_fn)
    chunks: dict[str, list[tuple[int, np.ndarray]]] = defaultdict(list)
    with torch.inference_mode():
        for batch_data, mask_list, _, _, chunk_starts in loader:
            for modality, values, mask in zip(base_config["modality_types"], batch_data, mask_list):
                embeddings = base_model(values.float(), mask.bool())[1].numpy()
                for start, embedding in zip(chunk_starts, embeddings):
                    chunks[modality].append((int(start), embedding))
    sequences = [np.concatenate([a for _, a in sorted(chunks[m], key=lambda item: item[0])], axis=0)
                 for m in base_config["modality_types"]]
    steps = min(map(len, sequences))
    x = torch.from_numpy(np.stack([s[:steps] for s in sequences])).unsqueeze(0).float()
    params = _load_json(checkpoints / "model_sleep_staging" / "config.json")["model_params"]
    head = staging_model_class(**params)
    head.load_state_dict(_state_dict_for_model(torch.load(checkpoints / "model_sleep_staging" / "best.pth", map_location="cpu"), head))
    head.eval()
    captured = {}
    hook = head.lstm.register_forward_hook(lambda module, inputs, output: captured.__setitem__("lstm", output[0]))
    with torch.inference_mode():
        logits, _ = head(x, torch.zeros((1, x.shape[1], steps), dtype=torch.bool))
    hook.remove()
    return {"features": captured["lstm"][0].numpy(), "probabilities": torch.softmax(logits, dim=-1)[0].numpy(),
            "fc_weight": head.fc.weight.detach().numpy(), "fc_bias": head.fc.bias.detach().numpy()}
