"""Scoped prediction caching and GUI cache cleanup for OPTICS."""

import json
import os
from pathlib import Path
import stat
import tempfile


CACHE_DIRECTORIES = {
    "Predictions": "cached_predictions",
    "BLAST results": "cached_blastp_analysis",
    "Downloaded sequences": "cached_seqs",
    "Structures": "cached_structures",
}
ALL_CACHES = "All caches"


def load_prediction_cache(path):
    """Treat missing or malformed prediction caches as empty."""
    try:
        with open(path, encoding="utf-8") as stream:
            cache = json.load(stream)
        return cache if isinstance(cache, dict) else {}
    except (FileNotFoundError, json.JSONDecodeError, UnicodeDecodeError):
        return {}


def save_prediction_cache(path, cache):
    """Replace a cache atomically so an interrupted write cannot truncate it."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=path.parent, suffix=".tmp", delete=False
        ) as stream:
            temporary_path = Path(stream.name)
            json.dump(cache, stream, indent=4)
        os.replace(temporary_path, path)
    finally:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()


def prepare_prediction_cache(path, sequences, ignore_cache=False, delete_cached_predictions=False):
    """Optionally evict only the supplied, prepared protein sequences.

    Ignore alone performs no cache I/O. Combined with deletion it evicts the
    requested entries, while callers still skip saving fresh predictions.
    """
    if ignore_cache and not delete_cached_predictions:
        return {}, 0
    cache = load_prediction_cache(path)
    removed = 0
    if delete_cached_predictions:
        for sequence in sequences:
            if sequence in cache:
                del cache[sequence]
                removed += 1
        if removed:
            save_prediction_cache(path, cache)
    return ({} if ignore_cache else cache), removed


def _is_link(path):
    """Also reject Windows junctions/reparse points, including on Python 3.11."""
    info = path.lstat()
    return stat.S_ISLNK(info.st_mode) or bool(
        getattr(info, "st_file_attributes", 0)
        & getattr(stat, "FILE_ATTRIBUTE_REPARSE_POINT", 0)
    )


def cache_files(workspace, category):
    """Inventory only known cache directories; never traverse filesystem links."""
    categories = list(CACHE_DIRECTORIES) if category == ALL_CACHES else [category]
    if any(name not in CACHE_DIRECTORIES for name in categories):
        raise ValueError(f"Unknown cache category: {category}")
    workspace = Path(workspace).resolve()
    data = workspace / "data"
    if not data.exists():
        return []
    if _is_link(data):
        raise ValueError(f"Refusing to clear caches through a linked data directory: {data}")
    files = []
    for name in categories:
        root = data / CACHE_DIRECTORIES[name]
        if not root.exists() and not root.is_symlink():
            continue
        pending = [root]
        while pending:
            path = pending.pop()
            if _is_link(path) or not path.resolve().is_relative_to(root):
                raise ValueError(f"Refusing to clear a linked or out-of-scope cache path: {path}")
            if path.is_dir():
                pending.extend(path.iterdir())
            elif path.is_file():
                files.append(path)
    return files


def clear_cache(workspace, category):
    """Clear the selected GUI category, retaining directories for cache writers.

    No CLI exposes this operation. Validate the complete inventory before
    deleting anything, and only unlink files (never recursively delete a tree).
    """
    files = cache_files(workspace, category)
    for path in files:
        path.unlink()
    return len(files)
