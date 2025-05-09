import json
from pathlib import Path
import re
import shutil
import ast  # safe way to parse string literals

# Create a new index.json file from the metadata.json files in each entry folder. Note: signatures are stringified for readability.
def rebuild_index(entries_dir=Path(r"cnap_db/patterns"), index_path=Path(r"cnap_db/index.json")):
    """
    Scan all entry folders and rebuild index.json.

    Args:
        entries_dir (str or Path): Path to folder containing all entry subfolders
        index_path (str or Path): Output path for the rebuilt index
    """
    entries_dir = Path(entries_dir)
    index = {}

    for entry_folder in sorted(entries_dir.iterdir(), key=extract_numeric_id):
        if not entry_folder.is_dir():
            continue

        metadata_file = entry_folder / "metadata.json"
        if not metadata_file.exists():
            print(f"Skipping {entry_folder.name}: no metadata.json")
            continue

        with open(metadata_file) as f:
            metadata = json.load(f)

        entry_id = entry_folder.name
        index[entry_id] = {
            "name": metadata.get("name", ""),
            "tags": metadata.get("tags", ""),
            #"counts": metadata.get("counts", ""),
            "folder": str(entry_folder.relative_to(entries_dir.parent)),
            "cna_sequence": str(metadata.get("cna_sequence", [])),
            "frequencies": str(metadata.get("frequencies", [])),
            "coordination": metadata.get("coordination", 0),
        }

    with open(index_path, "w") as f:
        json.dump(index, f, indent=2)
    print(f"Rebuilt index with {len(index)} entries at {index_path}")

# Update metadata.json files in each entry folder with data from index.json.
def sync_index_to_metadata(index_path=Path(r"cnap_db/index.json"), entries_dir=Path(r"cnap_db/patterns")):
    """
    Overwrite metadata.json files in each entry folder with updated data from index.json.

    Args:
        index_path (str or Path): Path to updated index.json
        entries_dir (str or Path): Path to entries folder
    """
    index = json.load(open(index_path))
    entries_dir = Path(entries_dir)

    for entry_id, data in index.items():
        entry_path = entries_dir / entry_id
        metadata_path = entry_path / "metadata.json"

        if not entry_path.exists():
            print(f"Entry folder {entry_id} missing, skipping.")
            continue

        if metadata_path.exists():
            old_metadata = json.load(open(metadata_path))
        else:
            old_metadata = {}

        # Merge updated fields from index
        updated_metadata = {**old_metadata, **data}

        with open(metadata_path, "w") as f:
            json.dump(updated_metadata, f, indent=2)

        print(f"Updated {metadata_path}")

# Batch edit metadata.json files in all entry folders.
def batch_edit_metadata(entries_dir=Path(r"cnap_db/patterns"), remove_fields=None, rename_fields=None, field_add_callbacks=None, field_edit_callbacks=None, entry_filter=None):
    """
    Edits metadata.json files in all entry folders.

    Args:
        entries_dir (str or Path): path to entries folder
        remove_fields (list[str]): list of field names to delete
        rename_fields (dict): keys are old field names, values are new names
        field_add_callbacks (dict): keys are field names, values are functions to call to generate the field value
        field_edit_callbacks (dict): keys are field names, values are functions to apply to the field value
    """
    entries_dir = Path(entries_dir)
    remove_fields = remove_fields or []
    rename_fields = rename_fields or {}

    for entry_folder in sorted(entries_dir.iterdir()):
        metadata_path = entry_folder / "metadata.json"
        if not metadata_path.exists():
            continue

        with open(metadata_path) as f:
            metadata = json.load(f)

        if entry_filter is not None and not entry_filter(metadata):
            soft_delete_entry(entry_folder, trash_root="cnap_db/trash")
            continue

        # Remove specified fields
        for key in remove_fields:
            metadata.pop(key, None)

        # Rename fields
        for old_key, new_key in rename_fields.items():
            if old_key in metadata:
                metadata[new_key] = metadata.pop(old_key)
        
        # Add missing fields
        if field_add_callbacks:
            for key, func in field_add_callbacks.items():
                if key not in metadata:
                    metadata[key] = func(metadata)

        # Apply per-field edit callbacks
        if field_edit_callbacks:
            for key, func in field_edit_callbacks.items():
                if key in metadata:
                    metadata[key] = func(metadata[key])

        # Save
        with open(metadata_path, "w") as f:
            json.dump(metadata, f, indent=2)

        print(f"Updated {entry_folder.name}")

# Reindex entries by renaming folders and updating metadata.json files.
def reindex_entries(entries_dir=Path(r"cnap_db/patterns"), index_path=Path(r"cnap_db/index.json"), prefix="CNAP-"):
    entries_dir = Path(entries_dir)
    all_folders = sorted(
    [f for f in entries_dir.iterdir() if f.is_dir()],
    key=extract_numeric_id
    )

    new_index = {}
    rename_map = []

    for i, old_path in enumerate(all_folders):
        if not old_path.is_dir():
            continue

        new_id = f"{prefix}{i}"
        new_path = entries_dir / new_id

        with open(old_path / "metadata.json") as f:
            metadata = json.load(f)

        metadata["id"] = new_id
        with open(old_path / "metadata.json", "w") as f:
            json.dump(metadata, f, indent=2)

        new_index[new_id] = {
            "name": metadata.get("name", ""),
            "tag": metadata.get("tag", ""),
            "source": metadata.get("source", ""),
            "folder": str(new_path.relative_to(entries_dir.parent))
        }

        if new_path != old_path:
            rename_map.append((old_path, new_path))

    # Now do the renaming after all paths are computed
    for old_path, new_path in rename_map:
        if new_path.exists():
            raise FileExistsError(f"Destination {new_path} already exists!")
        old_path.rename(new_path)

    with open(index_path, "w") as f:
        json.dump(new_index, f, indent=2)

    print(f"Reindexed {len(rename_map)} entries.")

# Splits the CNA-signature field into two separate fields: cna_sequence and frequencies. Only needed when converting from old format.
def split_cna_signature_field(entries_dir=Path(r"cnap_db/patterns")):
    entries_dir = Path(entries_dir)
    for entry_folder in sorted(entries_dir.iterdir()):
        metadata_path = entry_folder / "metadata.json"
        if not metadata_path.exists():
            continue

        with open(metadata_path) as f:
            metadata = json.load(f)

        sig_str = metadata.get("CNA-signature")
        if sig_str is None:
            continue

        try:
            parsed = ast.literal_eval(sig_str)
            metadata["cna_sequence"] = [tpl[0] for tpl in parsed]
            metadata["frequencies"] = [tpl[1] for tpl in parsed]
            del metadata["CNA-signature"]

            with open(metadata_path, "w") as f:
                json.dump(metadata, f, indent=2)
            print(f"Updated {entry_folder.name}")

        except Exception as e:
            print(f"Failed on {entry_folder.name}: {e}")

# Extract numeric ID from folder name or ID
def extract_numeric_id(folder):
    # Extract the numeric part after the last '-' (e.g., CNAP-42 → 42)
    try:
        return int(str(folder.name).split("-")[-1])
    except ValueError:
        return float('inf')  # send malformed names to the end

# Soft-deletes an entry by moving it to a trash folder
def soft_delete_entry(entry_path, trash_root=Path(r"cnap_db/trash")):
    trash_root = Path(trash_root)
    trash_root.mkdir(parents=True, exist_ok=True)
    trash_path = trash_root / entry_path.name
    if trash_path.exists():
        print(f"Entry {trash_path.name} already in trash, skipping.")
        return
    shutil.move(str(entry_path), str(trash_path))
    print(f"Moved {entry_path.name} → trash/")

# Filter function to match only patterns in counts that contain info about the run they are from
def filter_counts(counts_dict):
    pattern = re.compile(r"D\d+\([^,]+,\s*[^,)]+\)")
    return {k: v for k, v in counts_dict.items() if pattern.match(k)}

# Use this function to convert from the old database format to the new one
def convert_from_old_db():
    """
    This function looks at all entries in the db and converts them to the new format, if needed.
    It also removes any entries that are empty or have no counts.
    Note: counts are converted and filtered to only include those that contain info about the run they are from. Then the code checks their lenght after the process
    and can remove some data you may need.
    Caution: this function will overwrite the metadata.json files in the entries folder. It will also reindex the db after the conversion. It is very destructive!
    """
    rebuild_index()
    batch_edit_metadata(rename_fields={"pattern_name":"name", "pattern_id":"id", "sigs":"CNA-signatute"}, remove_fields=["sources"], field_edit_callbacks={"counts": filter_counts})
    batch_edit_metadata(entry_filter=lambda x: len(x.get("counts", {}))>0)
    reindex_entries()
    split_cna_signature_field()
    batch_edit_metadata(field_add_callbacks={"coordination": lambda x: sum(x.get("frequencies", []))})
    rebuild_index()

if __name__ == "__main__":
    sync_index_to_metadata()
    rebuild_index()

    pass


