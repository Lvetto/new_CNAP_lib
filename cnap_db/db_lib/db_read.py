from pathlib import Path
import json
import numpy as np
from collections import defaultdict
import ast

"""
    
    Objects have a index member that is not used and is not well defined a lot of the time!!!!
    
"""


# A class to represent a CNAP entry.
# It exposes a dictionary-like interface to the metadata and provides methods to load data from the entry folder.
class CNAPEntry:
    def __init__(self, path: Path):
        self.path = Path(path)
        self.metadata_path = self.path / "metadata.json"

        if not self.metadata_path.exists():
            raise FileNotFoundError(f"Missing metadata: {self.metadata_path}")
        
        with open(self.metadata_path) as f:
            self._metadata = json.load(f)

        self._parse_metadata()

        self._cache = {}

    def _parse_metadata(self):
        for key in ("cna_sequence", "frequencies"):
            value = self._metadata.get(key)
            if isinstance(value, str):
                try:
                    parsed = ast.literal_eval(value)
                    if isinstance(parsed, list):
                        self._metadata[key] = parsed
                except Exception as e:
                    print(f"Failed to parse {key} in {self['id']}: {e}")

    # Dict-like access
    def __getitem__(self, key):
        return self._metadata.get(key)

    def get(self, key, default=None):
        return self._metadata.get(key, default)

    def __setitem__(self, key, value):
        self._metadata[key] = value

    def __delitem__(self, key):
        del self._metadata[key]

    def __iter__(self):
        return iter(self._metadata)

    def __len__(self):
        return len(self._metadata)

    def __contains__(self, key):
        return key in self._metadata

    # Data loading helpers
    def load_numpy(self, name: str):
        if name in self._cache:
            return self._cache[name]
        path = self.path / f"{name}.npy"
        if not path.exists():
            raise FileNotFoundError(f"Missing: {path}")
        arr = np.load(path)
        self._cache[name] = arr
        return arr

    # Save metadata back to disk
    def save_metadata(self):
        with open(self.metadata_path, "w") as f:
            json.dump(self._metadata, f, indent=2)

    # Display
    def summary(self):
        return {
            "Id": self["id"],
            "Signature": (self.get("cna_sequence", []), self.get("frequencies", [])),
            "Name": self.get("name", ""),
            "Tags": self.get("tag", ""),
        }

    def __repr__(self):
        return f"<CNAPEntry {self['id']} — {self.get('tag', '?')} — {self.get('name', '')}>"

# Load the index file for the db
def load_index(index_path=Path(r"cnap_db/index.json")):
    with open(index_path) as f:
        return json.load(f)

# Load all entries in the db
def load_entries_list(index=Path(r"cnap_db/index.json"), entries_dir=Path(r"cnap_db/patterns")):
    return [
        CNAPEntry(Path(entries_dir) / entry_id)
        for entry_id in index
        if (Path(entries_dir) / entry_id).exists()
    ]

# A class to represent the CNAP database.
# Its just a small wrapper around a list of CNAPEntry objects.
class CNAPDatabase:
    def __init__(self, index_path=Path(r"cnap_db/index.json"), entries_dir=Path(r"cnap_db/patterns")):
        self.index = load_index(index_path)
        self.entries = load_entries_list(self.index, entries_dir)

    # Overload getters and setters to allow list-like access to the entries
    def __getitem__(self, key):
        return self.entries[key]

    def __setitem__(self, key, value):
        self.entries[key] = value

    def __delitem__(self, key):
        del self.entries[key]

    def __iter__(self):
        return iter(self.entries)

    def __len__(self):
        return len(self.entries)

    def filter(self, condition):
        """
        Return a new CNAPDatabase with only entries matching the condition.
        
        Args:
            condition (Callable[[CNAPEntry], bool])
        
        Returns:
            CNAPDatabase: new filtered database
        """
        new_db = CNAPDatabase.__new__(CNAPDatabase)  # bypass __init__
        new_db.index = self.index
        new_db.entries = self.entries
        new_db.entries = [e for e in self.entries if condition(e)]
        return new_db

    def sort(self, key_func, reverse=False):
        new_db = CNAPDatabase.__new__(CNAPDatabase)
        new_db.index = self.index
        new_db.entries = self.entries
        # Maintain order using a list of tuples
        sorted_items = sorted(self.entries, key=lambda x: key_func(x), reverse=reverse)
        new_db.entries = sorted_items
        return new_db
    
    def group_by(self, key_func):
        groups = defaultdict(list)
        for entry in self.entries:
            group_key = key_func(entry)
            groups[group_key].append(entry)

        dbs = []
        for _, group_entries in groups.items():
            new_db = CNAPDatabase.__new__(CNAPDatabase)
            new_db.index = self.index
            new_db.entries = group_entries
            dbs.append(new_db)

        return dbs

    def map(self, func):
        return [func(entry) for entry in self.entries]

    def save(self):
        """
        Save metadata of all entries in the current DB to disk.
        Only calls the entry's save_metadata() method.
        """
        for entry in self.entries:
            try:
                entry.save_metadata()
            except Exception as e:
                print(f"Failed to save {entry['id']}: {e}")



if __name__ == "__main__":
    db = CNAPDatabase()
    db = db.filter(lambda e: [4,2,1] in e["cna_sequence"])

    for entry in db:
        print(entry.summary())


    pass
