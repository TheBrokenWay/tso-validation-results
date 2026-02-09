import hashlib
import json
import os
from datetime import datetime

class MerkleTree:
    """
    Elite forensic verification for Predator X logs.
    Implements a binary Merkle Tree to ensure log integrity.
    """
    def __init__(self, data_list):
        self.leaves = [self._hash(str(data)) for data in data_list]
        self.root = self._build_tree(self.leaves)

    def _hash(self, data):
        return hashlib.sha256(data.encode('utf-8')).hexdigest()

    def _build_tree(self, nodes):
        if not nodes:
            return None
        if len(nodes) == 1:
            return nodes[0]
        
        new_level = []
        for i in range(0, len(nodes), 2):
            left = nodes[i]
            right = nodes[i+1] if i+1 < len(nodes) else left
            new_level.append(self._hash(left + right))
        
        return self._build_tree(new_level)

    def get_root(self):
        return self.root

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

class LogForensics:
    """
    Manages log verification and distributed audit replication.
    """
    def __init__(self, log_dir=None):
        self.log_dir = log_dir or os.path.join(_REPO_ROOT, "PX_Warehouse", "logs")
        self.forensic_db = os.path.join(self.log_dir, "forensic_manifest.json")

    def seal_log(self, log_filename):
        """
        Generates a Merkle Root for a log file and seals it.
        """
        log_path = os.path.join(self.log_dir, log_filename)
        if not os.path.exists(log_path):
            return None

        with open(log_path, 'r') as f:
            lines = f.readlines()

        if not lines:
            return None

        tree = MerkleTree(lines)
        root_hash = tree.get_root()

        manifest = self._load_manifest()
        manifest[log_filename] = {
            "root_hash": root_hash,
            "line_count": len(lines),
            "sealed_at": datetime.now().isoformat(),
            "status": "VERIFIED"
        }
        self._save_manifest(manifest)
        return root_hash

    def verify_log(self, log_filename):
        """
        Verifies if the current log matches the sealed root hash.
        """
        manifest = self._load_manifest()
        if log_filename not in manifest:
            return False, "Log not sealed."

        log_path = os.path.join(self.log_dir, log_filename)
        with open(log_path, 'r') as f:
            lines = f.readlines()

        tree = MerkleTree(lines)
        current_root = tree.get_root()
        sealed_root = manifest[log_filename]["root_hash"]

        if current_root == sealed_root:
            return True, "Integrity verified."
        else:
            return False, f"HASH MISMATCH: Expected {sealed_root}, got {current_root}"

    def _load_manifest(self):
        if os.path.exists(self.forensic_db):
            with open(self.forensic_db, 'r') as f:
                return json.load(f)
        return {}

    def _save_manifest(self, manifest):
        with open(self.forensic_db, 'w') as f:
            json.dump(manifest, f, indent=4)

if __name__ == "__main__":
    forensics = LogForensics()
    # Example usage: seal the audit log
    root = forensics.seal_log("predator_x_v3_audit.log")
    print(f"Log Sealed. Merkle Root: {root}")
    
    verified, msg = forensics.verify_log("predator_x_v3_audit.log")
    print(f"Verification: {verified} - {msg}")
