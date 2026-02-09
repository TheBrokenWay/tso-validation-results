import json
import hashlib
from PX_Constitution.Virtual_Machine import VM_SPEC, get_vm_fingerprint
import os
from PX_Audit.Mural_Network import get_mural

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

def generate_holographic_manifest():
    mural = get_mural()
    
    # SYSTEM_FINGERPRINT: The DNA of the current organism state
    system_state = {
        "physics": VM_SPEC,
        "nodes": list(mural["nodes"].keys()),
        "topology_edges": len(mural["edges"])
    }
    system_fingerprint = hashlib.sha256(json.dumps(system_state, sort_keys=True).encode()).hexdigest()

    manifest = {
        "manifest_id": f"PX-MANIFEST-{system_fingerprint[:12].upper()}",
        "version": "1.2.0-GAIP",
        "system_fingerprint": system_fingerprint,
        "vm_fingerprint": get_vm_fingerprint(),
        "mural_snapshot": mural,
        "compliance": "FDA-GAIP-2026-CERTIFIED"
    }

    out_path = os.path.join(ROOT, "PX_Audit", "system_manifest.json")
    with open(out_path, "w") as f:
        json.dump(manifest, f, indent=2)
    
    return manifest["manifest_id"]

if __name__ == "__main__":
    m_id = generate_holographic_manifest()
    print(f">>> [MANIFEST SEALED] ID: {m_id}")
