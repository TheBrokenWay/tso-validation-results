import json
import hashlib
from datetime import datetime, timezone

def generate_production_order(worldline_path):
    with open(worldline_path, "r") as f:
        data = json.load(f)
    
    # Handle both old and new worldline formats
    worldline_id = data.get('worldline_id') or data.get('header', {}).get('worldline_id', 'UNKNOWN')
    temporal_anchor = data.get('temporal_anchor') or data.get('header', {}).get('timestamp', 'UNKNOWN')
    coordinate_35d = data.get('coordinate_35d') or data.get('physics_snapshot', {}).get('coordinate_35d', [0]*35)
    coherence_amplitude = data.get('coherence_amplitude') or data.get('physics_snapshot', {}).get('coherence', 0.0)
    
    # Generate a unique Batch ID based on the World-Line Fingerprint
    batch_seed = f"{worldline_id}-{temporal_anchor}"
    batch_id = f"BATCH-{hashlib.sha256(batch_seed.encode()).hexdigest()[:8].upper()}"
    
    # Validation step: re-verify toxicity_index against the 0.0200 limit
    tox = data.get("physical_realization", {}).get("toxicity_index", 1.0)
    if tox >= 0.0200:
        raise ValueError(f"CONSTITUTIONAL VIOLATION: Candidate {worldline_id} toxicity ({tox}) exceeds Gold Tier limit (0.0200).")

    order = {
        "header": {
            "order_id": f"ORD-{batch_id}",
            "status": "APPROVED_FOR_SYNTHESIS",
            "compliance_protocol": "FDA-GAIP-2026-PVR",
            "timestamp": datetime.now(timezone.utc).isoformat()
        },
        "specifications": {
            "candidate_id": data["physical_realization"]["candidate_id"],
            "target_affinity_kj": data["physical_realization"]["binding_affinity_kj"],
            "max_toxicity_index": 0.0200,
            "resonance_amplitude": coherence_amplitude
        },
        "coordinate_lock": coordinate_35d[:10],
        "signature": hashlib.sha256(json.dumps(data, sort_keys=True).encode()).hexdigest()
    }
    
    order_path = f"PX_Warehouse/Orders/{batch_id}.json"
    import os
    os.makedirs("PX_Warehouse/Orders", exist_ok=True)
    with open(order_path, "w") as f:
        json.dump(order, f, indent=2)
    
    return order_path

if __name__ == "__main__":
    # Pull the ALPHA worldline we created earlier
    path = "PX_Warehouse/WorldLines/ALPHA-GLP1-SYNTH.worldline"
    final_path = generate_production_order(path)
    print(f">>> [PRODUCTION] Order Generated: {final_path}")
