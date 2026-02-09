import json
import hashlib
import os
from datetime import datetime, timezone

def generate_production_order(worldline_path):
    if not os.path.exists(worldline_path):
        print(f"❌ ERROR: WorldLine file not found: {worldline_path}")
        return

    with open(worldline_path, "r") as f:
        data = json.load(f)
    
    # GAIP SCHEMA COMPATIBILITY LAYER
    # Handle both Flat (Legacy) and Nested (GAIP Spec) structures
    if "header" in data:
        # GAIP Spec Compliant
        wl_id = data["header"]["worldline_id"]
        t_anchor = data["header"]["temporal_anchor"]
        cand_id = data["physical_realization"]["candidate_id"]
        affinity = data["physical_realization"]["binding_affinity_kj"]
        coherence = data.get("physics_snapshot", {}).get("coherence", 0.95) # Fallback if missing
        coords = data.get("physics_snapshot", {}).get("coordinate_35d", [])
    else:
        # Legacy / Flat Structure
        wl_id = data["worldline_id"]
        t_anchor = data["temporal_anchor"]
        cand_id = data["physical_realization"]["candidate_id"]
        affinity = data["physical_realization"]["binding_affinity_kj"]
        coherence = data.get("coherence_amplitude", 0.0)
        coords = data.get("coordinate_35d", [])

    # Generate a unique Batch ID based on the World-Line Fingerprint
    batch_seed = f"{wl_id}-{t_anchor}"
    batch_id = f"BATCH-{hashlib.sha256(batch_seed.encode()).hexdigest()[:8].upper()}"
    
    order = {
        "header": {
            "order_id": f"ORD-{batch_id}",
            "status": "APPROVED_FOR_SYNTHESIS",
            "compliance_protocol": "FDA-GAIP-2026-PVR",
            "timestamp": datetime.now(timezone.utc).isoformat()
        },
        "specifications": {
            "candidate_id": cand_id,
            "target_affinity_kj": affinity,
            "max_toxicity_index": 0.0200,
            "resonance_amplitude": coherence
        },
        "coordinate_lock": coords[:10] if coords else [],
        "signature": hashlib.sha256(json.dumps(data, sort_keys=True).encode()).hexdigest()
    }
    
    # Ensure directory exists relative to root
    output_dir = os.path.join("PX_Warehouse", "Orders")
    os.makedirs(output_dir, exist_ok=True)
    
    order_path = os.path.join(output_dir, f"{batch_id}.json")
    
    with open(order_path, "w") as f:
        json.dump(order, f, indent=2)
        
    print(f"✅ BATCH ORDER GENERATED: {order_path}")
    return order_path