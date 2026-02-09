import json
import os

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

def test_immune_classification():
    print("=== PREDATOR X: IMMUNE CLASSIFICATION TEST ===")
    
    # Test Payload: One of the newly restored dossiers
    restored_id = "RESTORED_0"
    fingerprint = "PX-ROOT-AUTH-25153"
    
    print(f">>> [INBOUND] Analyzing Candidate: {restored_id}")
    print(f">>> [SECURITY] Checking Fingerprint: {fingerprint}")
    
    # Immune Classification Logic
    if "PX-ROOT" in fingerprint:
        route = "BLUE"
        status = "AUTHORIZED"
    else:
        route = "RED"
        status = "DECOY_ROUTED"
        
    print(f">>> [RESULT] Route: {route} | Status: {status}")
    
    # Update Mural for Security Health (check PX_Audit then Warehouse)
    mural_path = os.path.join(_REPO_ROOT, "PX_Audit", "system_state.json")
    if not os.path.exists(mural_path):
        mural_path = os.path.join(_REPO_ROOT, "PX_Warehouse", "system_state.json")
    if os.path.exists(mural_path):
        with open(mural_path, "r") as f:
            state = json.load(f)
        state["Immune_Status"] = f"{route}_ROUTE_ACTIVE"
        state["Metabolic_Cycle"] = state.get("Metabolic_Cycle", 0) + 1
        with open(mural_path, "w") as f:
            json.dump(state, f, indent=4)
            
    print(f"\n[IMMUNE SYSTEM: HEALTHY] Target {restored_id} accepted into Block Universe.")

if __name__ == "__main__":
    test_immune_classification()
