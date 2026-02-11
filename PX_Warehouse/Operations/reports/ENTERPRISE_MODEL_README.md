# PX_Warehouse — Enterprise Edition

**Single source of truth for paths:** `PX_Warehouse/warehouse_layout.py`

## Canonical structure (current)

| Zone | Purpose |
|------|--------|
| **Feeder** | Loading dock for queue files (e2e_test_queue.json, etc.). Resolve via `get_queue_path()`. |
| **Calibration_Molecules** | QA vault for reference trial dossiers; LiveRuns/run_* under it for live orchestrator output. |
| **Prv_Dossiers** | Repurposed PRV dossiers; subdirs **Diamond**, **Gold**, **Silver**, **Bronze** (flat). |
| **Novel_Dossiers** | Novel PRV dossiers; same four tiers (flat). |
| **Learning_Material** | Trial simulations and general learning output (flat). |
| **WorldLines** | Reference library; optional tier subdirs. Use `WorldLine_Database.path` or `DEFAULT_WORLDLINES_PATH`. |
| **Operations** | Scripts, **Inputs** (fallback for queue files), Control_Scripts. |

**Terminology:** Diamond, Gold, Silver, Bronze; Prv_Dossiers, Novel_Dossiers, Learning_Material. No legacy names (Commercial_Dossiers, Dossier_Final, etc.) for new code.

**ADMET risk_level (4-tier):** TOXICITY_DIAMOND | TOXICITY_GOLD | TOXICITY_SILVER | TOXICITY_FAILURE. See governance and `PX_Engine/operations/ADMET.py`.

Existing folders and historical content may not match this layout; new code and scripts **must** use `warehouse_layout` and the governed structure.
