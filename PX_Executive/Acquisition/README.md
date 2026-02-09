# PX_Executive/Acquisition

Harvest and Bridge tools for trial acquisition and learning material ingestion.

**Clean Desk (PX_ENTERPRISE_LAYOUT):** Moved from repo root per Filing Protocol.

| Script | Purpose |
|--------|---------|
| `harvest_failed_trials.py` | Harvest failed trials from ClinicalTrials.gov API v2 → Olympus_Research |
| `olympus_bridge_v2.py` | NCT raw API → PX_Warehouse/Learning_Material |
| `requirements.txt` | Dependencies for harvest (requests) |

**Run from repo root:**
```bash
python PX_Executive/Acquisition/harvest_failed_trials.py [options]
python PX_Executive/Acquisition/olympus_bridge_v2.py [options]
```
