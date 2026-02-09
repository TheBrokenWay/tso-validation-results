"""
PX_Executive/run_calibration_ingest.py
THE UNIVERSAL TEACHER: Ingests ChEMBL, PubChem, ZINC, NCATS, ClinicalTrials.
Anchors the 35D Manifold with known truth from E:/foundation/PX_Data.
"""
import os
import sys
import csv
import logging
from pathlib import Path

# --- PATH FIX: FORCE ROOT INCLUSION ---
# This ensures we can find PX_Engine regardless of where the command is run
current_file = Path(__file__).resolve()
repo_root = current_file.parents[1]  # Should be E:\foundation
sys.path.insert(0, str(repo_root))   # Force it to the top of the list

# Also add the Current Working Directory to be safe
if os.getcwd() not in sys.path:
    sys.path.insert(0, os.getcwd())

# --- IMPORTS (canonical: OPE + ADMET; no Physics_Engine in this repo) ---
try:
    from PX_Engine.operations.OPE import run_ope
    from PX_Engine.operations.ADMET import run_admet
    from PX_Warehouse.WorldLine_Database import WorldLineDatabase
except ImportError as e:
    print(f"❌ CRITICAL IMPORT ERROR: {e}")
    print(f"   Debug: Repo Root calculated as: {repo_root}")
    print(f"   Debug: sys.path is: {sys.path}")
    sys.exit(1)

# --- CONFIGURATION ---
SOURCE_DIR = Path(r"E:\foundation\PX_Data")
DB = WorldLineDatabase(str(repo_root / "PX_Warehouse" / "WorldLines"))

# Setup Logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger("CALIBRATION")

def get_column_map(headers):
    """
    Dynamically maps the CSV headers to our standardized 'name' and 'smiles'.
    Returns (name_col, smiles_col) or None.
    """
    h_lower = [h.lower() for h in headers]
    
    # Map for SMILES
    if 'canonical_smiles' in h_lower:
        s_col = headers[h_lower.index('canonical_smiles')]
    elif 'smiles' in h_lower:
        s_col = headers[h_lower.index('smiles')]
    else:
        return None # No structure found

    # Map for NAME
    if 'pref_name' in h_lower:
        n_col = headers[h_lower.index('pref_name')]
    elif 'intervention_name' in h_lower:
        n_col = headers[h_lower.index('intervention_name')]
    elif 'name' in h_lower:
        n_col = headers[h_lower.index('name')]
    elif 'chembl_id' in h_lower: # Fallback for ChEMBL if name is missing
        n_col = headers[h_lower.index('chembl_id')]
    elif 'nct_id' in h_lower: # Fallback for ClinicalTrials
        n_col = headers[h_lower.index('nct_id')]
    elif 'cid' in h_lower: # Fallback for PubChem
        n_col = headers[h_lower.index('cid')]
    elif 'zinc_id' in h_lower: # Fallback for ZINC
        n_col = headers[h_lower.index('zinc_id')]
    else:
        return None

    return n_col, s_col

def ingest_standards():
    logger.info(f"🎓 STARTING MASSIVE CALIBRATION from {SOURCE_DIR}")
    logger.info(f"   Root Path: {repo_root}")
    
    if not SOURCE_DIR.exists():
        logger.error(f"❌ Source directory not found: {SOURCE_DIR}")
        return

    csv_files = list(SOURCE_DIR.glob("*.csv"))
    if not csv_files:
        logger.error("❌ No .csv files found in PX_Data.")
        return

    total_learned = 0

    for csv_file in csv_files:
        logger.info(f"📂 Processing Library: {csv_file.name}")
        
        try:
            # encoding='utf-8-sig' handles BOM if present
            with open(csv_file, 'r', encoding='utf-8-sig', errors='replace') as f:
                # Read headers first to map columns
                sample_line = f.readline()
                if not sample_line: continue
                f.seek(0)
                
                reader = csv.DictReader(f)
                if not reader.fieldnames:
                    logger.warning(f"   ⚠️ Empty headers in {csv_file.name}")
                    continue

                mapping = get_column_map(reader.fieldnames)
                if not mapping:
                    logger.warning(f"   ⚠️ Could not map columns for {csv_file.name}. Skipping.")
                    continue
                
                name_col, smiles_col = mapping
                logger.info(f"   -> Mapped: Name='{name_col}', SMILES='{smiles_col}'")

                count = 0
                for row in reader:
                    name = row.get(name_col, "Unknown")
                    smiles = row.get(smiles_col)

                    # Data Cleaning
                    if not smiles or len(str(smiles)) < 5: continue
                    if "smiles" in str(smiles).lower(): continue # Skip header rows
                    if not name or str(name).strip() == "":
                        # Try fallback IDs if primary name is empty
                        if 'chembl_id' in row: name = row['chembl_id']
                        elif 'nct_id' in row: name = row['nct_id']
                        elif 'cid' in row: name = f"CID_{row['cid']}"
                        elif 'zinc_id' in row: name = f"ZINC_{row['zinc_id']}"
                        else: continue

                    # --- THE LEARNING PROCESS (OPE + ADMET → 35D for WorldLine) ---
                    try:
                        # 1. OPE + ADMET (canonical engines; build 35D from these)
                        ope = run_ope(smiles)
                        admet = run_admet(smiles, ope)
                        tox = (admet.get("toxicity") or {})
                        toxicity_index = float(tox.get("toxicity_index", 0.02))
                        safety_margins = admet.get("safety_margins") or {}
                        safety_margin = float(safety_margins.get("safety_margin", 0.0))
                        # Build 35D coordinate for manifold memory (same shape as record_attempt)
                        coords = [
                            float(ope.get("logp", 2.0)),
                            float(ope.get("molecular_weight", 350.0)),
                            float(ope.get("tpsa", 70.0)),
                            float(ope.get("ec50", 1.0)),
                            toxicity_index,
                            float(tox.get("herg_risk", 0.2)),
                            float(tox.get("hepatotoxicity_risk", 0.2)),
                            float(tox.get("ames_mutagenicity", 0.05)),
                            70.0, 50.0, safety_margin, toxicity_index,
                        ]
                        while len(coords) < 35:
                            coords.append(0.0)
                        coords = coords[:35]
                        coherence = max(0.0, 1.0 - toxicity_index)
                        risk_level = tox.get("risk_level", "TOXICITY_GOLD")

                        # 2. Burn into Memory
                        lab_results = {
                            "status": risk_level,
                            "toxicity_index": toxicity_index,
                            "safety_margin": safety_margin,
                            "outcome": "Success",
                            "stage": "CALIBRATED_STANDARD",
                        }

                        # Sanitize ID for filename
                        safe_name = str(name).replace(" ", "_").replace("/", "-").replace("\\", "-")
                        safe_id = f"STD_{safe_name[:50]}" # Limit length
                        
                        DB.record_materialization(
                            task_id=safe_id,
                            block=coords,
                            coherence=coherence,
                            lab_results=lab_results,
                            route="CALIBRATION",
                            origin=f"Standard: {name} ({csv_file.name})"
                        )
                        
                        total_learned += 1
                        count += 1
                        if count % 100 == 0:
                            print(f"   ... Learned {count} molecules from {csv_file.name}")

                    except Exception as e:
                        # Silent fail on bad SMILES to keep speed up
                        pass

        except Exception as e:
            logger.error(f"   ❌ Critical error reading {csv_file.name}: {e}")

    logger.info(f"✅ CALIBRATION COMPLETE. Learned {total_learned} universal standards.")

if __name__ == "__main__":
    ingest_standards()