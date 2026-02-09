# PX_WAREHOUSE_REGENERATION — FINAL EXECUTION SEQUENCE

The definitive instructions for regenerating every dossier and retesting the entire warehouse.

---

## 1. ARCHIVE OLD DOSSIERS (DO NOT DELETE)

The IDE must **move**, not delete, the old dossiers.

**PowerShell:**
```powershell
mkdir PX_Warehouse\Archive_Primary
mkdir PX_Warehouse\Archive_Novel

move PX_Warehouse\Prv_Dossiers PX_Warehouse\Archive_Primary
move PX_Warehouse\Novel_Dossiers PX_Warehouse\Archive_Novel
```

### Recreate empty folders with tier subdirectories

**PowerShell:**
```powershell
mkdir PX_Warehouse\Prv_Dossiers
mkdir PX_Warehouse\Prv_Dossiers\Diamond
mkdir PX_Warehouse\Prv_Dossiers\Gold
mkdir PX_Warehouse\Prv_Dossiers\Silver
mkdir PX_Warehouse\Prv_Dossiers\Bronze

mkdir PX_Warehouse\Novel_Dossiers
mkdir PX_Warehouse\Novel_Dossiers\Diamond
mkdir PX_Warehouse\Novel_Dossiers\Gold
mkdir PX_Warehouse\Novel_Dossiers\Silver
mkdir PX_Warehouse\Novel_Dossiers\Bronze
```

**Why:** The pipeline writes into tier-specific subfolders. If they don't exist, the engines cannot write new dossiers.

---

## 2. REGENERATE THE ENTIRE WAREHOUSE

This rebuilds all dossiers from scratch using the real engines, real physics, real governance, and real grading.

### Full lifecycle regeneration

**PowerShell:**
```powershell
just cycle
```

If `just` is not in PATH (e.g. Nix/direnv not loaded), run the equivalent from repo root:
```powershell
python PX_Executive/run_one_cycle_test.py
```

This runs: Feed Engine → Novel Engine → Repurposed Engine → Finalization Engine.

It produces: new Prv_Dossiers, new Novel_Dossiers, new Finalized_Dossiers/\<tier\>, all with valid JSON, full ADMET, PK/PD, Virtual Efficacy, Trial Engine, GradingEngine, Evidence Package, WorldLine entry, and Zeus authorization.

---

## 3. VALIDATE THE NEW WAREHOUSE

After regeneration, run:

**PowerShell:**
```powershell
python PX_Validation/tests/run_all_tests.py
python run_e2e_layers.py
```

**Both must pass.** This confirms: all dossiers valid, complete, graded, lifecycle-sorted, governance-compliant; physics, trial simulation, evidence generation; warehouse deterministic and sealed.

---

## 4. FINAL STATE (WHAT YOU END UP WITH)

| Location | Contents |
|----------|----------|
| **PX_Warehouse/Archive_Primary/** | Old Prv_Dossiers (preserved) |
| **PX_Warehouse/Archive_Novel/** | Old Novel_Dossiers (preserved) |
| **PX_Warehouse/Prv_Dossiers/\<tier\>/** | Fresh repurposed dossiers |
| **PX_Warehouse/Novel_Dossiers/\<tier\>/** | Fresh novel dossiers |
| **PX_Warehouse/Finalized_Dossiers/Diamond|Gold|Silver|Bronze/** | Fully processed dossiers |

**Governance & validation:** All tests pass. All layers OK. Warehouse sealed.

---

This is the exact sequence the IDE must execute to rewrite every dossier, regenerate the warehouse, retest every dossier, and complete the final step of Predator X installation. Nothing else is needed.
