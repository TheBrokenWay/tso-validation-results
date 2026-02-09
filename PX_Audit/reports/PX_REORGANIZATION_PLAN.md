# PREDATOR X - REPOSITORY REORGANIZATION PLAN
**Date:** January 26, 2026  
**Status:** 🔴 CRITICAL ISSUES IDENTIFIED

---

## EXECUTIVE SUMMARY

Comprehensive audit identified **critical misorganization** from old numbered structure (01_Executive, 02_Audit, etc.) to new PX_ prefix. Key issues:

- **4 files with broken imports** causing import failures
- **6 stub engine files** in wrong location (foundation/engine vs PX_Engine)
- **10+ files** still referencing deleted directories
- **1 missing `__init__.py`** in PX_Laboratory
- **2 orphaned directories** with broken facades

---

## CRITICAL ISSUES (Priority 1 - Immediate Fix Required)

### 1. ❌ `PX_System/foundation/engines/__init__.py`
**Problem:** References deleted `05_Engine/engines/` directory  
**Impact:** Will cause import failures  
**Action:** DELETE entire `engines/` directory (orphaned facade)

### 2. ❌ `PX_System/foundation/Evidence_Package.py` (Line 18)
**Problem:** `sys.path.insert(0, str(_REPO_ROOT / "01_Executive"))`  
**Impact:** References deleted directory  
**Action:** Change to use `PX_System.foundation.Sovereign_Log_Chain`

### 3. ❌ `PX_System/foundation/Data_Sources.py` (Line 9)
**Problem:** `sys.path.insert(0, str(_REPO_ROOT / "05_Engine"))`  
**Impact:** References deleted directory  
**Action:** Update path or remove if unused

### 4. ❌ `PX_System/foundation/discovery/__init__.py` (Line 9)
**Problem:** `sys.path.insert(0, str(_REPO_ROOT / "05_Engine"))`  
**Impact:** References deleted directory  
**Action:** Remove or update path

---

## STRUCTURAL ISSUES (Priority 2 - Reorganization)

### Engine Files Duplication

**Current State:**
```
PX_Engine/                          (REAL ENGINES - 6 files)
├── Vector_Core.py                  ✅ Real implementation
├── Metabolism.py                   ✅ Real implementation
├── Trajectory_Predictor.py         ✅ Real implementation
├── Block_Orchestrator.py           ✅ Real implementation
├── Engine_Orchestrator.py          ✅ Real implementation
└── Stress_Test.py                  ✅ Real implementation

PX_System/foundation/engine/        (STUBS - 6 files)
├── OBE.py                          ⚠️ Stub (3 lines)
├── OCE.py                          ⚠️ Stub (3 lines)
├── OLE.py                          ⚠️ Stub (3 lines)
├── OME.py                          ⚠️ Stub (3 lines)
├── OPE.py                          ⚠️ Stub (3 lines - but core.py has QSAR constants for it)
└── OSE.py                          ⚠️ Stub (3 lines)
```

**Analysis:**
- `PX_Engine/` has REAL implementations (Vector, Metabolism, Trajectory)
- `foundation/engine/` has STUBS for different engines (Operational engines: OBE, OCE, etc.)
- These serve **different purposes** but naming is confusing

**Recommendation:**
- **Option A:** Move `foundation/engine/` stubs to `PX_Engine/operations/` for clarity
- **Option B:** Implement the operational engines properly in `PX_Engine/`
- **Option C:** Delete stubs if not needed and remove from imports

---

### Discovery Files

**Current State:**
```
PX_System/foundation/discovery/
├── candidate_discovery_engine.py   ⚠️ Stub (raises NotImplementedError)
└── __init__.py                     ❌ References deleted 05_Engine
```

**Recommendation:**
- **Option A:** Create top-level `PX_Discovery/` directory
- **Option B:** Move to `PX_Engine/discovery/` if engine-related
- **Option C:** Implement properly or remove if unused

---

## MISSING FILES

### 1. Missing `__init__.py`
```
PX_Laboratory/                      ❌ Missing __init__.py
├── Manufacturing_Manifest.py
├── Simulation_Engine.py
└── Synthetic_Expansion.py
```

**Action:** Create `PX_Laboratory/__init__.py`

---

## FILES REFERENCING DELETED DIRECTORIES

### References to `01_Executive`:
1. `PX_System/foundation/Evidence_Package.py` (line 18) ❌
2. `PX_System/foundation/__init__.py` (line 15) ⚠️

### References to `02_Audit`:
1. `PX_System/foundation/__init__.py` (line 14) ⚠️

### References to `05_Engine`:
1. `PX_System/foundation/engines/__init__.py` (line 22) ❌
2. `PX_System/foundation/Data_Sources.py` (line 9) ❌
3. `PX_System/foundation/discovery/__init__.py` (line 9) ❌
4. `PX_System/config/__init__.py` (line 8) ⚠️
5. `PX_System/Restore_All_Systems.py` (multiple) ⚠️
6. `PX_System/Run_Olympus.py` (line 17) ⚠️
7. `PX_System/update_logic.py` (multiple) ⚠️
8. `PX_System/main.py` (multiple) ⚠️

### References to `08_Laboratory`:
1. `PX_System/foundation/__init__.py` (line 16) ⚠️

---

## ORPHANED DIRECTORIES

1. **`PX_System/foundation/engines/`** - Only contains broken `__init__.py`, DELETE
2. **`PX_System/foundation/quint/brains/`** - Empty (only `__init__.py`), evaluate need

---

## REORGANIZATION PLAN

### Phase 1: Critical Fixes (Execute Now)

1. ✅ Delete `PX_System/foundation/engines/` directory
2. ✅ Fix `Evidence_Package.py` import (line 18)
3. ✅ Fix `Data_Sources.py` import (line 9)
4. ✅ Fix `discovery/__init__.py` import (line 9)
5. ✅ Create `PX_Laboratory/__init__.py`

### Phase 2: Engine Reorganization (Next)

**Proposal:** Move operational engine stubs to `PX_Engine/operations/`

```
PX_Engine/
├── operations/              NEW DIRECTORY
│   ├── __init__.py
│   ├── OBE.py              Move from foundation/engine/
│   ├── OCE.py              Move from foundation/engine/
│   ├── OLE.py              Move from foundation/engine/
│   ├── OME.py              Move from foundation/engine/
│   ├── OPE.py              Move from foundation/engine/
│   └── OSE.py              Move from foundation/engine/
├── Vector_Core.py          Keep
├── Metabolism.py           Keep
├── Trajectory_Predictor.py Keep
├── Block_Orchestrator.py   Keep
├── Engine_Orchestrator.py  Keep
└── Stress_Test.py          Keep
```

**Then delete:** `PX_System/foundation/engine/` directory

### Phase 3: Discovery Reorganization (Next)

**Proposal:** Create top-level `PX_Discovery/` directory

```
PX_Discovery/                NEW DIRECTORY
├── __init__.py
└── candidate_discovery_engine.py  Move from foundation/discovery/
```

**Then delete:** `PX_System/foundation/discovery/` directory

### Phase 4: Update Foundation Package (Next)

**Update `PX_System/foundation/__init__.py`:**
- Remove EXPORTS references to deleted directories
- Update to reference PX_Executive, PX_Audit, PX_Laboratory

### Phase 5: Code Cleanup (Final)

Update all remaining files with old references:
1. `PX_System/config/__init__.py`
2. `PX_System/Restore_All_Systems.py`
3. `PX_System/Run_Olympus.py`
4. `PX_System/update_logic.py`
5. `PX_System/main.py`

---

## FINAL STRUCTURE (Proposed)

```
e:\foundation\
├── PX_Audit/                ✅ Correct location
├── PX_Constitution/         ✅ Correct location
├── PX_Discovery/            🆕 NEW - Move from foundation/discovery/
├── PX_Engine/               ✅ Correct location
│   └── operations/          🆕 NEW - Move from foundation/engine/
├── PX_Executive/            ✅ Correct location
├── PX_Laboratory/           ✅ Correct location (add __init__.py)
├── PX_Security/             ✅ Correct location
├── PX_System/               ✅ Correct location
│   └── foundation/          ✅ Keep for core functionality
│       ├── core.py          ✅ QSAR constants
│       ├── api.py           ✅ API layer
│       ├── integrations/    ✅ Keep (net_policy, smiles_security, retry)
│       ├── quint/           ✅ Keep (evaluate brains/ subdirectory)
│       ├── ZeusLaws.py      ✅ Constitutional governance
│       ├── Sovereign_Log_Chain.py  ✅ Audit trail
│       ├── Evidence_Package.py     ✅ FDA compliance
│       ├── Emergency_Stop.py       ✅ Safety switch
│       └── [other core files]      ✅ System core
├── PX_Validation/           ✅ Correct location
└── PX_Warehouse/            ✅ Correct location
```

---

## EXECUTION CHECKLIST

### Immediate (Phase 1):
- [ ] Delete `PX_System/foundation/engines/`
- [ ] Fix `Evidence_Package.py` import
- [ ] Fix `Data_Sources.py` import
- [ ] Fix `discovery/__init__.py` import
- [ ] Create `PX_Laboratory/__init__.py`
- [ ] Test all imports with `PX_System_Test.py`

### Short-term (Phase 2-3):
- [ ] Create `PX_Engine/operations/`
- [ ] Move 6 operational engine stubs
- [ ] Update engine imports
- [ ] Create `PX_Discovery/`
- [ ] Move discovery files
- [ ] Delete old directories
- [ ] Test all imports

### Medium-term (Phase 4-5):
- [ ] Update `foundation/__init__.py`
- [ ] Update 5 files with old references
- [ ] Update documentation
- [ ] Full system test
- [ ] Update `PX_ADJUSTMENT_REPORT.md`

---

**Report Status:** 🔴 CRITICAL FIXES REQUIRED  
**Estimated Effort:** 2-3 hours for complete reorganization  
**Risk:** Medium (import failures if not done carefully)
