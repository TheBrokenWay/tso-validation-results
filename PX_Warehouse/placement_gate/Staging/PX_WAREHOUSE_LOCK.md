# PX_WAREHOUSE LOCK

**Status:** LOCKED  
**As of:** 2026-02-03 (Post–Full Rebuild Phase 2)

The warehouse has been rebuilt to the canonical structure per **PX_WAREHOUSE_FULL_REBUILD_DIRECTIVE.md**.

- **Canonical structure** (§1 of directive) is enforced.
- **Any new folder** requires amendment to PX_WAREHOUSE_ENTERPRISE_MODEL.md and migration instruction.
- **Any drift** triggers governance failure and generation of WAREHOUSE_DRIFT_REPORT_<timestamp>.md.
- **Drift checks** run on every commit via the **pre-commit hook** (`.git/hooks/pre-commit`), which runs `python run_warehouse_simulation.py --enforce`. If structural drift is detected (non-canonical top-level folders), the commit is blocked.
- **After clone:** run `python scripts/install_warehouse_pre_commit_hook.py` once to install the hook. To run the check manually: `python run_warehouse_simulation.py --enforce`. See also WAREHOUSE_DRIFT_REPORT_EXECUTION_LOGIC.md.

**Latest drift report:** WAREHOUSE_DRIFT_REPORT_20260203_012324.md (PASS)
