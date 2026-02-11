# PX_Warehouse/Operations/Inputs

Input data for batch operations (reprocess candidates, queue, intake policy). **Filing Protocol:** Batch inputs moved here from root.

| File | Purpose |
|------|---------|
| `intake_policy.json` | Intake policy: sources (ChEMBL, DrugBank, ClinicalTrials, FDA_OrangeBook, PubChem, NCATS_Repurposing), filters, scientific_filters, output.queue_file, governance. |
| `prv_24h_queue.json` | PRV 24h queue (candidates with id, type R/N, smiles, name). Default queue file when intake_policy.output.queue_file not set. |
| `reprocess_candidates.json` | SMILES → name map for reprocess pipeline (extract scripts write here). |

Scripts that read/write these use this path when run from repo root.
