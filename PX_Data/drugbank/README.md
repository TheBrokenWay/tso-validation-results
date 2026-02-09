# DrugBank data for intake harvester

DrugBank is **license-restricted**. You must download the data yourself and place it here.

1. **Create an account:** https://go.drugbank.com/public_users/sign_up  
2. **Download:** https://go.drugbank.com/releases — choose CSV (e.g. "drugs" or "structures" export).  
3. **Place your file** as `drugbank.csv` in this folder, **or** place the export elsewhere and set `intake_policy.json` → `sources.DrugBank.local_file` to that path (e.g. `PX_Data/drugbank/drugbank_export.csv`).

## Required CSV columns for the harvester

The intake harvester expects a CSV with at least these columns (names can vary; the converter maps them):

| Our column    | DrugBank / alternative names |
|---------------|------------------------------|
| drugbank_id   | drugbank_id, id              |
| name          | name, title                  |
| status        | groups, state (e.g. "approved", "investigational") |
| smiles        | smiles, canonical_smiles     |
| atc_codes     | atc_codes (optional)        |
| targets       | targets, polypeptides (optional) |

If your export uses different column names, run the converter:

```bash
python PX_Warehouse/Operations/scripts/convert_drugbank_to_intake.py path/to/your_drugbank_export.csv
```

That writes `drugbank.csv` in this folder in the format the harvester expects.
