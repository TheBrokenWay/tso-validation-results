# DATA ENGINE: SPONSOR-BASED COMPANY FILING

**Constitutional rule for discontinued/withdrawn assets from designated sponsors.**

Authority: PX_Constitution. Applies to all writers placing assets under `PX_Warehouse/CommercialAssets/`.

---

## Rule (plain language)

```
IF asset.sponsor in [Pfizer, Janssen, Johnson & Johnson, Eli Lilly, Merck, Novartis, Roche, GSK, Sanofi, AstraZeneca, Takeda, Bayer, BMS, Amgen]
   AND asset.status in ["discontinued", "withdrawn", "terminated", "failed_phase"]
THEN
   file under PX_Warehouse/CommercialAssets/<Company>/<AssetID>/

ELSE
   file under Diamond / Gold / Silver / Bronze according to commercial score (normal tiered process).
```

- **Company folders** (supplement tiered structure; do not replace it):
  - `Pfizer/` — Pfizer
  - `Janssen_JnJ/` — Janssen or Johnson & Johnson
  - `Eli_Lilly/` — Eli Lilly
  - `Merck/` — Merck
  - `Novartis/` — Novartis
  - `Roche/` — Roche
  - `GSK/` — GSK (GlaxoSmithKline)
  - `Sanofi/` — Sanofi
  - `AstraZeneca/` — AstraZeneca
  - `Takeda/` — Takeda
  - `Bayer/` — Bayer
  - `BMS/` — Bristol-Myers Squibb
  - `Amgen/` — Amgen

- **Status values** that trigger company filing: `discontinued`, `withdrawn`, `terminated`, `failed_phase`.

- **All other assets** (any sponsor, any other status) use normal tiered routing: Diamond, Gold, Silver, Bronze by commercial score.

---

## Reference implementation

- **PX_Constitution/Filing_Rules.py**:
  - `get_filing_location(asset_metadata)` → `("company", "Pfizer")` or `("tier", "Gold")`
  - `canonical_asset_path(repo_root, tier_or_company, asset_id)` supports both tier and company folder names.
