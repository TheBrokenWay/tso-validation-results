"""
One-time enrichment script: adds molecular_constraints, target_profile,
clinical_requirements, and regulatory nested sections to all existing
disease constraint JSON files.

Preserves all existing flat fields (QUINT evaluator compatibility).
Derives nested values from existing flat fields + category defaults.

Usage: python scripts/enrich_constraints.py
"""

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DISEASE_DIR = REPO_ROOT / "PX_Domain" / "PRV_Diseases"
MANIFEST_PATH = DISEASE_DIR / "manifest.json"


def load_manifest_categories():
    """Build disease_id -> category mapping from manifest."""
    with open(MANIFEST_PATH, "r", encoding="utf-8") as f:
        manifest = json.load(f)

    categories = {}
    for d in manifest.get("prv_diseases", []):
        disease_id = d["name"].lower().replace(" ", "_")
        categories[disease_id] = d.get("category", "tropical")
    return categories


# Category-specific defaults for new sections
TISSUE_TARGETS = {
    ("tropical", "virus"):    ["blood", "liver", "respiratory"],
    ("tropical", "bacteria"): ["skin", "blood", "liver"],
    ("tropical", "parasite"): ["blood", "liver", "skin"],
    ("rare_pediatric", "genetic"): ["cns", "muscle", "liver"],
    ("mcm", "virus"):        ["skin", "blood", "liver"],
    ("mcm", "bacteria"):     ["blood", "respiratory", "liver"],
    ("mcm", "physical"):     ["blood", "bone_marrow"],
}

MECHANISM_MAP = {
    "virus":    ["antiviral"],
    "bacteria": ["antibiotic"],
    "parasite": ["antiparasitic"],
    "genetic":  ["gene_therapy_modifier"],
    "physical": ["radioprotector"],
}

CLINICAL_DEFAULTS = {
    ("tropical", "virus"): {
        "route": ["injectable", "oral"],
        "patient_population": "Adults and children in endemic regions",
        "treatment_setting": "hospital",
    },
    ("tropical", "bacteria"): {
        "route": ["oral", "injectable"],
        "patient_population": "Adults and children in endemic regions",
        "treatment_setting": "field",
    },
    ("tropical", "parasite"): {
        "route": ["oral"],
        "patient_population": "Adults and children in endemic regions",
        "treatment_setting": "outpatient",
    },
    ("rare_pediatric", "genetic"): {
        "route": ["oral", "injectable"],
        "patient_population": "Pediatric patients, weight-based dosing",
        "treatment_setting": "hospital",
    },
    ("mcm", "virus"): {
        "route": ["oral", "injectable"],
        "patient_population": "All ages, emergency/stockpile setting",
        "treatment_setting": "hospital",
    },
    ("mcm", "bacteria"): {
        "route": ["oral", "injectable"],
        "patient_population": "All ages, emergency/stockpile setting",
        "treatment_setting": "hospital",
    },
    ("mcm", "physical"): {
        "route": ["injectable", "oral"],
        "patient_population": "All ages, radiation exposure setting",
        "treatment_setting": "hospital",
    },
}


def enrich_file(filepath, category):
    """Add nested sections to a single constraint file."""
    with open(filepath, "r", encoding="utf-8") as f:
        data = json.load(f)

    pathogen_type = data.get("pathogen_type", "virus")
    key = (category, pathogen_type)

    # Build molecular_constraints from existing flat fields
    data["category"] = category
    data["molecular_constraints"] = {
        "mw_min": 150,
        "mw_max": data.get("molecular_weight_max", 500.0),
        "logp_min": data.get("logp_min", -1.0),
        "logp_max": data.get("logp_max", 5.0),
        "hbd_max": data.get("hbd_max", 5),
        "hba_max": data.get("hba_max", 10),
        "tpsa_max": 140,
        "rotatable_bonds_max": 10,
        "toxicity_max": data.get("toxicity_threshold", 0.0200),
    }

    # Build target_profile from existing fields
    data["target_profile"] = {
        "pathogen_type": pathogen_type,
        "mechanism_classes": MECHANISM_MAP.get(pathogen_type, ["unknown"]),
        "tissue_targets": TISSUE_TARGETS.get(key, ["blood", "liver"]),
    }

    # Build clinical_requirements from category defaults
    clinical = CLINICAL_DEFAULTS.get(key, {
        "route": ["oral", "injectable"],
        "patient_population": "Adults and children",
        "treatment_setting": "hospital",
    })
    data["clinical_requirements"] = clinical

    # Build regulatory section
    data["regulatory"] = {
        "fda_pathway": "PRV",
        "guidance_document": None,
        "orphan_eligible": True,
    }

    # Bump version
    data["constraint_version"] = "2.0.0"

    with open(filepath, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
        f.write("\n")

    return data["disease_id"]


def main():
    categories = load_manifest_categories()

    # Extra files not in manifest get default categories
    categories["cholera"] = "tropical"
    categories["yellow_fever"] = "tropical"

    enriched = []
    for filepath in sorted(DISEASE_DIR.glob("*.json")):
        if filepath.stem == "manifest":
            continue

        disease_id = filepath.stem
        category = categories.get(disease_id, "tropical")
        enrich_file(filepath, category)
        enriched.append(disease_id)
        print(f"  Enriched: {disease_id} ({category})")

    print(f"\nTotal enriched: {len(enriched)} files")


if __name__ == "__main__":
    main()
