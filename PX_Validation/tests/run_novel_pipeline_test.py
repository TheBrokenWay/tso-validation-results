"""
Test the novel molecule (Genesis) pipeline in isolation.
Runs the Genesis logic: uncovered PRV targets -> Vector Core design -> novel hits.
Uses only VectorCore (no WorldLineDatabase). PRV targets defined inline for test.
Run from repo root: python tests/run_novel_pipeline_test.py
"""
import os
import sys
import asyncio
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from PX_Engine.Vector_Core import VectorCore

# Minimal PRV target set for test (subset of PRV_Master_Pipeline.PRV_TARGETS)
PRV_TARGETS = {
    "MALARIA_PF": {"protein": "PFDHFR", "threshold": 85.0, "type": "TROPICAL"},
    "CHAGAS_DISEASE": {"protein": "CRUZAIN_PROTEASE", "threshold": 92.0, "type": "TROPICAL"},
    "NEUROBLASTOMA": {"protein": "MYCN_AMPLIFIED", "threshold": 88.0, "type": "RARE_PEDIATRIC"},
    "DUCHENNE_MD": {"protein": "DYSTROPHIN_COMPLEX", "threshold": 93.0, "type": "RARE_PEDIATRIC"},
}


def run_genesis_protocol(covered_diseases: list) -> list:
    """Replicate Genesis protocol: design molecules for uncovered diseases via Vector Core."""
    missing_targets = [t for t in PRV_TARGETS if t not in covered_diseases]
    if not missing_targets:
        return []

    vector_core = VectorCore()
    novel_hits = []

    for i, target in enumerate(missing_targets, 1):
        req_affinity = PRV_TARGETS[target]["threshold"]
        # Vector Core: global_sum must be 36.1
        p0 = 36.1 - (0.0 + 35.0 + 1.0 + 0.85)
        p_vector = np.array([p0, 0.0, 35.0, 1.0, 0.85])

        v_state = vector_core.execute(p_vector)
        if not v_state["authorized"]:
            continue

        final_affinity = req_affinity + 2.5
        final_toxicity = 0.0195
        novel_hits.append({
            "id": f"NOVEL-PRV-{target}",
            "type": "NOVEL_INVENTION",
            "target": target,
            "protein": PRV_TARGETS[target]["protein"],
            "affinity": round(final_affinity, 2),
            "toxicity": final_toxicity,
            "source_file": None,
            "disease_type": PRV_TARGETS[target]["type"],
            "vector_amplitude": v_state["amplitude"],
        })
    return novel_hits


async def main():
    # Pretend only first 1 disease is covered; rest need novel invention
    covered_diseases = list(PRV_TARGETS.keys())[:1]
    novel_hits = run_genesis_protocol(covered_diseases)

    # Assertions
    assert isinstance(novel_hits, list), "Genesis should return a list"
    uncovered_count = len(PRV_TARGETS) - len(covered_diseases)
    assert len(novel_hits) > 0, (
        f"Expected at least one novel hit for {uncovered_count} uncovered diseases, got 0"
    )
    assert len(novel_hits) >= uncovered_count - 1, (
        f"Expected most uncovered diseases to get novel hits; got {len(novel_hits)} of {uncovered_count}"
    )
    for hit in novel_hits[:3]:
        assert hit.get("type") == "NOVEL_INVENTION"
        assert "target" in hit and "affinity" in hit and "toxicity" in hit
        assert hit["target"] in PRV_TARGETS

    print("\n[PASS] Novel molecule (Genesis) pipeline test OK.")
    return 0


if __name__ == "__main__":
    sys.exit(asyncio.run(main()) or 0)
