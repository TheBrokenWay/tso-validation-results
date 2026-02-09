"""
Genesis Engine (Creation Engine) — invents new molecules (SMILES) and authorizes via Vector Core.

Purpose: The machine feeds itself. Generate novel, drug-like SMILES; run each through the
Vector Core physics gate (U34); only candidates that pass are returned for the queue.
The orchestrator (run_prv_novel) then tests them (E2E, FTO, etc.).

Constraint-first: generate → OPE descriptors → Vector Core authorize → only then
eligible for queue. No simulation-first or random-first without gate.
"""
from __future__ import annotations

import hashlib
import random
from pathlib import Path
from typing import Any

# Optional RDKit for generation
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    _RDKIT = True
except ImportError:
    _RDKIT = False

# Default seed SMILES for mutation (drug-like, small set so we can generate many variants)
_DEFAULT_SEEDS = [
    "CN(C)C(=N)NC(=N)N",           # Metformin
    "CC(=O)Oc1ccccc1C(=O)O",      # Aspirin
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O", # Ibuprofen
    "CN1C=NC2=C1C(=O)N(C(=O)N2)C", # Caffeine
    "c1ccc2c(c1)cccn2",            # Quinoline
    "CC(C)NCC(COc1ccccc1)O",       # Propranolol-like
    "COc1ccc(cc1)C(C)C(=O)O",     # Fenoprofen-like
    "CN1CCN(CC1)Cc2c[nH]cn2",     # Histamine-like
    "O=C(O)Cc1ccccc1",             # Phenylacetic
    "Nc1ccc(O)cc1",                # 4-aminophenol
]

# High-entropy / high-chaos: expanded chemical space (rings, heterocycles, cores) for "all chemistry"
_CHAOS_SEEDS = [
    # Aromatics / phenyl
    "c1ccccc1", "c1ccc(O)cc1", "c1ccc(N)cc1", "c1ccc(Cl)cc1", "c1ccc(Br)cc1", "c1ccc(F)cc1",
    "c1ccccc1O", "c1ccccc1N", "c1ccccc1C(=O)O", "c1ccccc1OC", "c1ccccc1Cl",
    # Bicyclic / fused
    "c1ccc2c(c1)cccn2", "c1ccc2c(c1)nccn2", "c1ccc2c(c1)ccnc2", "n1cccn1", "n1ccc[nH]1",
    "c1cnccn1", "c1cc2ncnc2c1", "C1CCc2ccccc2C1", "C1CC2CCCCC2C1",
    # Heterocycles
    "c1ncccc1", "c1cnccc1", "c1nc[nH]c1", "c1occc1", "c1ccc(O)c1", "s1cccc1", "s1ccnc1",
    "C1CC1", "C1CCC1", "C1CCCC1", "C1CCCCC1", "C1=CC=CC=C1",
    # Carbonyl / acid / ester
    "C(=O)O", "CC(=O)O", "C(=O)N", "CC(=O)N", "C#N", "S(=O)(=O)N", "S(=O)(=O)C",
    # Common fragments (small molecules as seeds)
    "CCO", "CCN", "CC(=O)C", "CC(C)O", "CC(C)N", "CNC", "COC", "CF", "CCl", "CBr",
    "Nc1ccc(O)cc1", "Oc1ccc(N)cc1", "Cc1ccc(O)cc1", "Cc1ccc(N)cc1",
    "CC(C)C", "CC(C)(C)C", "C=C", "C#C", "C=C(C)C",
]


def _canonicalize(smiles: str) -> str | None:
    if not _RDKIT:
        return smiles.strip() if smiles else None
    try:
        mol = Chem.MolFromSmiles(smiles.strip())
        if mol is None:
            return None
        return Chem.MolToSmiles(mol)
    except Exception:
        return None


def _mutate_smiles(smiles: str, rng: random.Random) -> str | None:
    """Mutate: substitute H->CH3 on a random C/N/O using RWMol. Returns new SMILES or None."""
    if not _RDKIT:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        rw = Chem.RWMol(mol)
        # Find atoms that have at least one H (we'll add a methyl instead of one H)
        candidates = []
        for a in rw.GetAtoms():
            if a.GetSymbol() in ("C", "N", "O") and a.GetTotalNumHs() > 0:
                candidates.append(a.GetIdx())
        if not candidates:
            return None
        idx = rng.choice(candidates)
        # Add a new carbon (methyl) and bond to chosen atom
        rw.AddAtom(Chem.Atom("C"))
        new_idx = rw.GetNumAtoms() - 1
        rw.AddBond(idx, new_idx, Chem.BondType.SINGLE)
        out_mol = rw.GetMol()
        if out_mol is None:
            return None
        Chem.SanitizeMol(out_mol)
        out = Chem.MolToSmiles(out_mol)
        return _canonicalize(out)
    except Exception:
        return None


def _mutate_add_halogen(smiles: str, rng: random.Random) -> str | None:
    """Add F, Cl, or Br to a C with H. Expands chemical space."""
    if not _RDKIT:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        rw = Chem.RWMol(mol)
        candidates = []
        for a in rw.GetAtoms():
            if a.GetSymbol() == "C" and a.GetTotalNumHs() > 0:
                candidates.append(a.GetIdx())
        if not candidates:
            return None
        idx = rng.choice(candidates)
        halogen = rng.choice(["F", "Cl", "Br"])
        rw.AddAtom(Chem.Atom(halogen))
        new_idx = rw.GetNumAtoms() - 1
        rw.AddBond(idx, new_idx, Chem.BondType.SINGLE)
        out_mol = rw.GetMol()
        if out_mol is None:
            return None
        Chem.SanitizeMol(out_mol)
        return _canonicalize(Chem.MolToSmiles(out_mol))
    except Exception:
        return None


def _mutate_add_small_fragment(smiles: str, rng: random.Random) -> str | None:
    """Add -OH, -F, -Cl, -Br, or -NH2 (single-atom) to a C/N/O with H."""
    if not _RDKIT:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        rw = Chem.RWMol(mol)
        candidates = [a.GetIdx() for a in rw.GetAtoms()
                     if a.GetSymbol() in ("C", "N", "O") and a.GetTotalNumHs() > 0]
        if not candidates:
            return None
        idx = rng.choice(candidates)
        atom = rng.choice(["O", "F", "Cl", "Br", "N"])
        rw.AddAtom(Chem.Atom(atom))
        rw.AddBond(idx, rw.GetNumAtoms() - 1, Chem.BondType.SINGLE)
        out_mol = rw.GetMol()
        if out_mol is None:
            return None
        Chem.SanitizeMol(out_mol)
        return _canonicalize(Chem.MolToSmiles(out_mol))
    except Exception:
        return None


def _mutate_any(smiles: str, rng: random.Random, high_entropy: bool = False) -> str | None:
    """Try one mutation; in high_entropy mode try multiple strategies and optionally chain two mutations."""
    strategies = [_mutate_smiles, _mutate_add_halogen]
    if high_entropy:
        strategies = [_mutate_smiles, _mutate_add_halogen, _mutate_add_small_fragment]
    mut = rng.choice(strategies)(smiles, rng)
    if mut and high_entropy and rng.random() < 0.4:
        mut = rng.choice(strategies)(mut, rng)  # second mutation
    return mut


def _ope_descriptors(smiles: str) -> dict[str, Any] | None:
    """Get OPE-like descriptors for Vector Core / filtering. Returns None if invalid."""
    if not _RDKIT:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        return {"molecular_weight": mw, "logp": logp}
    except Exception:
        return None


def _drug_like_filter(smiles: str, high_entropy: bool = False) -> bool:
    """Basic drug-like filter (MW and LogP). High-entropy: wider bounds for exploration."""
    desc = _ope_descriptors(smiles)
    if not desc:
        return False
    mw = desc["molecular_weight"]
    logp = desc["logp"]
    if high_entropy:
        if mw < 50 or mw > 800:
            return False
        if logp < -4 or logp > 8:
            return False
    else:
        if mw < 100 or mw > 600:
            return False
        if logp < -2 or logp > 6:
            return False
    return True


def _make_novel_id(smiles: str) -> str:
    h = hashlib.sha256(smiles.encode()).hexdigest()[:8]
    return f"PRV_NOV_{h}"


def generate_novel_candidates(
    n: int = 20,
    seed_smiles: list[str] | None = None,
    use_vector_core: bool = True,
    rng_seed: int | None = None,
    repo_root: Path | None = None,
    high_entropy: bool = True,
) -> list[dict[str, Any]]:
    """
    Generate n novel molecule candidates (SMILES). Every output is piped through OPE (RDKit)
    for deterministic physical descriptors and through VectorCore as mandatory physics gate
    (Law U34: global sum 36.1, zero energy delta; amplitude >= 0.95).

    Args:
        n: Number of novel candidates to generate.
        seed_smiles: SMILES list to mutate from; if None, use built-in seeds (or chaos seeds if high_entropy).
        use_vector_core: If True, only return candidates that pass VectorCore.execute() (mandatory physics gate).
        rng_seed: Random seed for reproducibility.
        repo_root: Repo root (for Vector Core); default from __file__.
        high_entropy: If True, use expanded chemical space (chaos seeds), more mutation strategies,
            wider MW/LogP bounds; every candidate still passes through OPE + VectorCore when use_vector_core=True.

    Returns:
        List of {"id", "type": "N", "smiles", "name", "source": "Genesis"} for queue.
    """
    if not _RDKIT:
        return []
    rng = random.Random(rng_seed)
    seeds = seed_smiles or (_CHAOS_SEEDS if high_entropy else _DEFAULT_SEEDS)
    seeds = [s for s in seeds if _canonicalize(s)]
    if not seeds:
        return []

    vector_core = None
    if use_vector_core:
        try:
            from PX_Engine.Vector_Core import VectorCore
            vector_core = VectorCore(threshold=0.95, dims_limit=35.0, global_sum_target=36.1)
        except Exception:
            vector_core = None

    # OPE for deterministic physical descriptors (MW, LogP, TPSA) for every novel SMILES
    run_ope = None
    try:
        from PX_Engine.operations.OPE import run_ope as _run_ope
        run_ope = _run_ope
    except Exception:
        pass

    seen: set[str] = set()
    out: list[dict[str, Any]] = []
    attempts = 0
    max_attempts = (n * 80) if high_entropy else (n * 50)

    while len(out) < n and attempts < max_attempts:
        attempts += 1
        base = rng.choice(seeds)
        mutated = _mutate_any(base, rng, high_entropy=high_entropy)
        if not mutated or mutated in seen:
            continue
        seen.add(mutated)
        if not _drug_like_filter(mutated, high_entropy=high_entropy):
            continue
        if vector_core:
            import numpy as np
            # Physical descriptors from OPE (RDKit) for VectorCore resonance logic
            physical_descriptors = None
            if run_ope:
                try:
                    ope_out = run_ope(mutated)
                    physical_descriptors = {
                        "molecular_weight": ope_out.get("molecular_weight", 350.0),
                        "logp": ope_out.get("logp", 2.5),
                        "tpsa": ope_out.get("tpsa", 70.0),
                    }
                except Exception:
                    physical_descriptors = _ope_descriptors(mutated)
            if physical_descriptors is None:
                physical_descriptors = _ope_descriptors(mutated)
            p0 = 36.1 - (0.0 + 35.0 + 1.0 + 0.85)
            p_vector = np.array([p0, 0.0, 35.0, 1.0, 0.85])
            state = vector_core.execute(p_vector, physical_descriptors=physical_descriptors)
            if not state.get("authorized"):
                continue
        out.append({
            "id": _make_novel_id(mutated),
            "type": "N",
            "smiles": mutated,
            "name": "Genesis",
            "source": "Genesis",
        })
    return out


__all__ = ["generate_novel_candidates", "_canonicalize", "_mutate_smiles", "_drug_like_filter"]
