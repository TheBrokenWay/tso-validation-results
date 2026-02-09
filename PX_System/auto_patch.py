import os
import re
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))

def main():
    print(">>> [OLYMPUS] AUTO-PATCHER PROTOCOL ENGAGED")
    print(f"    Target Root: {ROOT}")

    # PROTOCOL 1: PACKAGE UNIFICATION (__init__.py)
    # Recursively ensures every folder is a Python package
    init_count = 0
    for dirpath, dirnames, filenames in os.walk(ROOT):
        if "__pycache__" in dirpath or ".git" in dirpath:
            continue
            
        init_path = os.path.join(dirpath, "__init__.py")
        if not os.path.exists(init_path):
            try:
                with open(init_path, "w", encoding="utf-8") as f:
                    f.write("# OLYMPUS PACKAGE NODE\n")
                init_count += 1
                # print(f"    [FIX] +Package: {os.path.relpath(dirpath, ROOT)}")
            except Exception as e:
                print(f"    [ERR] {e}")
    
    print(f">>> [STATUS] Structure Normalized. Added {init_count} __init__ nodes.")

    # PROTOCOL 2: IMPORT HARMONIZATION
    # Scans for broken relative imports (e.g., 'from ..sentinels') and routes them to 'foundation'
    # strict_regex matches "from . " or "from .. " 
    patch_count = 0
    
    # We only patch if we are sure. For now, we fix the known "Sentinels" issue often found in dumps.
    for dirpath, _, filenames in os.walk(ROOT):
        for file in filenames:
            if not file.endswith(".py"): continue
            
            path = os.path.join(dirpath, file)
            with open(path, "r", encoding="utf-8") as f:
                code = f.read()
            
            # Specific fixes for the dump files
            original_code = code
            
            # Fix 1: Route explicit sub-module relative imports to foundation
            # Matches: from foundation import sentinels -> from foundation import sentinels
            code = re.sub(r"from\s+\.\.quint\s+", "from foundation ", code)
            
            # Fix 2: Generic relative catch-all (Conservative)
            # Only applies if standard relative import fails context
            # (Skipped to prevent over-patching in this iteration)

            if code != original_code:
                with open(path, "w", encoding="utf-8") as f:
                    f.write(code)
                patch_count += 1

    print(f">>> [STATUS] Import Logic Patched in {patch_count} files.")
    
    # PROTOCOL 3: CLASS GUARDS
    # Ensures critical classes exist even if the file was truncated
    _ensure_class("08_Laboratory/Neural_ADMET.py", "class NeuralADMET",
                  "\n\nclass NeuralADMET:\n    def predict(self, m):\n        raise NotImplementedError('NeuralADMET guard class: implementation required')\n")

def _ensure_class(rel_path, class_sig, guard_code):
    path = os.path.join(ROOT, rel_path)
    if os.path.exists(path):
        with open(path, "r", encoding="utf-8") as f:
            code = f.read()
        if class_sig not in code:
            with open(path, "a", encoding="utf-8") as f:
                f.write(guard_code)
            print(f">>> [FIX] Added guard class in {rel_path}")

if __name__ == "__main__":
    main()
