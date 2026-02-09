#!/usr/bin/env bash
# Run from repo root in WSL: bash scripts/fix_git_and_nix_wsl.sh
# Fixes: (1) Git index corruption, (2) CRLF in .git/hooks, (3) pre-commit hook path.
set -e
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

echo "[1/4] Removing Git index.lock if present..."
rm -f .git/index.lock

echo "[2/4] Rebuilding Git index from HEAD (fixes 'invalid data in index')..."
rm -f .git/index
git reset HEAD

echo "[3/4] Fixing CRLF in Git hooks (shebang must be python3, not python3\\\\r)..."
for f in .git/hooks/pre-commit .git/hooks/pre-push; do
  [ -f "$f" ] && sed -i 's/\r$//' "$f" && chmod +x "$f" && echo "  Fixed: $f"
done

echo "[4/4] Pre-commit hook: ensure it runs PX_Warehouse script from repo root..."
if [ -f .git/hooks/pre-commit ]; then
  if grep -q 'run_warehouse_simulation.py' .git/hooks/pre-commit && ! grep -q 'PX_Warehouse/run_warehouse_simulation' .git/hooks/pre-commit; then
    sed -i 's|"run_warehouse_simulation.py"|"PX_Warehouse/run_warehouse_simulation.py"|' .git/hooks/pre-commit
    echo "  Updated pre-commit to call PX_Warehouse/run_warehouse_simulation.py"
  fi
fi

echo "Done. You can now: git add flake.nix .envrc && git commit -m 'Add flake and envrc'"
echo "Then: direnv allow && direnv reload"
