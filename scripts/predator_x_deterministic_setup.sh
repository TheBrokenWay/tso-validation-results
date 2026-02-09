#!/usr/bin/env bash
###############################################
# Predator X — Full Deterministic Integration #
###############################################
# Run inside WSL (Ubuntu). Step 1: from PowerShell run wsl --install
# After step 2 (Nix install), open a new terminal or run: . ~/.nix-profile/etc/profile.d/nix.sh

# --- Phase 1: Foundation (WSL2, Nix, Direnv, Just, Rye) ---

# 1. Install WSL2 (PowerShell only; run this outside WSL)
# wsl --install

# 2. Inside WSL: Install Nix
sh <(curl -L https://nixos.org/nix/install) --daemon

# 3. Enable Flakes
mkdir -p ~/.config/nix
echo "experimental-features = nix-command flakes" >> ~/.config/nix/nix.conf

# 4. Install Direnv, Just, and Rye through Nix
nix-env -iA nixpkgs.direnv nixpkgs.just nixpkgs.rye

# 5. Hook Direnv to your shell
echo 'eval "$(direnv hook bash)"' >> ~/.bashrc
source ~/.bashrc


# --- Phase 2: Governance & Code Quality (Ruff, OPA, Hydra) ---

# Install governance + static analysis tools
nix-env -iA nixpkgs.ruff nixpkgs.opa nixpkgs.python311Packages.hydra-core

# Verify
ruff --version
opa version


# --- Phase 3: Data Historian (DVC) ---

# Install DVC for Physics Map versioning
nix-env -iA nixpkgs.dvc


# --- Phase 4: Build Engine (Bazel 7) ---

# Bazel will be pinned inside your flake.nix
echo "Bazel 7 will be version-locked via flake.nix"


###############################################
# End of Predator X Full Stack Setup Script   #
###############################################
