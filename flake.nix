{
  description = "Predator X: Deterministic Industrial Laboratory Environment";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-23.11";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs { inherit system; };
      in
      {
        devShells.default = pkgs.mkShell {
          buildInputs = with pkgs; [

            # 1. Foundation & Build
            bazel_7
            just
            direnv

            # 2. Python Architect & Quality
            rye
            ruff
            pyright
            python311

            # 3. Governance & Data
            dvc
            python311Packages.hydra-core

            # 4. System Dependencies for RDKit/Physics
            zlib
            expat
          ];

          shellHook = ''
            echo "--- PREDATOR X: DETERMINISTIC ENVIRONMENT LOADED ---"
            echo "Status: 1.0.0 | Toolchain: Pinned | Metal: Local"
            bazel version
            rye --version
          '';
        };
      });
}
