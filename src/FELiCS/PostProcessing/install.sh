#!/usr/bin/env bash
#
# install.sh -- put the FELiCS FSI post-processing tools on your PATH, so you can
#               run them from ANY case folder while editing them in ONE place.
#
#   cd <FELiCS_root>/FELiCS/PostProcessing
#   ./install.sh
#
# This creates SYMLINKS (not copies) in ~/.local/bin:
#
#     fsi-modes   -> .../FELiCS/PostProcessing/plot_fsi_modes.py
#     fsi-growth  -> .../FELiCS/PostProcessing/plot_fsi_growth.py
#
# Because they are symlinks, editing the scripts here immediately affects every
# case folder -- nothing to re-copy, nothing to re-install.
#
# Uninstall:  rm ~/.local/bin/fsi-modes ~/.local/bin/fsi-growth
#
set -euo pipefail

SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN_DIR="${1:-$HOME/.local/bin}"

echo "Source : $SRC_DIR"
echo "Target : $BIN_DIR"
mkdir -p "$BIN_DIR"

chmod +x "$SRC_DIR/plot_fsi_modes.py" "$SRC_DIR/plot_fsi_growth.py"

ln -sf "$SRC_DIR/plot_fsi_modes.py"  "$BIN_DIR/fsi-modes"
ln -sf "$SRC_DIR/plot_fsi_growth.py" "$BIN_DIR/fsi-growth"

echo
echo "Installed:"
echo "  fsi-modes  -> $SRC_DIR/plot_fsi_modes.py"
echo "  fsi-growth -> $SRC_DIR/plot_fsi_growth.py"

# is BIN_DIR on PATH?
case ":$PATH:" in
  *":$BIN_DIR:"*)
    echo
    echo "$BIN_DIR is already on your PATH. You're done."
    ;;
  *)
    echo
    echo "NOTE: $BIN_DIR is NOT on your PATH. Add this to ~/.bashrc:"
    echo
    echo "    export PATH=\"\$HOME/.local/bin:\$PATH\""
    echo
    echo "then:  source ~/.bashrc"
    ;;
esac

echo
echo "Usage from any case folder:"
echo "    fsi-modes  fsi_structural_modes.csv --n-modes 5"
echo "    fsi-growth fsi_structural_modes.csv --t-end 150 --degrees"
