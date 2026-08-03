#!/usr/bin/env bash
# Rebuild PITA's RNAddG4 binary (64-bit) from the ViennaRNA-1.6 sources in tools/.
#
# WHY THIS EXISTS
#   RNAddG4 computes PITA's target-site accessibility term (dG0/dG1 -> dGopen).
#   PITA scores ddG = dGduplex - dGopen. The RNAddG4 executable originally shipped
#   in tools/ (md5 bd08d3ab...) SEGFAULTS before emitting any output on this base
#   image. Nothing checks for that: RNAddG_compute.pl doesn't test $? , run_pita()
#   discards stderr, and the Perl computes `ddG = dGduplex + (dG1 - dG0)` where the
#   missing dG0/dG1 are undef -> 0. The result is a well-formed .tab file with
#   dGopen blank and ddG silently equal to dGduplex, which makes PITA over-predict
#   ~20-40x while looking perfectly healthy.
#
#   The shipped tools/lib/libRNA.a is 32-bit and cannot link a 64-bit binary, so
#   this script builds the ViennaRNA-1.6 library from source (64-bit) first, then
#   compiles and links RNAddG4 against it.
#
# VERIFIED
#   The resulting binary produces byte-identical output to isoTar v1's working
#   RNAddG4 (md5 4901a2f2...) across test sequences, and an end-to-end PITA run
#   reproduces v1's target set exactly. Like v1's binary, it prints its results and
#   then aborts on exit -- a latent double-free in the 2006 RNAddG4.c source. That
#   abort is harmless because stdout is flushed first; the broken build differs in
#   that it dies BEFORE printing.
#
# USAGE
#   Run inside the isotar image (it has gcc 5.4 and the sources under /opt):
#     docker run --rm -v "$PWD:/repo" --entrypoint bash frankfeng78/isotar-v2:<tag> \
#       /repo/scripts/rebuild_RNAddG4.sh /repo/tools/PITA64bit/Bin/ViennaRNA/ViennaRNA-1.6/Progs/RNAddG4
#
#   $1 = optional output path. Defaults to writing next to the sources.
set -euo pipefail

VRNA="${VRNA_DIR:-/opt/PITA64bit/Bin/ViennaRNA/ViennaRNA-1.6}"
OUT="${1:-$VRNA/Progs/RNAddG4}"
BUILD="$(mktemp -d)"
trap 'rm -rf "$BUILD"' EXIT

[ -f "$VRNA/Progs/RNAddG4.c" ] || { echo "ERROR: $VRNA/Progs/RNAddG4.c not found" >&2; exit 1; }

echo "Building 64-bit libRNA from $VRNA/lib ..."
cd "$BUILD"
for f in "$VRNA"/lib/*.c; do
    gcc -O2 -I"$VRNA/H" -I"$VRNA" -c "$f" -o "$(basename "$f" .c).o"
done
ar rcs libRNA64.a ./*.o
echo "  archive: $(ar t libRNA64.a | wc -l) members"

echo "Compiling and linking RNAddG4 ..."
gcc -O2 -I"$VRNA/H" -I"$VRNA" -c "$VRNA/Progs/RNAddG4.c" -o RNAddG4.o
gcc -O2 -o "$BUILD/RNAddG4.new" RNAddG4.o libRNA64.a -lm

# Smoke test: must print a values line. A binary that prints nothing is the
# failure mode we are fixing, so treat empty output as a hard error.
echo "Smoke testing ..."
RESULT="$(echo "GGGGAAAACCCCAAAAGGGGUUUU" | "$BUILD/RNAddG4.new" -u 0 -s 5 -f 25 -t 25 2>/dev/null || true)"
if [ -z "$RESULT" ]; then
    echo "ERROR: rebuilt RNAddG4 produced no output -- refusing to install." >&2
    exit 1
fi
echo "  ok: ${RESULT:0:70}..."

install -m 755 "$BUILD/RNAddG4.new" "$OUT"
echo "Installed -> $OUT"
md5sum "$OUT"
