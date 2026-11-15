#!/usr/bin/env bash
#
# Cut a new release of the isoTar-v2 backend.
#
#   scripts/release.sh <major|minor|patch|X.Y.Z>
#
# Bumps the VERSION file, prepends a CHANGELOG entry (with the commit subjects
# since the last tag), commits as "chore(release): vX.Y.Z" and creates an
# annotated git tag. Prints the docker build/push and git push commands to run
# next. Dependency-free (bash + git only).
#
# Run it; do not `source` it. Sourced, its `set -e` and `exit` would apply to
# your interactive shell and close the terminal on the first error.
_sourced=0
if [ -n "${ZSH_VERSION:-}" ]; then
    case "${ZSH_EVAL_CONTEXT:-}" in *:file*) _sourced=1 ;; esac
elif [ -n "${BASH_VERSION:-}" ]; then
    [ "${BASH_SOURCE[0]}" != "$0" ] && _sourced=1
fi
if [ "$_sourced" = 1 ]; then
    echo "error: run this script, do not source it:  scripts/release.sh <major|minor|patch|X.Y.Z>" >&2
    unset _sourced
    return 1
fi
unset _sourced

set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

VERSION_FILE="VERSION"
CHANGELOG="CHANGELOG.md"
IMAGE="frankfeng78/isotar-v2"

die() { echo "error: $*" >&2; exit 1; }

[ -f "$VERSION_FILE" ] || die "$VERSION_FILE not found"
[ $# -eq 1 ] || die "usage: scripts/release.sh <major|minor|patch|X.Y.Z>"

# Refuse to release with uncommitted changes to TRACKED files, so the release
# commit is clean and reviewable and the tag describes what is actually in it.
#
# Untracked files are deliberately allowed (-uno). The release commit adds only
# VERSION and CHANGELOG by name, so an untracked file can never leak into it,
# and local run artifacts (out/, work/, scratch tarballs) would otherwise block
# every release.
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
    die "tracked files have uncommitted changes; commit or stash them first
$(git status --short --untracked-files=no)"
fi

current="$(tr -d '[:space:]' < "$VERSION_FILE")"
echo "$current" | grep -Eq '^[0-9]+\.[0-9]+\.[0-9]+$' || die "current VERSION '$current' is not X.Y.Z"
IFS='.' read -r major minor patch <<EOF
$current
EOF

case "$1" in
    major) new="$((major + 1)).0.0" ;;
    minor) new="${major}.$((minor + 1)).0" ;;
    patch) new="${major}.${minor}.$((patch + 1))" ;;
    [0-9]*.[0-9]*.[0-9]*)
        echo "$1" | grep -Eq '^[0-9]+\.[0-9]+\.[0-9]+$' || die "invalid version '$1'"
        new="$1" ;;
    *) die "bump must be major|minor|patch or an explicit X.Y.Z" ;;
esac

tag="v${new}"
git rev-parse "$tag" >/dev/null 2>&1 && die "tag $tag already exists"

echo "Releasing $current -> $new"

# --- collect changes since the last tag ---
last_tag="$(git describe --tags --abbrev=0 2>/dev/null || true)"
if [ -n "$last_tag" ]; then
    range="${last_tag}..HEAD"
else
    range=""
fi
# shellcheck disable=SC2086
notes="$(git log $range --no-merges --pretty=format:'- %s' 2>/dev/null || true)"
[ -n "$notes" ] || notes="- (no changes recorded since ${last_tag:-repo start})"

# --- build the new CHANGELOG in a temp file first ---
# Nothing in the working tree is touched until this succeeds, so a failure here
# leaves no half-done release behind.
#
# The entry is inserted immediately before the first existing "## " release
# heading, or appended to the end if there is none yet. The notes are passed to
# printf as an argument, never embedded in an awk program or `awk -v` string:
# macOS (BSD) awk rejects a -v value that contains a newline, which is every
# release with more than one commit.
date="$(date +%Y-%m-%d)"
tmp="$(mktemp)"
trap 'rm -f "$tmp"' EXIT
if grep -q '^## ' "$CHANGELOG"; then
    {
        awk '/^## / { exit } { print }' "$CHANGELOG"
        printf '## %s - %s\n\n%s\n\n' "$new" "$date" "$notes"
        awk 'found { print; next } /^## / { found = 1; print }' "$CHANGELOG"
    } > "$tmp"
else
    {
        cat "$CHANGELOG"
        printf '\n## %s - %s\n\n%s\n' "$new" "$date" "$notes"
    } > "$tmp"
fi

# --- update VERSION + CHANGELOG, commit, tag ---
printf '%s\n' "$new" > "$VERSION_FILE"
cp "$tmp" "$CHANGELOG"

git add "$VERSION_FILE" "$CHANGELOG"
git commit -m "chore(release): ${tag}"
git tag -a "$tag" -m "$tag"

cat <<MSG

Released ${tag}.

Next steps:
  docker build -t ${IMAGE}:${new} --build-arg VERSION=${new} -f Dockerfile .
  docker push ${IMAGE}:${new}
  git push --follow-tags

Then bump the image tag in the deployment's docker-compose.yml to ${new}.
MSG
