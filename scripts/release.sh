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
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

VERSION_FILE="VERSION"
CHANGELOG="CHANGELOG.md"
IMAGE="frankfeng78/isotar-v2"

die() { echo "error: $*" >&2; exit 1; }

[ -f "$VERSION_FILE" ] || die "$VERSION_FILE not found"
[ $# -eq 1 ] || die "usage: scripts/release.sh <major|minor|patch|X.Y.Z>"

# Refuse to release a dirty tree (so the release commit is clean and reviewable).
if [ -n "$(git status --porcelain)" ]; then
    die "working tree is dirty; commit or stash changes first"
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

# --- update VERSION ---
printf '%s\n' "$new" > "$VERSION_FILE"

# --- prepend CHANGELOG entry below the top heading block ---
date="$(date +%Y-%m-%d)"
tmp="$(mktemp)"
entry="## ${new} - ${date}\n\n${notes}\n"
# Insert the new entry immediately before the first existing "## " release heading,
# or append to the end if there is none yet.
if grep -qn '^## ' "$CHANGELOG"; then
    awk -v entry="$entry" '
        !done && /^## / { printf entry "\n"; done=1 }
        { print }
    ' "$CHANGELOG" > "$tmp"
else
    cp "$CHANGELOG" "$tmp"
    printf '\n%b\n' "$entry" >> "$tmp"
fi
mv "$tmp" "$CHANGELOG"

# --- commit + tag ---
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
