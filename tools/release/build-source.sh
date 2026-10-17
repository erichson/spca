#!/usr/bin/env bash
# Build once, then retain these exact bytes for checks and submission.
# Usage: R_BUILD=/path/to/R bash tools/release/build-source.sh /new/output/dir
set -euo pipefail

package_dir=$(cd "$(dirname "$0")/../.." && pwd)
release_dir=${1:?Supply a new output directory outside the package checkout}
r_build=${R_BUILD:-R}
mkdir -p "$release_dir"
release_dir=$(cd "$release_dir" && pwd)
case "$release_dir/" in
  "$package_dir/"*) echo 'Use an output directory outside the checkout.' >&2; exit 1 ;;
esac
if compgen -G "$release_dir/sparsepca_*.tar.gz" > /dev/null; then
  echo 'An archive already exists here; preserve it and choose a new directory.' >&2
  exit 1
fi
if test -n "$(git -C "$package_dir" status --porcelain)"; then
  echo 'Commit the reviewed source before recording its release identity.' >&2
  exit 1
fi

git -C "$package_dir" rev-parse HEAD > "$release_dir/source-commit.txt"
git -C "$package_dir" rev-parse 'HEAD^{tree}' > "$release_dir/source-tree.txt"
git -C "$package_dir" remote get-url origin > "$release_dir/source-repository.txt"
date -u '+%Y-%m-%dT%H:%M:%SZ' > "$release_dir/build-time-utc.txt"
cd "$release_dir"
"$r_build" CMD build "$package_dir" 2>&1 | tee build.log
archives=(sparsepca_*.tar.gz)
test "${#archives[@]}" -eq 1
sha256sum "${archives[0]}" > SHA256SUMS
tar -tzvf "${archives[0]}" > archive-contents.txt
"$r_build" --vanilla --slave -e \
  'writeLines(capture.output(sessionInfo()), "build-session.txt")'
python3 - "${archives[0]}" <<'PYTHON'
import hashlib
import sys
import tarfile

with tarfile.open(sys.argv[1], "r:gz") as archive:
    files = sorted((item for item in archive.getmembers() if item.isfile()),
                   key=lambda item: item.name)
    with open("archive-files-sha256.tsv", "w", encoding="utf-8") as output:
        output.write("path\tsha256\tbytes\tmode\n")
        for item in files:
            content = archive.extractfile(item).read()
            output.write(f"{item.name}\t{hashlib.sha256(content).hexdigest()}\t"
                         f"{item.size}\t{item.mode:o}\n")
PYTHON
sha256sum --check SHA256SUMS
