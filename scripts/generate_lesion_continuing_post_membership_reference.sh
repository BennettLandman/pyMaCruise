#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/.." && pwd)"
source_file="$repo_root/tests/java_reference/LesionContinuingPostMembershipReferenceHarness.java"
output_file="$repo_root/tests/data/synthetic/lesion_continuing_post_membership/java_reference.json"
build_dir="$(mktemp -d)"
trap 'rm -rf "$build_dir"' EXIT

javac -d "$build_dir" "$source_file"
java -cp "$build_dir" LesionContinuingPostMembershipReferenceHarness \
  | python3 -c 'import hashlib,json,pathlib,sys; p=pathlib.Path(sys.argv[1]); d=json.load(sys.stdin); d["harness_sha256"]=hashlib.sha256(pathlib.Path(sys.argv[2]).read_bytes()).hexdigest(); p.parent.mkdir(parents=True,exist_ok=True); p.write_text(json.dumps(d,sort_keys=True,separators=(",",":"))+"\n")' \
    "$output_file" "$source_file"
