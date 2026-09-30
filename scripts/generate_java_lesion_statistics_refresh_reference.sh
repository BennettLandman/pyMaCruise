#!/bin/sh
set -eu
repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
classes="$repo_root/build/java-reference/lesion-statistics-refresh-classes"
output="$repo_root/tests/data/synthetic/lesion_statistics_refresh/java_reference.json"
power_source="$repo_root/tests/java_reference/upstream/edu/jhmi/rad/medic/utilities/PowerTable.java"
mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" "$power_source" \
  "$repo_root/tests/java_reference/LesionStatisticsRefreshReferenceHarness.java"
java -cp "$classes" LesionStatisticsRefreshReferenceHarness > "$output"
python -c 'import hashlib,json,pathlib,sys; out=pathlib.Path(sys.argv[1]); ref=json.loads(out.read_text()); ref["harness_sha256"]=hashlib.sha256(pathlib.Path(sys.argv[2]).read_bytes()).hexdigest(); out.write_text(json.dumps(ref,separators=(",",":"))+"\n")' \
  "$output" "$repo_root/tests/java_reference/LesionStatisticsRefreshReferenceHarness.java"
