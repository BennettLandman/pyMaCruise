#!/bin/sh
set -eu
repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
classes="$repo_root/build/java-reference/lesion-covariance-classes"
output="$repo_root/tests/data/synthetic/lesion_covariance/java_reference.json"
mkdir -p "$classes" "$(dirname -- "$output")"
javac -d "$classes" "$repo_root/tests/java_reference/LesionCovarianceReferenceHarness.java"
java -cp "$classes" LesionCovarianceReferenceHarness > "$output"
python -c 'import hashlib,json,pathlib,sys; out=pathlib.Path(sys.argv[1]); ref=json.loads(out.read_text()); ref["harness_sha256"]=hashlib.sha256(pathlib.Path(sys.argv[2]).read_bytes()).hexdigest(); ref["input_sha256"]=hashlib.sha256(pathlib.Path(sys.argv[3]).read_bytes()).hexdigest(); out.write_text(json.dumps(ref,separators=(",",":"))+"\n")' \
  "$output" "$repo_root/tests/java_reference/LesionCovarianceReferenceHarness.java" \
  "$repo_root/tests/data/synthetic/lesion_covariance/covariance_inputs.nii.gz"
