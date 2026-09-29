#!/usr/bin/env bash
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/.." && pwd)"
source_file="${repo_root}/../toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic/methods/LesionToadSegmentation.java"
expected_sha="545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
actual_sha="$(shasum -a 256 "${source_file}" | awk '{print $1}')"
if [[ "${actual_sha}" != "${expected_sha}" ]]; then
  echo "Unexpected selected source hash: ${actual_sha}" >&2
  exit 1
fi
classes="$(mktemp -d)"
trap 'rm -rf "${classes}"' EXIT
javac -d "${classes}" "${repo_root}/tests/java_reference/LesionToadsBlackholeOutmaskSweepReferenceHarness.java"
java -cp "${classes}" LesionToadsBlackholeOutmaskSweepReferenceHarness > "${repo_root}/tests/data/synthetic/lesion_toads_blackhole_outmask_sweep/java_reference_blackhole_first.json"
java -cp "${classes}" LesionToadsBlackholeOutmaskSweepReferenceHarness outmask-first > "${repo_root}/tests/data/synthetic/lesion_toads_blackhole_outmask_sweep/java_reference_outmask_first.json"
