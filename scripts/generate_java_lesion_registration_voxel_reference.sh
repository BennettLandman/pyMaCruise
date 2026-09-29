#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="${repo_root}/build/java-reference/lesion-registration-voxel"
mkdir -p "${build_dir}"
javac -d "${build_dir}" \
  "${repo_root}/tests/java_reference/LesionRegistrationVoxelReferenceHarness.java"
java -cp "${build_dir}" LesionRegistrationVoxelReferenceHarness \
  > "${repo_root}/tests/data/synthetic/lesion_registration_voxel/java_reference.json"
