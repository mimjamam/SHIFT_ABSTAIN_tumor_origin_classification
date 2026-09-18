"""
Download files from the cBioPortal datahub GitHub repo (github.com/cBioPortal/datahub).

The datahub stores its data files (data_mutations.txt, data_clinical_sample.txt,
etc.) as Git LFS pointers. This repo's .lfsconfig points at a custom LFS server
(not github.com's own LFS storage), so a plain `git lfs pull` needs that server
reachable and the git-lfs binary installed. Neither is guaranteed in a minimal
environment, so this script re-implements the small part of the LFS protocol we
need: read the pointer file's oid/size, ask the LFS batch API for a presigned
download URL, then curl it directly.

Usage:
    python3 fetch_cbioportal_lfs.py <repo_path> <local_dest> [ref]

    repo_path  path inside the datahub repo, e.g. public/msk_impact_2017/data_mutations.txt
    local_dest local file path to write to
    ref        git ref, defaults to "master"
"""
import sys
import json
import subprocess
import os

LFS_BATCH_URL = "https://nsssw8k94d.execute-api.us-east-1.amazonaws.com/objects/batch"
RAW_BASE = "https://raw.githubusercontent.com/cBioPortal/datahub"


def get_pointer(repo_path: str, ref: str = "master"):
    """Fetch an LFS pointer file and return (oid, size), or None if the file
    at repo_path isn't an LFS pointer (e.g. missing, or a small plain-text file)."""
    url = f"{RAW_BASE}/{ref}/{repo_path}"
    out = subprocess.run(["curl", "-s", "--max-time", "30", url], capture_output=True, text=True)
    text = out.stdout
    if not text.startswith("version https://git-lfs"):
        return None
    oid = None
    size = None
    for line in text.splitlines():
        if line.startswith("oid sha256:"):
            oid = line.split(":", 1)[1].strip()
        if line.startswith("size "):
            size = int(line.split(" ", 1)[1].strip())
    if oid is None or size is None:
        return None
    return oid, size


def resolve_url(oid: str, size: int) -> str:
    """Ask the datahub's LFS batch API for a presigned download URL for one object."""
    body = json.dumps({
        "operation": "download",
        "transfers": ["basic"],
        "objects": [{"oid": oid, "size": size}],
    })
    out = subprocess.run(
        ["curl", "-s", "--max-time", "30", "-X", "POST", LFS_BATCH_URL,
         "-H", "Accept: application/vnd.git-lfs+json",
         "-H", "Content-Type: application/vnd.git-lfs+json",
         "-d", body],
        capture_output=True, text=True,
    )
    data = json.loads(out.stdout)
    return data["objects"][0]["actions"]["download"]["href"]


def fetch(repo_path: str, dest: str, ref: str = "master") -> bool:
    pointer = get_pointer(repo_path, ref)
    if pointer is None:
        print(f"SKIP (not an LFS pointer, or missing): {repo_path}")
        return False
    oid, size = pointer
    href = resolve_url(oid, size)
    dest_dir = os.path.dirname(dest)
    if dest_dir:
        os.makedirs(dest_dir, exist_ok=True)
    subprocess.run(["curl", "-s", "--max-time", "300", "-o", dest, href])
    actual = os.path.getsize(dest) if os.path.exists(dest) else 0
    ok = actual == size
    status = "OK" if ok else "MISMATCH"
    print(f"{status} {repo_path} -> {dest} ({actual}/{size} bytes)")
    return ok


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    repo_path_arg = sys.argv[1]
    dest_arg = sys.argv[2]
    ref_arg = sys.argv[3] if len(sys.argv) > 3 else "master"
    success = fetch(repo_path_arg, dest_arg, ref_arg)
    sys.exit(0 if success else 1)
