"""Install the validated project-local CARD ontology into argNorm."""

import argparse
import hashlib
import os
import shutil
import tempfile
from importlib.util import find_spec
from pathlib import Path

EXPECTED_SHA256 = "0b132e8fb454b45158793fe691f2d4f4b772bb13e3896d2d1dcfd7db3a2a1318"
SOURCE = Path(__file__).resolve().parent / "data" / "aro.obo"


def checksum(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=["check", "install"])
    args = parser.parse_args()

    if not SOURCE.is_file():
        raise FileNotFoundError(f"Missing project ontology: {SOURCE}")

    if checksum(SOURCE) != EXPECTED_SHA256:
        raise ValueError("Project ontology checksum mismatch")

    spec = find_spec("argnorm")
    if spec is None or spec.origin is None:
        raise RuntimeError("argNorm is not installed in this Python environment")

    destination = Path(spec.origin).resolve().parent / "data" / "aro.obo"

    print(f"Project ontology: {SOURCE}")
    print(f"argNorm ontology: {destination}")

    if destination.is_file() and checksum(destination) == EXPECTED_SHA256:
        print("OK: argNorm already has the validated ontology.")
        return

    if args.action == "check":
        raise SystemExit(
            "ERROR: argNorm ontology differs from the validated project ontology. "
            "Run the installation script with 'install' to fix it."
        )

    if destination.is_file():
        old_checksum = checksum(destination)
        backup = destination.with_name(f"aro.obo.backup_{old_checksum[:12]}")
        if backup.exists():
            if checksum(backup) != old_checksum:
                raise ValueError(f"Backup checksum mismatch: {backup}")
        else:
            shutil.copy2(destination, backup)
        print(f"Original ontology backed up: {backup}")

    with tempfile.NamedTemporaryFile(
        dir=destination.parent, prefix=".aro.obo.", delete=False
    ) as handle:
        temporary = Path(handle.name)

    try:
        shutil.copy2(SOURCE, temporary)
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)

    if checksum(destination) != EXPECTED_SHA256:
        raise RuntimeError("Installed ontology failed checksum verification")

    print("OK: validated CARD ontology installed successfully.")


if __name__ == "__main__":
    main()
