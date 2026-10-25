"""Synchronize the standalone HD skill runtime from src/ezst."""
import argparse
from pathlib import Path
import shutil

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--check', action='store_true')
args = parser.parse_args()
destination = root/'plugins/spatial-transcriptomics/skills/spatial-visium-hd'
mismatches = []
for source in sorted((root/'src/ezst').glob('*.py')):
    target = destination/source.name
    if args.check:
        if not target.exists() or source.read_bytes() != target.read_bytes():
            mismatches.append(source.name)
    else:
        shutil.copyfile(source, target)
if mismatches:
    raise SystemExit('Plugin runtime differs: ' + ', '.join(mismatches))
print('Plugin runtime is synchronized')
