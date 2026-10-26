"""Download the official Atera 1.0.4 human kidney analysis example (1.86 GB)."""
import argparse
import hashlib
from pathlib import Path
import shutil
import subprocess
import json

PREFIX = 'Atera_v1_WTA_Kidney_non_diseased_FFPE'
BASE = f'https://cf.10xgenomics.com/samples/atera/1.0.4/{PREFIX}/'
FILES = {
    'experiment.spatial': (2896, '06b21516b0cc4c69b6e48f31cffa8e291638ebd369283e60a1f0b6b4b76256d0'),
    'panel_config.json': (11439806, '98c98cc6b6731ae212949b1c077842c38ca73b16a27d5c4f47c9e136363f3e9a'),
    'cell_feature_matrix.zarr.zip': (1340432777, '62dd63ea5ad3da2989d0384e873d60781c153351176aefd14cf1ede42b3c1b7e'),
    'cells.zarr.zip': (512038495, '82618b41ece8bedc25aada992774c94d3765a7d9c404bdf0acee8a093914adee'),
}


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as src:
        while block := src.read(8*1024*1024):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--proxy', help='Optional HTTP(S) proxy URL')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    curl = shutil.which('curl')
    records = []
    for name, (size, digest) in FILES.items():
        target = args.output/name
        url = BASE+PREFIX+'_'+name
        if target.exists():
            if target.stat().st_size != size or sha256(target) != digest:
                raise RuntimeError(f'Existing file failed validation: {target}; move it aside before retrying')
        else:
            if curl is None:
                raise RuntimeError('curl is required to download the example')
            partial = target.with_suffix(target.suffix+'.partial')
            command = [curl, '--fail', '--location', '--retry', '3', '--retry-all-errors',
                       '--continue-at', '-', '--output', str(partial)]
            if args.proxy:
                command += ['--proxy', args.proxy]
            subprocess.run(command+[url], check=True)
            if partial.stat().st_size != size or sha256(partial) != digest:
                raise RuntimeError(f'Download failed size/SHA256 validation: {partial}')
            partial.replace(target)
        records.append(dict(file=name, bytes=size, sha256=digest, url=url))
        print(f'Validated {name}: {size:,} bytes', flush=True)
    (args.output/'download_manifest.json').write_text(json.dumps(records, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
