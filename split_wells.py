"""Write one SpatialData zarr per well, zipped, because the CDN can't serve a single 111 GB file.

Usage: python split_wells.py <full.zarr> <out dir>
"""

import shutil
import subprocess
import sys
from pathlib import Path

import spatialdata as sd

import prep


def main(full: str, out: str) -> None:
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    sdata = sd.read_zarr(full)
    wells = sorted(sdata.tables["global-table-imaging"].obs["Well"].astype(str).unique())
    for well in wells:
        name = f"full_{well}.zarr"
        if (out / f"{name}.zip").exists():
            continue
        prep.subset(sdata, [well]).write(out / name, overwrite=True)
        # stored (-0), not deflated: zarr chunks are already compressed
        subprocess.run(["zip", "-0", "-r", "-q", f"{name}.zip", name], cwd=out, check=True)
        shutil.rmtree(out / name)
        print(f"{name}.zip  {(out / f'{name}.zip').stat().st_size / 1e9:.1f} GB", flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:])
