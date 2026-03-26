import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("SUNPY_CONFIGDIR", str(Path("/tmp/solpolpy-sunpy-config")))

Path(os.environ["SUNPY_CONFIGDIR"]).mkdir(parents=True, exist_ok=True)
