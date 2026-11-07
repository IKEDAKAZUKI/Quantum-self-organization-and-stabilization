"""Keep generated outputs separate from the distributed scientific inputs."""
from pathlib import Path
from .io import ROOT


def output_directory(path: str | Path) -> Path:
    path=Path(path).resolve()
    # Results under the package must be under results/. Else use any external folder.
    if path.is_relative_to(ROOT) and not path.is_relative_to(ROOT/'results'):
        raise ValueError('Choose results/ within the project or an external output directory.')
    path.mkdir(parents=True,exist_ok=True)
    return path
