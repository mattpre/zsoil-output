# zsoil-output

Python SDK for reading ZSoil FEM result files: mesh (`C_Mesh`), solver history
(`C_HistoryOfExecution`, `C_Rcf_info`), and per-element-type results readers
(`C_NodalResults`, `C_BeamResults`, `C_TrussResults`, `C_ShellResults`,
`C_Contact_2D_Results`, `C_Contact_3D_Results`, `C_ContinuumResults`, ...).

Formerly shipped only as a folder inside the ZSoil installation
(`Tools v<year>/ZSoilPy3`) and added to `sys.path` by hand; now also
installable as a regular Python package.

## Install

```bash
pip install git+https://github.com/ZSoil-FEM/zsoil-output.git
```

## Usage

Every module installs at the top level (no `zsoil_output.` prefix), matching
how the SDK has always been imported:

```python
from C_Mesh import Mesh
from C_HistoryOfExecution import HistoryOfExecution
from C_NodalResults import NodalResults
```

See [zsoil-output-examples](https://github.com/ZSoil-FEM/zsoil-output-examples)
for worked examples reading real ZSoil projects with this SDK.
