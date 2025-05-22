# PGRN Disulfide Bridges Analysis Tool

## Overview
This tool analyzes and visualizes disulfide bridges in protein structures, with a specific focus on Progranulin (PGRN). Disulfide bonds are crucial structural elements in proteins that contribute to their stability and functional properties. The tool offers both 2D schematic representations and interactive 3D visualizations of these important structural features.

## Features
- **Automated Detection**: Identifies disulfide bridges using precise geometric criteria (S-S distance and dihedral angles)
- **2D Visualization**: Creates detailed diagrams showing the positions and properties of disulfide bonds
- **3D Interactive Model**: Generates browser-based 3D visualizations with interactive controls (requires py3Dmol)
- **Comprehensive Analysis**: Provides distances and dihedral angles for each detected bond
- **Structure Modification**: Adds detected bonds back to the protein structure for further analysis

## Requirements
- Python 3.6+
- Required packages:
  - numpy
  - matplotlib
  - biotite
  - py3Dmol (optional, for 3D visualization)

## Installation
```bash
pip install numpy matplotlib biotite
pip install py3Dmol  # Optional, for 3D visualization
```

## Usage
Run the script directly:
```bash
python disulfide_bridges_analysis_pgrn.py
```

The script will:
1. Fetch the PGRN structure from the PDB database
2. Detect disulfide bonds using specified criteria
3. Generate visualizations
4. Save modified structure with disulfide bonds

## Output Files
- **pgrn_disulfide_bonds.svg**: 2D diagram of disulfide bridges (vector format)
- **pgrn_disulfide_bonds.png**: 2D diagram of disulfide bridges (raster format)
- **pgrn_disulfide_bonds_3d.html**: Interactive 3D visualization (requires web browser)
- **pgrn_with_disulfide_bonds.bcif**: Modified structure file with disulfide bonds included

## Visualization Details
### 2D Visualization
- Displays protein sequence with highlighted cysteine residues
- Shows disulfide bridges as arcs connecting cysteine pairs
- Color-coded by bond distance
- Includes labels with bond metrics (distance and dihedral angle)

### 3D Visualization
- Interactive model viewable in any modern web browser
- Color-coded backbone with highlighted cysteine residues
- Disulfide bonds shown as orange cylinders
- Interactive controls for rotation, zoom, and feature toggling

## Customization
You can modify the detection criteria by changing parameters in the `detect_disulfide_bonds` function:
- `distance`: Target S-S distance (default: 2.05 Å)
- `distance_tol`: Distance tolerance (default: 0.05 Å)
- `dihedral`: Target dihedral angle (default: 90°)
- `dihedral_tol`: Dihedral angle tolerance (default: 10°)

## Scientific Background
Progranulin (PGRN) is a glycoprotein with important roles in various biological processes, including inflammation, wound healing, and neurodegeneration. The disulfide bonds in PGRN contribute significantly to its structural stability and functional domains. Recent research has highlighted PGRN's protective role in cardiovascular health through various mechanisms.

## Citation
If you use this tool in your research, please cite:

Qiao, G., Lu, Y., Wu, J., Ren, C., Lin, R., & Zhang, C. (2025). Progranulin's Protective Mechanisms: A New Frontier in Cardiovascular Disease. *Cells*, Special Issue "Molecular Pathogenesis of Cardiovascular Diseases".

```bibtex
@article{qiao2025progranulin,
  title={Progranulin's Protective Mechanisms: A New Frontier in Cardiovascular Disease},
  author={Qiao, Gan and Lu, Yongxiang and Wu, Jianping and Ren, Chunyang and Lin, Roudian and Zhang, Chunxiang},
  journal={Cells},
  year={2025},
  note={Special Issue "Molecular Pathogenesis of Cardiovascular Diseases"}
}
```
## Contact
For questions or support, please contact:
- Gan Qiao (dqz377977905@swmu.edu.cn)
- Chunxiang Zhang (zhangchunxiang@swmu.edu.cn)
