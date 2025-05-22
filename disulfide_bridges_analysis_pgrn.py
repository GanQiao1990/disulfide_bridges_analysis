#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Detection of disulfide bonds in PGRN

This script detects disulfide bridges in protein structures, especially for PGRN.
The detected disulfide bonds are visualized and added to the 'bonds' attribute
of the AtomArray.

The employed criteria for disulfide bonds are quite simple:
the S_γ atoms of two cystein residues must be in a vicinity
of 2.05 ± 0.05 Å and the dihedral angle of
C_β - S_γ - S′_γ - C′_β must be 90 ± 10°.
"""

# Code source: Patrick Kunzmann (original biotite example)
# Modified for PGRN analysis

import io
import os
from tempfile import gettempdir
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
import biotite.sequence as seq
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import biotite.database.rcsb as rcsb

# Optional imports for 3D visualization
try:
    import py3Dmol
    PY3DMOL_AVAILABLE = True
except ImportError:
    PY3DMOL_AVAILABLE = False


def detect_disulfide_bonds(structure, distance=2.05, distance_tol=0.05,
                          dihedral=90, dihedral_tol=10):
    """
    Detect disulfide bonds in a protein structure.
    
    Parameters
    ----------
    structure : AtomArray
        The protein structure to analyze.
    distance : float, optional
        The optimal S-S distance in Å. Default is 2.05.
    distance_tol : float, optional
        The tolerance for the S-S distance in Å. Default is 0.05.
    dihedral : float, optional
        The optimal dihedral angle in degrees. Default is 90.
    dihedral_tol : float, optional
        The tolerance for the dihedral angle in degrees. Default is 10.
    
    Returns
    -------
    ndarray
        Array of index pairs, where each pair represents a disulfide bond.
    """
    # Array where detected disulfide bonds are stored
    disulfide_bonds = []
    # A mask that selects only S-gamma atoms of cysteins
    sulfide_mask = (structure.res_name == "CYS") & \
                   (structure.atom_name == "SG")
    # sulfides in adjacency to other sulfides are detected in an
    # efficient manner via a cell list
    cell_list = struc.CellList(
        structure,
        cell_size=distance+distance_tol,
        selection=sulfide_mask
    )
    # Iterate over every index corresponding to an S-gamma atom
    for sulfide_i in np.where(sulfide_mask)[0]:
        # Find indices corresponding to other S-gamma atoms,
        # that are adjacent to the position of structure[sulfide_i]
        # We use the faster 'get_atoms_in_cells()' instead of
        # `get_atoms()`, as precise distance measurement is done
        # afterwards anyway
        potential_bond_partner_indices = cell_list.get_atoms_in_cells(
            coord=structure.coord[sulfide_i]
        )
        # Iterate over every index corresponding to an S-gamma atom
        # as bond partner
        for sulfide_j in potential_bond_partner_indices:
            if sulfide_i == sulfide_j:
                # A sulfide cannot create a bond with itself:
                continue
            # Create 'Atom' instances
            # of the potentially bonds S-gamma atoms
            sg1 = structure[sulfide_i]
            sg2 = structure[sulfide_j]
            # For dihedral angle measurement the corresponding
            # C-beta atoms are required, too
            cb1 = structure[
                (structure.chain_id == sg1.chain_id) &
                (structure.res_id == sg1.res_id) &
                (structure.atom_name == "CB")
            ]
            cb2 = structure[
                (structure.chain_id == sg2.chain_id) &
                (structure.res_id == sg2.res_id) &
                (structure.atom_name == "CB")
            ]
            # Measure distance and dihedral angle and check criteria
            bond_dist = struc.distance(sg1, sg2)
            bond_dihed = np.abs(np.rad2deg(struc.dihedral(cb1, sg1, sg2, cb2)))
            if bond_dist  > distance - distance_tol and \
               bond_dist  < distance + distance_tol and \
               bond_dihed > dihedral - dihedral_tol and \
               bond_dihed < dihedral + dihedral_tol:
                    # Atom meet criteria -> we found a disulfide bond
                    # -> the indices of the bond S-gamma atoms
                    # are put into a tuple with the lower index first
                    bond_tuple = sorted((sulfide_i, sulfide_j))
                    # Add bond to list of bonds, but each bond only once
                    if bond_tuple not in disulfide_bonds:
                        disulfide_bonds.append(bond_tuple)
    return np.array(disulfide_bonds, dtype=int)


def plot_disulfide_bonds(structure, disulfide_bonds, rows=6, output_file='disulfide_bonds_plot.svg'):
    """
    Create a visualization of disulfide bonds.
    
    Parameters
    ----------
    structure : AtomArray
        The protein structure.
    disulfide_bonds : ndarray
        Array of disulfide bond index pairs.
    rows : int, optional
        Number of rows for the plot. Default is 6.
    output_file : str, optional
        Output file name. Default is 'disulfide_bonds_plot.svg'.
    """
    # Create a sequence object for each residue in the structure
    sequence = seq.ProteinSequence(structure.res_name[structure.atom_name == "CA"])

    # Calculate the number of residues per row
    total_residues = len(sequence)
    residues_per_row = (total_residues + rows - 1) // rows  # Ceiling division to ensure all residues are shown

    # Create the figure with specified rows - increase the height for more space
    fig, axs = plt.subplots(rows, 1, figsize=(20, rows * 3))
    fig.suptitle("Disulfide Bridges Visualization", fontsize=16, fontweight='bold')
    plt.subplots_adjust(hspace=0.6, top=0.95)  # Increased spacing between rows

    # Create a custom colormap for the disulfide bonds based on distance
    colors = [(0.0, 'blue'), (0.5, 'green'), (1.0, 'red')]
    cmap = LinearSegmentedColormap.from_list('disulfide_cmap', colors, N=100)
    
    # Group disulfide bonds by residue to identify potential overlaps
    disulfide_by_residue = {}
    for sg1_index, sg2_index in disulfide_bonds:
        sg1_res_id = structure.res_id[sg1_index]
        sg2_res_id = structure.res_id[sg2_index]
        if sg1_res_id not in disulfide_by_residue:
            disulfide_by_residue[sg1_res_id] = []
        if sg2_res_id not in disulfide_by_residue:
            disulfide_by_residue[sg2_res_id] = []
        disulfide_by_residue[sg1_res_id].append((sg1_index, sg2_index))
        disulfide_by_residue[sg2_res_id].append((sg1_index, sg2_index))
    
    # Calculate distances and angles for all disulfide bonds
    bond_distances = []
    bond_dihedrals = []
    for sg1_index, sg2_index in disulfide_bonds:
        sg1 = structure[sg1_index]
        sg2 = structure[sg2_index]
        cb1 = structure[
            (structure.chain_id == sg1.chain_id) &
            (structure.res_id == sg1.res_id) &
            (structure.atom_name == "CB")
        ][0]  # Get the first (and should be only) atom
        cb2 = structure[
            (structure.chain_id == sg2.chain_id) &
            (structure.res_id == sg2.res_id) &
            (structure.atom_name == "CB")
        ][0]  # Get the first (and should be only) atom
        bond_dist = struc.distance(sg1, sg2)
        bond_dihed = np.abs(np.rad2deg(struc.dihedral(cb1, sg1, sg2, cb2)))
        bond_distances.append(bond_dist)
        bond_dihedrals.append(bond_dihed)
    
    # Normalize distances for coloring
    if bond_distances:
        min_dist = min(bond_distances)
        max_dist = max(bond_distances)
        norm_distances = [(d - min_dist) / (max_dist - min_dist) if max_dist > min_dist else 0.5 for d in bond_distances]
    else:
        norm_distances = []
    
    MARGIN = 0.2
    # Track bonds that have been annotated to avoid duplicates
    annotated_bonds = set()

    for row in range(rows):
        ax = axs[row] if rows > 1 else axs
        start_residue = row * residues_per_row + 1
        end_residue = min((row + 1) * residues_per_row, total_residues)
        
        # Set more vertical space for bridges and annotations
        ax.set_xlim(start_residue - MARGIN, end_residue + MARGIN)
        ax.set_ylim(-0.2, 2.0)  # Increased vertical space
        ax.set_xticks(np.arange(start_residue, end_residue + 1))
        
        # Highlight cysteine residues
        cys_positions = [i+1 for i, aa in enumerate(sequence) if aa == 'C']
        sequence_labels = list(str(sequence)[start_residue-1:end_residue])
        
        # Use the original sequence without adding residue numbers to cysteines
        ax.set_xticklabels(sequence_labels, fontdict={'weight': 'bold'}, rotation=0)
        
        # Highlight cysteine residues
        for pos in cys_positions:
            if start_residue <= pos <= end_residue:
                ax.axvspan(pos-0.5, pos+0.5, alpha=0.2, color='yellow')
        
        ax.set_title(f"Residues {start_residue}-{end_residue}", fontsize=10)
        ax.yaxis.set_tick_params(left=False, right=False, labelleft=False, labelright=False)
        ax.xaxis.set_tick_params(bottom=True, top=False, labelbottom=True, labeltop=False, width=0)
        ax.set_frame_on(True)
        ax.grid(axis='x', linestyle='--', alpha=0.3)
        
        # Draw the protein backbone as a thick line
        ax.axhline(y=0, xmin=0, xmax=1, color='black', linewidth=2)
        
        # Create a dictionary to store label positions for this row to avoid overlaps
        label_positions = {}
        
        # Step 1: First pass - add all bridges
        for i, (sg1_index, sg2_index) in enumerate(disulfide_bonds):
            sg1_res_id = structure.res_id[sg1_index]
            sg2_res_id = structure.res_id[sg2_index]
            sg1_chain = structure.chain_id[sg1_index]
            sg2_chain = structure.chain_id[sg2_index]
            
            # Get bond color based on distance
            bond_color = cmap(norm_distances[i] if i < len(norm_distances) else 0.5)
            
            # Calculate vertical offset for potentially overlapping bonds
            # More bonds on a residue = higher arc
            if sg1_res_id in disulfide_by_residue and sg2_res_id in disulfide_by_residue:
                count1 = len(disulfide_by_residue[sg1_res_id])
                count2 = len(disulfide_by_residue[sg2_res_id])
                vertical_offset = 0.8 + 0.3 * max(count1, count2) / 3  # Scale the height based on bond count
            else:
                vertical_offset = 0.8
                
            if start_residue <= sg1_res_id <= end_residue or start_residue <= sg2_res_id <= end_residue:
                ellipse_center = (sg1_res_id + sg2_res_id) / 2
                ellipse_width = abs(sg2_res_id - sg1_res_id)
                
                # Create a bond key for annotation tracking
                bond_key = tuple(sorted([sg1_res_id, sg2_res_id]))
                
                if start_residue <= ellipse_center <= end_residue and start_residue <= sg1_res_id <= end_residue and start_residue <= sg2_res_id <= end_residue:
                    # Full bridge within the row - draw as a semi-transparent bridge-like curve
                    # The further apart the cysteines, the higher the arc
                    height_factor = min(1.5, 0.3 + ellipse_width / 50)  # Cap the height for very long spans
                    
                    # Create points for the bridge
                    x = np.linspace(sg1_res_id, sg2_res_id, 100)
                    # Parabolic arc
                    y = height_factor * vertical_offset * 4 * (x - sg1_res_id) * (x - sg2_res_id) / -(ellipse_width**2)
                    
                    # Draw the bridge with gradient color and higher line width
                    ax.plot(x, y, color=bond_color, linewidth=3, solid_capstyle='round', zorder=5)
                    
                    # Add markers at the cysteine positions
                    ax.scatter([sg1_res_id, sg2_res_id], [0, 0], color=bond_color, s=50, zorder=10)
                    
                    # Calculate the highest point of the bridge for text placement
                    max_y = height_factor * vertical_offset
                    
                    # Add bond information as text annotation (only once per bond)
                    if bond_key not in annotated_bonds:
                        dist = bond_distances[i] if i < len(bond_distances) else "?"
                        dihed = bond_dihedrals[i] if i < len(bond_dihedrals) else "?"
                        
                        # Create a more compact label format
                        label = f"SS-bond: C{sg1_res_id}{sg1_chain}-C{sg2_res_id}{sg2_chain}\n{dist:.2f}Å, {dihed:.1f}°"
                        
                        # Position the label to avoid overlap with other labels
                        # Adjust the y-position to be above the highest point of the bridge
                        if ellipse_center in label_positions:
                            # If we already have a label at this x-position, stack it higher
                            y_pos = label_positions[ellipse_center] + 0.3  # Increase y-position
                        else:
                            y_pos = max_y + 0.3  # Default position above the bridge
                        
                        label_positions[ellipse_center] = y_pos  # Save position for future overlap checks
                        
                        ax.annotate(label, xy=(ellipse_center, max_y), xytext=(0, 15), 
                                   textcoords="offset points", ha='center', va='bottom',
                                   bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.9, ec=bond_color),
                                   arrowprops=dict(arrowstyle='->', color=bond_color, connectionstyle='arc3,rad=0.1'))
                        annotated_bonds.add(bond_key)
                else:
                    # Partial connections for bonds spanning multiple rows
                    if sg1_res_id >= start_residue and sg1_res_id <= end_residue:
                        # Start of the bond - draw as an up-arrow
                        ax.annotate("", xy=(sg1_res_id, 0), xytext=(sg1_res_id, 0.8), 
                                    arrowprops=dict(arrowstyle="->", color=bond_color, linewidth=2, mutation_scale=15))
                        
                        # Add a label indicating the other end
                        label_text = f"→C{sg2_res_id}{sg2_chain}"
                        ax.annotate(label_text, xy=(sg1_res_id, 0.85), xytext=(5, 0),
                                     textcoords="offset points", ha='left', va='center',
                                     bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.7, ec=bond_color))
                        
                    elif sg2_res_id >= start_residue and sg2_res_id <= end_residue:
                        # End of the bond - draw as an up-arrow
                        ax.annotate("", xy=(sg2_res_id, 0), xytext=(sg2_res_id, 0.8), 
                                    arrowprops=dict(arrowstyle="->", color=bond_color, linewidth=2, mutation_scale=15))
                        
                        # Add a label indicating the other end
                        label_text = f"←C{sg1_res_id}{sg1_chain}"
                        ax.annotate(label_text, xy=(sg2_res_id, 0.85), xytext=(-5, 0),
                                     textcoords="offset points", ha='right', va='center',
                                     bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.7, ec=bond_color))

    # Add a colorbar to show the distance scale
    sm = cm.ScalarMappable(cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axs, orientation='horizontal', pad=0.05, aspect=40, shrink=0.8)
    cbar.set_label('Normalized S-S Distance (Blue=Shortest, Red=Longest)')

    # Add title and legend
    plt.figtext(0.5, 0.01, "Note: Bridges show disulfide bonds between cysteine residues", 
                ha="center", fontsize=12, bbox={"facecolor":"white", "alpha":0.8, "pad":5})

    # Save the figure
    plt.savefig(output_file, format='svg', dpi=300, bbox_inches='tight')
    print(f"The plot has been saved as '{output_file}'")
    
    # Also save as PNG for easier viewing
    png_output = os.path.splitext(output_file)[0] + '.png'
    plt.savefig(png_output, format='png', dpi=300, bbox_inches='tight')
    print(f"The plot has also been saved as '{png_output}'")
    
    plt.close(fig)  # Close the figure to free up memory


def visualize_3d(structure, disulfide_bonds, output_file='disulfide_3d.html'):
    """
    Create an interactive 3D visualization of the protein structure and disulfide bonds
    using py3Dmol.
    
    Parameters
    ----------
    structure : AtomArray
        The protein structure.
    disulfide_bonds : ndarray
        Array of index pairs, where each pair represents a disulfide bond.
    output_file : str, optional
        Output HTML file name for the 3D visualization.
    """
    if not PY3DMOL_AVAILABLE:
        print("Warning: py3Dmol is not installed. Skipping 3D visualization.")
        print("To install py3Dmol, run: pip install py3Dmol")
        return
    
    try:
        # First try to use a PDB format which is more compatible with py3Dmol
        # Create a PDB file since it's more reliably parsed by py3Dmol
        import biotite.structure.io.pdb as pdb
        temp_file = os.path.join(gettempdir(), "temp_structure.pdb")
        pdb_file = pdb.PDBFile()
        pdb.set_structure(pdb_file, structure)
        with open(temp_file, "w") as f:
            pdb_file.write(f)
        
        # Create a py3Dmol view with larger dimensions
        view = py3Dmol.view(width=1500, height=1500)
        
        with open(temp_file, 'r') as f:
            view.addModel(f.read(), 'pdb')
    except Exception as e:
        print(f"Trying CIF format instead due to error: {e}")
        # Fallback to CIF if PDB writing fails
        temp_file = os.path.join(gettempdir(), "temp_structure.cif")
        temp_pdbx = pdbx.BinaryCIFFile()
        pdbx.set_structure(temp_pdbx, structure)
        with open(temp_file, "wb") as f:
            temp_pdbx.write(f)
        
        # Create a py3Dmol view with larger dimensions
        view = py3Dmol.view(width=1500, height=1500)
        
        try:
            with open(temp_file, 'r') as f:
                view.addModel(f.read(), 'cif')
        except UnicodeDecodeError:
            # If reading as text fails, try binary mode
            with open(temp_file, 'rb') as f:
                view.addModel(f.read().decode('latin1'), 'cif')
    
    # Set up a simplified visualization style focusing only on main backbone and disulfide bridges
    # Set all atoms to be invisible by default
    view.setStyle({}, {'cartoon': {'color': 'gray', 'opacity': 0.9}})
    
    # Special highlighting just for cysteines to make them stand out
    view.setStyle({'resn': 'CYS'}, {'cartoon': {'color': 'yellow', 'opacity': 0.8}})
    
    # Show only sulfur atoms of cysteines (for disulfide bridges)
    view.addStyle({'resn': 'CYS', 'atom': 'SG'}, {'sphere': {'radius': 0.6, 'color': 'yellow', 'opacity': 1.0}})
    
    # Hide all side chains and other non-backbone atoms
    view.removeAllSurfaces()
    view.removeAllLabels()
    
    # Add disulfide bonds visualization
    for sg1_index, sg2_index in disulfide_bonds:
        sg1 = structure[sg1_index]
        sg2 = structure[sg2_index]
        sg1_res_id = sg1.res_id
        sg2_res_id = sg2.res_id
        sg1_chain = sg1.chain_id
        sg2_chain = sg2.chain_id
        
        # Convert numpy values to native Python floats to prevent JSON serialization issues
        sg1_x, sg1_y, sg1_z = float(sg1.coord[0]), float(sg1.coord[1]), float(sg1.coord[2])
        sg2_x, sg2_y, sg2_z = float(sg2.coord[0]), float(sg2.coord[1]), float(sg2.coord[2])
        
        # Add a cylinder to represent the disulfide bond
        view.addCylinder({
            'start': {'x': sg1_x, 'y': sg1_y, 'z': sg1_z},
            'end': {'x': sg2_x, 'y': sg2_y, 'z': sg2_z},
            'radius': 0.2,
            'color': 'orange',
            'fromCap': 1,
            'toCap': 1
        })
        
        # Add labels for the cysteines, place them slightly offset for better visibility
        view.addLabel(f"C{sg1_res_id}{sg1_chain}-C{sg2_res_id}{sg2_chain}", 
                       {'position': {'x': (sg1_x+sg2_x)/2, 'y': (sg1_y+sg2_y)/2 + 0.5, 'z': (sg1_z+sg2_z)/2}, 
                       'backgroundColor': 'orange', 'fontColor': 'black', 'backgroundOpacity': 0.7, 'fontSize': 15})
    
    # Configure view settings
    view.zoomTo()
    view.setClickable({'clickable': True})
    view.setColorByElement({'carbon': 'gray'})
    view.rotate(30, 'y', animate=True)  # Rotates to see structure better
    
    # Generate comprehensive HTML with embedded styles and controls
    html_content = """<!DOCTYPE html>
<html>
<head>
    <title>PGRN Disulfide Bridges 3D Visualization</title>
    <script src="https://unpkg.com/jquery"></script>
    <script src="https://unpkg.com/3dmol/build/3Dmol-min.js"></script>
    <style>
        body { font-family: Arial, sans-serif; margin: 0; padding: 0; }
        #container { display: flex; flex-direction: column; height: 100vh; }
        #viewer { flex-grow: 1; width: 100%; height: 85vh; position: relative; }
        #controls { background-color: #f0f0f0; padding: 10px; }
        button { margin: 5px; padding: 5px 10px; background-color: #4CAF50; color: white; border: none; border-radius: 4px; cursor: pointer; }
        button:hover { background-color: #45a049; }
        #infoText { margin: 10px; padding: 10px; background-color: #e8f5e9; border-radius: 5px; }
    </style>
</head>
<body>
    <div id="container">
        <div id="controls">
            <button onclick="toggleSpin()">Toggle Spin</button>
            <button onclick="resetView()">Reset View</button>
            <button onclick="toggleLabels()">Toggle Labels</button>
            <button onclick="toggleBackbone()">Toggle Backbone</button>
            <button onclick="toggleBonds()">Toggle Disulfide Bonds</button>
        </div>
        <div id="infoText">
            <strong>PGRN Disulfide Bridges 3D Visualization</strong>: 
            Click and drag to rotate. Scroll to zoom. Click on atoms for info. Use buttons to control the view.
        </div>
        <div id="viewer"></div>
    </div>
    <script>
        let spinning = false;
        let viewer;
        let labelsVisible = true;
        let bondsVisible = true;
        let backboneVisible = true;
        
        $(document).ready(function() {
            viewer = $3Dmol.createViewer('viewer');
            
            // Load the data and setup from the embedded JSON
            let viewerData = VIEWERDATA;
            
            viewer.setBackgroundColor('white');
            viewer.fromJSON(viewerData);
            viewer.render();
            
            toggleSpin(); // Start with spinning on
        });
        
        function toggleSpin() {
            spinning = !spinning;
            viewer.spin(spinning);
        }
        
        function resetView() {
            viewer.zoomTo();
            viewer.render();
        }
        
        function toggleLabels() {
            labelsVisible = !labelsVisible;
            viewer.getLabels().forEach(function(label) {
                label.setVisibility(labelsVisible);
            });
            viewer.render();
        }
        
        function toggleBackbone() {
            backboneVisible = !backboneVisible;
            if (backboneVisible) {
                viewer.setStyle({}, {cartoon: {color: 'gray', opacity: 0.8}});
                viewer.setStyle({resn: 'CYS'}, {cartoon: {color: 'yellow'}});
            } else {
                viewer.setStyle({}, {cartoon: {opacity: 0}});
            }
            viewer.render();
        }
        
        function toggleBonds() {
            bondsVisible = !bondsVisible;
            // The cylinders are used for the disulfide bonds
            viewer.getCylinders().forEach(function(cyl) {
                cyl.setVisibility(bondsVisible);
            });
            viewer.render();
        }
    </script>
</body>
</html>
    """
    
    # Use the simpler write_html method that works correctly
    html_with_data = view.write_html()
    
    # Save the HTML file
    with open(output_file, 'w') as f:
        f.write(html_with_data)
    
    print(f"Interactive 3D visualization saved as '{output_file}'")
    print(f"Open the file directly in your browser to view the 3D visualization.")
    
    # Create a basic viewer opener script (Python) for convenience
    opener_script = f"""import webbrowser
import os

html_file = os.path.abspath("{output_file}")
webbrowser.open('file://' + html_file)
"""
    
    opener_file = os.path.splitext(output_file)[0] + '_open.py'
    with open(opener_file, 'w') as f:
        f.write(opener_script)
    
    print(f"Created {opener_file} - run this Python script to automatically open the visualization in a browser.")
    

def main():
    """
    Main function that loads the structure, detects disulfide bonds,
    visualizes them, and adds them back to the structure.
    """
    # Fetch the structure (using an AlphaFold model of PGRN)
    print("Fetching structure...")
    pdbx_file = pdbx.BinaryCIFFile.read(
        rcsb.fetch("AF_AFP28797F1", "bcif", gettempdir())
    )
    structure = pdbx.get_structure(pdbx_file, include_bonds=True, model=1)
    
    # Find and print existing disulfide bonds
    sulfide_indices = np.where(
        (structure.res_name == "CYS") & (structure.atom_name == "SG")
    )[0]
    
    # Remove existing disulfide bonds for testing purposes
    print("Original disulfide bonds:")
    for i, j, _ in structure.bonds.as_array():
        if i in sulfide_indices and j in sulfide_indices:
            print(structure[i])
            print(structure[j])
            print()
            structure.bonds.remove_bond(i, j)
    
    # Detect disulfide bonds
    print("Detecting disulfide bonds...")
    disulfide_bonds = detect_disulfide_bonds(structure)
    
    # Print detected disulfide bonds with distances and dihedral angles
    print("\nDetected disulfide bonds:")
    print(f"{'Cys 1':<15} {'Cys 2':<15} {'Distance (Å)':<15} {'Dihedral (°)':<15}")
    print("-" * 60)
    for sg1_index, sg2_index in disulfide_bonds:
        sg1 = structure[sg1_index]
        sg2 = structure[sg2_index]
        cb1 = structure[
            (structure.chain_id == sg1.chain_id) &
            (structure.res_id == sg1.res_id) &
            (structure.atom_name == "CB")
        ][0]  # Get the first (and should be only) atom
        cb2 = structure[
            (structure.chain_id == sg2.chain_id) &
            (structure.res_id == sg2.res_id) &
            (structure.atom_name == "CB")
        ][0]  # Get the first (and should be only) atom
        
        # Calculate distance and dihedral angle
        bond_dist = struc.distance(sg1, sg2)
        bond_dihed = np.abs(np.rad2deg(struc.dihedral(cb1, sg1, sg2, cb2)))
        
        cys1 = f"C{sg1.res_id}{sg1.chain_id}"
        cys2 = f"C{sg2.res_id}{sg2.chain_id}"
        print(f"{cys1:<15} {cys2:<15} {bond_dist:<15.2f} {bond_dihed:<15.1f}")
    
    # Create 2D visualization
    print("\nCreating 2D visualization...")
    plot_disulfide_bonds(structure, disulfide_bonds, rows=6, output_file='pgrn_disulfide_bonds.svg')
    
    # Create 3D visualization if py3Dmol is available
    print("\nCreating 3D visualization...")
    visualize_3d(structure, disulfide_bonds, output_file='pgrn_disulfide_bonds_3d.html')
    
    # Add detected bonds back to the structure
    for sg1_index, sg2_index in disulfide_bonds:
        structure.bonds.add_bond(sg1_index, sg2_index, struc.BondType.SINGLE)
    
    # Write the structure back to a file
    out_file = pdbx.BinaryCIFFile()
    pdbx.set_structure(out_file, structure)
    with open("pgrn_with_disulfide_bonds.bcif", "wb") as f:
        out_file.write(f)
    
    print("\nStructure with disulfide bonds saved as 'pgrn_with_disulfide_bonds.bcif'")
    print("\nVisualization files created:")
    print("  - pgrn_disulfide_bonds.svg (2D diagram)")
    print("  - pgrn_disulfide_bonds.png (2D diagram)")
    if PY3DMOL_AVAILABLE:
        print("  - pgrn_disulfide_bonds_3d.html (Interactive 3D visualization)")
    else:
        print("  - 3D visualization skipped (py3Dmol not installed)")


if __name__ == "__main__":
    main()
