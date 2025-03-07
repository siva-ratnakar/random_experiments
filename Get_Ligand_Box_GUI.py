import os
import tkinter as tk
from tkinter import filedialog, messagebox
import statistics

def is_ligand_atom(line):
    # Exclude standard amino acids & protein backbone atoms
    amino_acids = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
    if line.startswith("HETATM"):
        residue_name = line[17:20].strip()
        return residue_name not in amino_acids  # Keep only non-protein atoms
    return False

def get_ligand_coordinates(data, ext):
    coord = [[float(line[31:38]), float(line[39:46]), float(line[47:54])] 
             for line in data if is_ligand_atom(line)]
    
    if not coord:
        messagebox.showerror("Error", "No ligand detected in file!")
        return None
    
    xcoor, ycoor, zcoor = zip(*coord)
    return [list(xcoor), list(ycoor), list(zcoor)]

def min_max(coord):
    return [min(coord), max(coord)]

def center_XYZ(coord_range):
    return round(statistics.mean(coord_range), 3)

def length_WHD(coord_range, scale):
    return round(abs(coord_range[0] - coord_range[1]) * scale, 3)

def calculate_gridbox():
    file_path = file_entry.get()
    if not os.path.exists(file_path):
        messagebox.showerror("Error", "Invalid file path!")
        return
    
    try:
        with open(file_path, 'r') as f:
            data = f.readlines()
        ext = os.path.splitext(file_path)[-1]
        
        COOR = get_ligand_coordinates(data, ext)
        if COOR is None:
            return
        
        X, Y, Z = COOR
        ranges = [min_max(X), min_max(Y), min_max(Z)]
        center = [center_XYZ(r) for r in ranges]
        bxsize = [length_WHD(r, 1) for r in ranges]  # Scale = 1 for Angstrom units
        
        result_text.set(f"Grid Box Center: X={center[0]}, Y={center[1]}, Z={center[2]}\n"
                        f"Grid Box Size: X={bxsize[0]}, Y={bxsize[1]}, Z={bxsize[2]}")
    except Exception as e:
        messagebox.showerror("Error", f"Failed to process file: {e}")

def browse_file():
    file_path = filedialog.askopenfilename(filetypes=[("PDB files", "*.pdb"), ("PDBQT files", "*.pdbqt")])
    file_entry.delete(0, tk.END)
    file_entry.insert(0, file_path)

def quit_app():
    root.destroy()

# GUI Setup
root = tk.Tk()
root.title("Grid Box Calculator")
root.geometry("500x300")

tk.Label(root, text="Select Ligand PDB File:").pack(pady=5)
file_entry = tk.Entry(root, width=50)
file_entry.pack(pady=5)
browse_button = tk.Button(root, text="Browse", command=browse_file)
browse_button.pack(pady=5)

calc_button = tk.Button(root, text="Calculate Grid Box", command=calculate_gridbox)
calc_button.pack(pady=10)

result_text = tk.StringVar()
result_label = tk.Label(root, textvariable=result_text, justify="left")
result_label.pack(pady=10)

quit_button = tk.Button(root, text="Quit", command=quit_app)
quit_button.pack(pady=10)

root.mainloop()
