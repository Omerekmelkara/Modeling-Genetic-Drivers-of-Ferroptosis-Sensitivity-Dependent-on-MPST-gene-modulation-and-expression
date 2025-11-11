# -*- coding: utf-8 -*-
"""
Created on Mon May 26 15:19:24 2025

@author: Ekmel
"""

import pandas as pd

# --- Step 1: Load mapping and create dictionary ---
mapping_path = r"C:/Users/Ekmel/Desktop/CBS expression corallation with Ferroptosis Inhibitors in Erastin Resistant Cells/CSV and XLSX files/depmap data/uniprot_hugo_entrez_id_mapping.csv"
mapping_df = pd.read_csv(mapping_path)
mapping_dict = dict(zip(mapping_df["UniprotID"], mapping_df["Symbol"]))

# --- Step 2: Load Gygi data and replace Uniprot IDs with Symbols ---
gygi_path = r"C:/Users/Ekmel/Desktop/CBS expression corallation with Ferroptosis Inhibitors in Erastin Resistant Cells/CSV and XLSX files/depmap data/harmonized_MS_CCLE_Gygi.csv"
gygi_df = pd.read_csv(gygi_path)

# Replace column headers using mapping dict; keep original if no mapping
new_columns = [mapping_dict.get(col, col) for col in gygi_df.columns]
gygi_df.columns = new_columns

# Save intermediate result with symbol headers
gygi_with_symbols_path = r"C:/Users/Ekmel/Desktop/CBS expression corallation with Ferroptosis Inhibitors in Erastin Resistant Cells/CSV and XLSX files/depmap data/gygi_with_symbols.csv"
gygi_df.to_csv(gygi_with_symbols_path, index=False)

# --- Step 3: Reload to ensure fresh state and select columns ---
gygi_df = pd.read_csv(gygi_with_symbols_path)

# Define columns to keep: cell ID plus target proteins
target_genes = ["Unnamed: 0", "GPX4","AIFM2","SLC7A11","CARS1", "CBS", "TST", "MPST", ]

# Check which columns exist before selecting (avoid KeyError)
available_cols = [col for col in target_genes if col in gygi_df.columns]
subset_df = gygi_df[available_cols]

# Save selected genes subset
selected_genes_path = r"C:/Users/Ekmel/Desktop/CBS expression corallation with Ferroptosis Inhibitors in Erastin Resistant Cells/CSV and XLSX files/depmap data/gygi_with_selected_genes.csv"
subset_df.to_csv(selected_genes_path, index=False)

# --- Step 4: Load Erastin resistance data ---
erastin_path = r"C:/Users/Ekmel/Desktop/CBS expression corallation with Ferroptosis Inhibitors in Erastin Resistant Cells/CSV and XLSX files/depmap data/ERASTIN (BRD_BRD-A25004090-001-08-4) log2 fold change.csv"
erastin_df = pd.read_csv(erastin_path)

# --- Step 5: Normalize and rename ID columns for merging ---
subset_df.rename(columns={"Unnamed: 0": "CellID"}, inplace=True)
erastin_df.rename(columns={"model": "CellID"}, inplace=True)

# Optional: Strip whitespace and uppercase to avoid mismatches
subset_df["CellID"] = subset_df["CellID"].astype(str).str.strip().str.upper()
erastin_df["CellID"] = erastin_df["CellID"].astype(str).str.strip().str.upper()

# --- Step 6: Merge datasets on CellID (intersection only) ---
merged_df = pd.merge(subset_df, erastin_df, on="CellID", how="inner")

# --- Step 7: Save combined data ---
combined_path = r"C:/Users/Ekmel/Desktop/CBS expression corallation with Ferroptosis Inhibitors in Erastin Resistant Cells/combined_protein_erastin_resistance.csv"
merged_df.to_csv(combined_path, index=False)
