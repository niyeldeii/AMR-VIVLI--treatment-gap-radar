"""Central path config. Raw data lives ONE LEVEL UP from this repo and is never committed."""
from pathlib import Path

PKG_DIR = Path(__file__).resolve().parent.parent          # treatment_gap_radar/
DATA_ROOT = PKG_DIR.parent                                # challenge folder with the datasets
CONFIG_DIR = PKG_DIR / "config"
PROCESSED_DIR = PKG_DIR / "data_processed"
FIGURES_DIR = PKG_DIR / "figures"
PROCESSED_DIR.mkdir(exist_ok=True)
FIGURES_DIR.mkdir(exist_ok=True)

# Relative paths (from DATA_ROOT) to each raw dataset file.
RAW = {
    "ATLAS":            "ATLAS_Antibiotics/atlas_vivli_2004_2024.csv",
    "SOAR_201818":      "SOAR 201818/gsk_201818_published.csv",
    "SOAR_201910":      "SOAR 201910/GSK_SOAR_201910 raw data.xlsx",
    "SOAR_207965":      "SOAR 207965/SOAR 207965 Complete data set 04Sep25.xlsx",
    "INNOVIVA_ACINETO": "Surveillance of global clinical isolates of Acinetobacter baumannii-calcoaceticus complex collected from 2016-2021/IST-Entasis_Acinetobacter-Surveillance_2016-2021.xlsx",
    "DREAM_TB":         "Bedaquiline Drug Resistance Assessment in MDR-TB (DREAM)/BEDAQUILINE DREAM DATASET FOR VIVLI - 06-06-2022.xlsx",
    "KEYSTONE":         "KEYSTONE/Omadacycline_2015 to 2025_Surveillance_data.xlsx",
    "SIDERO_WT":        "SIDERO-WT/Updated_Shionogi Five year SIDERO-WT Surveillance data(without strain number)_Vivli_220409.xlsx",
    "GEARS":            "GEARS/Venatorx surveillance data_2024_06_06.xlsx",
    "PLEA_I":           "PLEA (Study I)/PLEA Study I (n=3150)_updated.xlsx",
    "PLEA_II":          "PLEA (Study II)/Study II (n=232).xlsx",
    "GASAR_III":        "GASAR (Study III)/GASAR Study III (n=494)_updated.xlsx",
    "SPIDAAR_ISOLATE":  "SPIDAAR RWE Study/spidaar_isolatedata.xls",
    "SPIDAAR_PATIENT":  "SPIDAAR RWE Study/spidaar_patientdata.xls",
    "RND_HUB":          "Projects.xlsx",
}

def raw_path(key: str) -> Path:
    return DATA_ROOT / RAW[key]
