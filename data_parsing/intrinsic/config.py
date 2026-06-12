"""Configuration for intrinsic-properties master build."""

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

# External physiology data (outside repo)
PHYSIOLOGY_DATA_ROOT = Path("/Users/rancze/Documents/Data/expresso_data/physiology")
DEFAULT_SOURCE = (
    PHYSIOLOGY_DATA_ROOT
    / "20260609_download_from_Drive"
    / "Intrinsic_properties_Analysis_Aug2024.xlsx"
)
DEFAULT_OUTPUT_DIR = PHYSIOLOGY_DATA_ROOT / "restructured"

PARAMETER_DATA_SHEETS = [
    "V2M_L5",
    "V1_L5",
    "V1_L2-3",
    "V2M_L2-3",
    "All Analysed data",
    "V2M_L5_2A_agonist",
    "V2M_L5_1A_Antagonist",
    "V2M_L5_MDL",
    "V2M_L5_TASK_Acidic_pH",
    "V2M_L5_Caesum",
]

DRUG_SHEETS = {
    "V2M_L5_2A_agonist": "5-HT2A agonist",
    "V2M_L5_1A_Antagonist": "5-HT1A antagonist",
    "V2M_L5_MDL": "MDL (5-HT2A antagonist)",
    "V2M_L5_TASK_Acidic_pH": "acidic pH (6.1)",
    "V2M_L5_Caesum": "cesium (intracellular)",
}

STANDARD_SHEETS = [
    "V2M_L5",
    "V1_L5",
    "V1_L2-3",
    "V2M_L2-3",
    "All Analysed data",
]

DEDUP_PRIORITY = [
    "All Analysed data",
    "V2M_L5",
    "V1_L5",
    "V1_L2-3",
    "V2M_L2-3",
    "V2M_L5_2A_agonist",
    "V2M_L5_1A_Antagonist",
    "V2M_L5_MDL",
    "V2M_L5_TASK_Acidic_pH",
    "V2M_L5_Caesum",
]

METADATA_SHEETS = {
    "Assumed_tlx": {"assumed_type": "Tlx"},
    "Assumed_PT_V2M": {"assumed_type": "ET"},
    "SC projecting cells": {"projection_target": "SC"},
}

ALL_CELLS_SHEET = "All cells"
CLUSTER_SHEET = "cluster_analysis_res"

# CCF code → broad area mapping
CCF_TO_AREA = {
    "VISp": "VISp",
    "VISam": "V2M",
    "VISpm": "V2M",
    "RSPagl": "V2M",
}

VALID_AREA_CCF = frozenset(CCF_TO_AREA)

# Known label merges (canonical ← variant)
EXPLICIT_LABEL_MERGES = {
    ("_short_depol", "First AP peak in RecordA )mV)1"): "First AP peak in RecordA (mV)",
    ("_short_depol", "Last AP peak in RecordA (mV)1"): "Last AP peak in RecordA (mV)",
    ("_chirp", "Phase lead integral"): "Phase lead integral (rad*Hz)",
    ("_chirp", "synchronous frequency (Hz)"): "Synchronous Frequency (Hz)",
}

SKIP_PARAM_LABELS = frozenset({"pH 6.1", "pH 7.3", "condition"})

MORPHOLOGY_ROW_LABELS = ("Morphology confirmed", "Morpholgy?")

METADATA_ROW_LABELS = {
    "Classic burster? 0/1": "classic_burster",
}

LABEL_MERGE_RATIO = 0.92
