# imports
import os
import csv
import sys

import pandas as pd
import sascorer
from RAscore import RAscore_XGB
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams

RDLogger.DisableLog("rdApp.*")

input_file = sys.argv[1]
output_file = sys.argv[2]

root = os.path.dirname(os.path.abspath(__file__))

with open(input_file, "r") as f:
    reader = csv.reader(f)
    next(reader)  # skip header
    smiles_list = [r[0] for r in reader]


def get_mol_from_smiles(smiles):
    return Chem.MolFromSmiles(smiles)

# Initialize PAINS filters (PAINS A, B, and C)
pains_params = FilterCatalogParams()
pains_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)  # just use .PAINS
pains_catalog = FilterCatalog.FilterCatalog(pains_params)

# Initialize Brenk alerts using RDKit's built-in Brenk filter catalog
brenk_params = FilterCatalogParams()
brenk_params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
brenk_catalog = FilterCatalog.FilterCatalog(brenk_params)


def check_pains(smile):
    """Check for presence of PAINS alert in molecule.

    :param smile: molecule to check
    :return: boolean indicating if any PAINS alerts are present
    """
    mol = get_mol_from_smiles(smile)
    if mol is None:
        return 0
    return 1 if pains_catalog.HasMatch(mol) else 0


def check_brenk(smile):
    """Check for presence of Brenk alert in molecule.

    :param smile: molecule to check
    :return: boolean indicating if any Brenk alerts are present
    """
    mol = get_mol_from_smiles(smile)
    if mol is None:
        return 0
    return 1 if brenk_catalog.HasMatch(mol) else 0






input_len = len(smiles_list)
output_len = len(outputs)
assert input_len == output_len

with open(output_file, "w") as f:
    writer = csv.writer(f)
    writer.writerow(["value"])  # header
    for o in outputs:
        writer.writerow([o])
