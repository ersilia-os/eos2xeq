import os
import csv
import sys
import pandas as pd
import warnings
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")
warnings.filterwarnings("ignore", category=RuntimeWarning, module="importlib._bootstrap")

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
    mol = get_mol_from_smiles(smile)
    if mol is None:
        return 0
    return 1 if pains_catalog.HasMatch(mol) else 0

def check_brenk(smile):
    mol = get_mol_from_smiles(smile)
    if mol is None:
        return 0
    return 1 if brenk_catalog.HasMatch(mol) else 0

has_pains = [check_pains(s) for s in smiles_list]
has_brenk = [check_brenk(s) for s in smiles_list]

# Check tanimoto similarity to known antibioitics
def get_fingerprint(smile, radius=2, nBits=2048):
    mol = get_mol_from_smiles(smile)
    if mol:
        return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits)
    else:
        return None

def max_tanimoto(fp, known_fps):
    if fp is None:
        return 0.0
    sims = DataStructs.BulkTanimotoSimilarity(fp, known_fps)
    return max(sims) if sims else 0.0

SMILES_COL = "SMILES"

known_abx = pd.read_csv(os.path.join(root, "..", "..", "checkpoints", "559_known_abx.csv"))

known_abx["fp"] = known_abx[SMILES_COL].apply(get_fingerprint)
known_fps = list(known_abx["fp"].dropna())

tan_similarity = [max_tanimoto(get_fingerprint(s), known_fps) for s in smiles_list]

TANIMOTO_THRESHOLD = 0.5

is_sim_known_ab = [1 if sim >= TANIMOTO_THRESHOLD else 0 for sim in tan_similarity]

# Substructure matches
NITROFURAN = "O=[N+](O)c1ccco1"
FLUOROQUINOLONE = "O=C(O)c2c[nH]c1ccc(F)cc1c2=O"
CARBEPENEM = "O=C(O)C1=CCC2CC(=O)N12"
BETALACTAM = "O=C1CCN1"

mols = [get_mol_from_smiles(s) for s in smiles_list]

normalizer = rdMolStandardize.Normalizer()
tautomer_enumerator = rdMolStandardize.TautomerEnumerator()
uncharger = rdMolStandardize.Uncharger()

def standardize(mol):
    if mol is None:
        return None
    try:
        mol = rdMolStandardize.Cleanup(mol)           # disconnect metals, normalize nitro, etc.
        mol = normalizer.normalize(mol)               # functional group normalization
        mol = uncharger.uncharge(mol)                 # neutralize charges
        mol = tautomer_enumerator.Canonicalize(mol)   # pick one tautomer
        Chem.SanitizeMol(mol)
        return mol
    except:
        return mol

mols = [standardize(mol) if standardize(mol) is not None else mol for mol in mols]

nitrofuran_query = standardize(Chem.MolFromSmiles(NITROFURAN))
fluoroquinolone_query = standardize(Chem.MolFromSmiles(FLUOROQUINOLONE))
carbepenem_query = standardize(Chem.MolFromSmiles(CARBEPENEM))
betalactam_query = standardize(Chem.MolFromSmiles(BETALACTAM))

nitrofuran_motif = [1 if mol and mol.HasSubstructMatch(nitrofuran_query) else 0 for mol in mols]
fluoroquinolone_motif = [1 if mol and mol.HasSubstructMatch(fluoroquinolone_query) else 0 for mol in mols]
carbepenem_motif = [1 if mol and mol.HasSubstructMatch(carbepenem_query) else 0 for mol in mols]
betalactam_motif = [1 if mol and mol.HasSubstructMatch(betalactam_query) else 0 for mol in mols]

# Put all outputs toghether

outputs = [
    has_pains,
    has_brenk,
    is_sim_known_ab,
    nitrofuran_motif,
    fluoroquinolone_motif,
    carbepenem_motif,
    betalactam_motif
]

outputs = list(map(list, zip(*outputs)))

input_len = len(smiles_list)
output_len = len(outputs)
assert input_len == output_len

with open(output_file, "w") as f:
    writer = csv.writer(f)
    writer.writerow(["has_pains", "has_brenk", "is_sim_known_ab", "nitrofuran_motif", "fluoroquinolone_motif", "carbepenem_motif", "betalactam_motif"])
    for o in outputs:
        writer.writerow(o)
