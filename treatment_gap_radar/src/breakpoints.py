"""CLSI M100 MIC breakpoints (curated subset) -> S/I/R for MIC-only datasets.

Values are (S_max, R_min) in ug/mL:  S if MIC <= S_max ; R if MIC >= R_min ; else I.
- S_max = None  -> no susceptible category (e.g. colistin: I<=2, R>=4).
- R_min = INF   -> no resistant category defined (only S/non-S); never called R.
CAVEAT: This is an approximate, simplified application of CLSI breakpoints for surveillance
analytics (non-meningitis/parenteral variants chosen; censored MICs treated at their bound).
It is NOT for clinical use. ATLAS uses its own native interpretation and bypasses this table.
"""
import numpy as np

INF = float("inf")

# pathogen -> organism group key
def group_for(pathogen, gram):
    if pathogen is None:
        return None
    p = pathogen.lower()
    if p.startswith("pseudomonas"):
        return "PSA"
    if p.startswith("acinetobacter"):
        return "ACI"
    if p.startswith("staphylococcus"):
        return "STA"
    if p.startswith("enterococcus"):
        return "ENC"
    if p == "streptococcus pneumoniae":
        return "SPN"
    if p.startswith("streptococcus"):
        return "STR"
    if p == "haemophilus influenzae":
        return "HINF"
    if p == "mycobacterium tuberculosis":
        return "MTB"
    # default: Enterobacterales (E.coli, Klebsiella, Enterobacter, Serratia, Citrobacter,
    # Proteus, Providencia, Salmonella, ...)
    if gram == "negative":
        return "ENT"
    return None

# (S_max, R_min) by group, canonical drug name
BP = {
 "ENT": {
   "Ampicillin": (8, 32), "Amoxicillin-clavulanate": (8, 32), "Ampicillin-sulbactam": (8, 32),
   "Piperacillin-tazobactam": (16, 128), "Ticarcillin-clavulanate": (16, 128),
   "Ceftriaxone": (1, 4), "Cefotaxime": (1, 4), "Ceftazidime": (4, 16), "Cefepime": (2, 16),
   "Cefixime": (1, 4), "Cefpodoxime": (2, 8), "Cefuroxime": (8, 32), "Cefoxitin": (8, 32),
   "Cefaclor": (8, 32), "Cefdinir": (1, 4), "Aztreonam": (4, 16),
   "Ertapenem": (0.5, 2), "Imipenem": (1, 4), "Meropenem": (1, 4), "Doripenem": (1, 4),
   "Meropenem-vaborbactam": (4, 16), "Imipenem-relebactam": (1, 4),
   "Ceftazidime-avibactam": (8, 16), "Ceftolozane-tazobactam": (2, 8), "Aztreonam-avibactam": (8, 16),
   "Cefiderocol": (4, 16),
   "Ciprofloxacin": (0.25, 1), "Levofloxacin": (0.5, 2), "Moxifloxacin": (0.25, 1),
   "Gentamicin": (2, 8), "Amikacin": (16, 64), "Tobramycin": (2, 8),
   "Trimethoprim-sulfamethoxazole": (2, 4), "Chloramphenicol": (8, 32),
   "Tetracycline": (4, 16), "Doxycycline": (4, 16), "Minocycline": (4, 16), "Tigecycline": (2, 8),
   "Colistin": (None, 4), "Polymyxin B": (None, 4), "Azithromycin": (16, 32),
 },
 "PSA": {
   "Piperacillin-tazobactam": (16, 128), "Ceftazidime": (8, 32), "Cefepime": (8, 32),
   "Aztreonam": (8, 32), "Imipenem": (2, 8), "Meropenem": (2, 8), "Doripenem": (2, 8),
   "Imipenem-relebactam": (2, 8), "Ceftazidime-avibactam": (8, 16), "Ceftolozane-tazobactam": (4, 16),
   "Cefiderocol": (4, 16), "Ciprofloxacin": (0.5, 2), "Levofloxacin": (1, 4),
   "Gentamicin": (4, 16), "Amikacin": (16, 64), "Tobramycin": (4, 16),
   "Colistin": (None, 4), "Polymyxin B": (None, 4), "Minocycline": (4, 16),
 },
 "ACI": {
   "Ampicillin-sulbactam": (8, 32), "Sulbactam": (4, 16), "Piperacillin-tazobactam": (16, 128),
   "Ceftazidime": (8, 32), "Cefepime": (8, 32), "Imipenem": (2, 8), "Meropenem": (2, 8),
   "Ciprofloxacin": (1, 4), "Levofloxacin": (2, 8), "Gentamicin": (4, 16), "Amikacin": (16, 64),
   "Tobramycin": (4, 16), "Trimethoprim-sulfamethoxazole": (2, 4),
   "Tetracycline": (4, 16), "Doxycycline": (4, 16), "Minocycline": (4, 16), "Tigecycline": (2, 8),
   "Colistin": (None, 4), "Polymyxin B": (None, 4), "Cefiderocol": (4, 16),
 },
 "STA": {
   "Oxacillin": (2, 4), "Penicillin": (0.12, 0.25), "Ceftaroline": (1, 4),
   "Ciprofloxacin": (1, 4), "Levofloxacin": (1, 4), "Moxifloxacin": (0.5, 2),
   "Erythromycin": (0.5, 8), "Azithromycin": (2, 8), "Clarithromycin": (2, 8), "Clindamycin": (0.5, 4),
   "Gentamicin": (4, 16), "Tetracycline": (4, 16), "Doxycycline": (4, 16), "Minocycline": (4, 16),
   "Omadacycline": (0.5, INF), "Tigecycline": (0.5, INF),
   "Vancomycin": (2, 16), "Teicoplanin": (8, 32), "Linezolid": (4, 8), "Daptomycin": (1, 2),
   "Trimethoprim-sulfamethoxazole": (2, 4), "Quinupristin-dalfopristin": (1, 4), "Rifampicin": (1, 4),
 },
 "ENC": {
   "Ampicillin": (8, 16), "Penicillin": (8, 16), "Vancomycin": (4, 32), "Teicoplanin": (8, 32),
   "Linezolid": (2, 8), "Daptomycin": (4, INF), "Ciprofloxacin": (1, 4), "Levofloxacin": (2, 8),
   "Tetracycline": (4, 16), "Doxycycline": (4, 16), "Minocycline": (4, 16),
   "Tigecycline": (0.25, INF), "Omadacycline": (0.12, INF), "Quinupristin-dalfopristin": (1, 4),
   "Gentamicin": (4, 16),
 },
 "SPN": {
   "Penicillin": (2, 8), "Amoxicillin": (2, 8), "Amoxicillin-clavulanate": (2, 8),
   "Ceftriaxone": (1, 4), "Cefotaxime": (1, 4), "Cefuroxime": (1, 4), "Cefaclor": (1, 4),
   "Cefdinir": (0.5, 2), "Cefpodoxime": (0.5, 2), "Cefepime": (1, 4),
   "Erythromycin": (0.25, 1), "Azithromycin": (0.5, 2), "Clarithromycin": (0.25, 1),
   "Clindamycin": (0.25, 1), "Levofloxacin": (2, 8), "Moxifloxacin": (1, 4),
   "Tetracycline": (2, 8), "Doxycycline": (0.25, 1), "Trimethoprim-sulfamethoxazole": (0.5, 4),
   "Vancomycin": (1, INF), "Linezolid": (2, INF),
 },
 "STR": {
   "Penicillin": (0.12, INF), "Ampicillin": (0.25, INF), "Amoxicillin": (0.25, INF),
   "Ceftriaxone": (0.5, INF), "Cefotaxime": (0.5, INF), "Cefuroxime": (0.5, INF),
   "Erythromycin": (0.25, 1), "Azithromycin": (0.5, 2), "Clarithromycin": (0.25, 1),
   "Clindamycin": (0.25, 1), "Levofloxacin": (2, 8), "Moxifloxacin": (1, 4),
   "Tetracycline": (2, 8), "Vancomycin": (1, INF),
 },
 # M. tuberculosis: WHO/CLSI critical concentrations (MGIT-oriented, approximate).
 # Encoded as (S_max=crit/2, R_min=crit): MIC >= critical concentration => resistant.
 "MTB": {
   "Isoniazid": (0.05, 0.1), "Rifampicin": (0.25, 0.5), "Ethambutol": (2.5, 5.0),
   "Levofloxacin": (0.5, 1.0), "Ofloxacin": (1.0, 2.0), "Moxifloxacin": (0.125, 0.25),
   "Bedaquiline": (0.5, 1.0), "Clofazimine": (0.5, 1.0), "Capreomycin": (1.25, 2.5),
   "Kanamycin": (1.25, 2.5), "Amikacin": (0.5, 1.0),
 },
 "HINF": {
   "Ampicillin": (1, 4), "Amoxicillin-clavulanate": (4, 8), "Ampicillin-sulbactam": (2, 4),
   "Ceftriaxone": (2, INF), "Cefotaxime": (2, INF), "Cefuroxime": (4, 16), "Cefaclor": (8, 32),
   "Cefdinir": (1, INF), "Cefixime": (1, INF), "Cefpodoxime": (2, INF),
   "Azithromycin": (4, INF), "Clarithromycin": (8, 32),
   "Levofloxacin": (2, INF), "Moxifloxacin": (1, INF), "Ciprofloxacin": (1, 4),
   "Trimethoprim-sulfamethoxazole": (0.5, 4), "Tetracycline": (2, 8), "Chloramphenicol": (2, 8),
 },
}


def interpret(pathogen, gram, drug, mic, mic_op=""):
    """Return 'S'/'I'/'R' or None if no breakpoint/MIC."""
    if mic is None or (isinstance(mic, float) and np.isnan(mic)):
        return None
    grp = group_for(pathogen, gram)
    if grp is None or grp not in BP:
        return None
    bp = BP[grp].get(drug)
    if bp is None:
        return None
    s_max, r_min = bp
    # censored handling: ">x"/">=x" -> true MIC >= x (use x as lower bound for R);
    # "<=x"/"<x" -> true MIC <= x (use x as upper bound for S)
    if r_min is not None and mic >= r_min:
        return "R"
    if s_max is not None and mic <= s_max:
        return "S"
    if s_max is None:          # colistin-style: not >=R_min -> Intermediate
        return "I"
    return "I"
