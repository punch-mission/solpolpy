"""Constants used in code"""
import astropy.units as u

VALID_KINDS = {"mzp": [["Bm", "Bz", "Bp", "alpha"], ["Bm", "Bz", "Bp"]],
               "bpb": [["B", "pB", "alpha"], ["B", "pB"]],
               "btbr": [["Bt", "Br", "alpha"], ["Bt", "Br"]],
               "stokes": [["Bi", "Bq", "Bu"], ["Bi", "Bq", "Bu", "alpha"]],
               "bp3": [["B", "pB", "pBp", "alpha"], ["B", "pB", "pBp"]],
               "bthp": [["B", "theta", "p"]],
               "fourpol": [["B0", "B90", "B45", "B135"]],  # fourpol comes before npol to force it instead of npol
               "npol": [["B*"]],  # star indicates angle in degrees, as many as desired are supported
               }

# offset angles come from https://www.sciencedirect.com/science/article/pii/S0019103515003620?via%3Dihub
STEREOA_OFFSET_ANGLE = 45.8 * u.degree
STEREOB_OFFSET_ANGLE = -18 * u.degree
