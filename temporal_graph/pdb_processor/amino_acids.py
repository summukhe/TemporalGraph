"""
This code base reports all supported amino acid residue types
and also supports their associated property index.
"""

__version__ = 1.0

__all__ = ['AminoAcid', 'get_amino', 'valid_amino_acids']


__amino_acids__ = {"A": "ALA", "C": "CYS", "D": "ASP",
                   "E": "GLU", "F": "PHE", "G": "GLY",
                   "H": "HIS", "I": "ILE", "K": "LYS",
                   "L": "LEU", "M": "MET", "N": "ASN",
                   "P": "PRO", "Q": "GLN", "R": "ARG",
                   "S": "SER", "T": "THR", "V": "VAL",
                   "W": "TRP", "Y": "TYR" }

__aa_mw__ =  {"A":  89.10, "C": 121.16, "D": 133.11,
              "E": 147.13, "F": 165.19, "G":  75.07,
              "H": 155.16, "I": 131.18, "K": 146.19,
              "L": 131.18, "M": 149.21, "N": 132.12,
              "P": 115.13, "Q": 146.15, "R": 174.20,
              "S": 105.09, "T": 119.12, "V": 117.15,
              "W": 204.23, "Y": 181.19}

__aa_pI__ =  {"A": 6.00, "C": 5.07, "D": 2.77,
              "E": 3.22, "F": 5.48, "G": 5.97,
              "H": 7.59, "I": 6.02, "K": 9.74,
              "L": 5.98, "M": 5.74, "N": 5.41,
              "P": 6.30, "Q": 5.65, "R": 10.76,
              "S": 5.68, "T": 5.60, "V": 5.96,
              "W": 5.89, "Y": 5.66}


# Average flexibility indices (Bhaskaran-Ponnuswamy, 1988)
__aa_flexibility__ = {"A": 0.357, "C": 0.346, "D": 0.511,
                    "E": 0.497, "F": 0.314, "G": 0.544,
                    "H": 0.323, "I": 0.462, "K": 0.466,
                    "L": 0.365, "M": 0.295, "N": 0.463,
                    "P": 0.509, "Q": 0.493, "R": 0.529,
                    "S": 0.507, "T": 0.444, "V": 0.386,
                    "W": 0.305, "Y": 0.420}


# Propensity to be buried inside (Wertz-Scheraga, 1978)
__aa_buried_propensity__ = {"A": 0.520, "C": 0.830, "D": 0.370,
                          "E": 0.380, "F": 0.870, "G": 0.410,
                          "H": 0.700, "I": 0.790, "K": 0.310,
                          "L": 0.770, "M": 0.760, "N": 0.420,
                          "P": 0.350, "Q": 0.350, "R": 0.490,
                          "S": 0.490, "T": 0.380, "V": 0.720,
                          "W": 0.860, "Y": 0.640}


# Hydropathy index (Kyte-Doolittle, 1982)
__aa_hydropathy__ = { "A":  1.800, "C":  2.500, "D": -3.500,
                      "E": -3.500, "F":  2.800, "G": -0.400,
                      "H": -3.200, "I":  4.500, "K": -3.900,
                      "L":  3.800, "M":  1.900, "N": -3.500,
                      "P": -1.600, "Q": -3.500, "R": -4.500,
                      "S": -0.800, "T": -0.700, "V":  4.200,
                      "W": -0.900, "Y": -1.300}


# Residue accessible surface area in tripeptide (Chothia, 1976)
__aa_sasa_free__ = {"A": 115.000, "C": 135.000, "D": 150.000,
                    "E": 190.000, "F": 210.000, "G":  75.000,
                    "H": 195.000, "I": 175.000, "K": 200.000,
                    "L": 170.000, "M": 185.000, "N": 160.000,
                    "P": 145.000, "Q": 180.000, "R": 225.000,
                    "S": 115.000, "T": 140.000, "V": 155.000,
                    "W": 255.000, "Y": 230.000}


# Residue accessible surface area in folded protein (Chothia, 1976)
__aa_sasa_folded__ = {"A": 25.000, "C": 19.000, "D": 50.000,
                      "E": 49.000, "F": 24.000, "G": 23.000,
                      "H": 43.000, "I": 18.000, "K": 97.000,
                      "L": 23.000, "M": 31.000, "N": 63.000,
                      "P": 50.000, "Q": 71.000, "R": 90.000,
                      "S": 44.000, "T": 47.000, "V": 18.000,
                      "W": 32.000, "Y": 60.000}

# Residue volume (Goldsack-Chalifoux, 1973)
__aa_volume__ = {"A":  88.300, "C": 112.400, "D": 110.800,
                 "E": 140.500, "F": 189.000, "G":  60.000,
                 "H": 152.600, "I": 168.500, "K": 175.600,
                 "L": 168.500, "M": 162.200, "N": 125.100,
                 "P": 122.200, "Q": 148.700, "R": 181.200,
                 "S":  88.700, "T": 118.200, "V": 141.400,
                 "W": 227.000, "Y": 193.000}


def valid_amino_acids(one_letter=True):
    if one_letter is True:
        return sorted(__amino_acids__.keys())
    else:
        return [__amino_acids__[a] for a in sorted(__amino_acids__.keys())]


class AminoAcid():
    def __init__(self, aa_type="G"):
        assert aa_type in __amino_acids__
        self.__aa = aa_type

    def __str__(self):
        return '%s' % __amino_acids__[self.__aa]

    def __eq__(self, other):
        if isinstance(other, AminoAcid):
            return self.__aa == other.__aa
        elif isinstance(other, str):
            return self.__aa == other
        else:
            return False

    def name(self, one_letter_code=False):
        if one_letter_code is True:
            return self.__aa
        else:
            return __amino_acids__[self.__aa]

    def molecular_weight(self):
        return __aa_mw__[self.__aa]

    def sasa_free(self):
        return __aa_sasa_free__[self.__aa]

    def sasa_folded(self):
        return __aa_sasa_folded__[self.__aa]

    def volume(self):
        return __aa_volume__[self.__aa]

    def flexibility_index(self):
        return __aa_flexibility__[self.__aa]

    def hydropathy_index(self):
        return __aa_hydropathy__[self.__aa]

    def buried_propensity(self):
        return __aa_buried_propensity__[self.__aa]

    def pIsoelectric(self):
        return __aa_pI__[self.__aa]


def get_amino(amino_name):
    assert isinstance(amino_name, str)
    amino_name = amino_name.strip()
    if len(amino_name) == 3:
        assert amino_name in __amino_acids__.values()
        for key, value in __amino_acids__.items():
            if amino_name == value:
                return AminoAcid(aa_type=key)
    else:
        assert len(amino_name) == 1 and amino_name in __amino_acids__
        return AminoAcid(aa_type=amino_name)
    raise KeyError("Invalid amino name: %s" % amino_name)

