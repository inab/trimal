/* *****************************************************************************

    trimAl v2.0: a tool for automated alignment trimming in large-scale
                 phylogenetics analyses.

    readAl v2.0: a tool for automated alignment conversion among different
                 formats.

    2009-2019
        Fernandez-Rodriguez V.  (victor.fernandez@bsc.es)
        Capella-Gutierrez S.    (salvador.capella@bsc.es)
        Gabaldon, T.            (tgabaldon@crg.es)

    This file is part of trimAl/readAl.

    trimAl/readAl are free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, the last available version.

    trimAl/readAl are distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with trimAl/readAl. If not, see <http://www.gnu.org/licenses/>.

***************************************************************************** */
#ifndef VALUES
#define VALUES

#include <string>

/* Characters used for different alignments type */

const std::string nucleotideResidues = "ACGTU";

const std::string aminoAcidResidues = "ARNDCQEGHILKMFPSTWYV";

const std::string degenerateNucleotideResidues = "ACGTURYKMSWBDHV";

/* Characters used to indicate indeterminations */
const std::string ambiguousAA = "BJXZ*";

/* Pyrrolysine:    'O' > 'TAG'  */
/* Selenocysteine: 'U' > 'TGA'  */
const std::string alternativeAminoAcidResidues = "UO";


#ifdef SIMMatrix
/* Default Identity Matrix for Canonical Nucleotides */
const std::array<std::array<float, 5>, 5> defaultNTMatrix = {{
  {1, 0, 0, 0, 0},
  {0, 1, 0, 0, 0},
  {0, 0, 1, 0, 0},
  {0, 0, 0, 1, 0},
  {0, 0, 0, 0, 1}
}};

const std::array<std::array<float, 15>, 15> defaultNTDegeneratedMatrix = {{
/* A: adenosine (A)        C: cytidine   (C)            G: guanine (G)            T: thymidine  (T)           U: uridine    (U)
 * R: purine    (G | A)    Y: pyrimidine (C | T/u)      K: keto    (G | T/u)      M: amino      (A | C)       S: strong     (G | C)
 * W: weak      (A | T/u)  B: not A      (G | C | T/u)  D: not C   (G | A | T/u)  H: not G      (A | C | T/u) V: not T/u    (G | C | A) */
  { 1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0},
  { 0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0},
  { 0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0},
  { 0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0},
  { 0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0},
  { 1/4., 0.0,  1/4., 0.0,  0.0,  1/4., 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0},
  { 0.0,  1/4., 0.0,  1/4., 1/4., 0.0,  1/4., 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0},
  { 0.0,  0.0,  1/4., 1/4., 1/4., 0.0,  0.0,  1/4., 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0},
  { 1/4., 1/4., 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1/4., 0.0,  0.0,  0.0,  0.0,  0.0,   0.0},
  { 0.0,  1/4., 1/4., 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1/4., 0.0,  0.0,  0.0,  0.0,   0.0},
  { 1/4., 0.0,  0.0,  1/4., 1/4., 0.0,  0.0,  0.0,  0.0,  0.0,  1/4., 0.0,  0.0,  0.0,   0.0},
  { 0.0,  1/6., 1/6., 1/6., 1/6., 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1/6., 0.0,  0.0,   0.0},
  { 1/6., 0.0,  1/6., 1/6., 1/6., 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1/6., 0.0,   0.0},
  { 1/6., 1/6., 0.0,  1/6., 1/6., 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1/6.,  0.0},
  { 1/6., 1/6., 1/6., 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1/6.}
}};

/* BLOSUM62 Similarity Matrix */
const std::array<std::array<float, 20>, 20> defaultAAMatrix = {{
  /* A: Alanine      R: Arginine      N: Asparagine      D: Aspartate      C: Cysteine       Q: Glutamine      E: Glutamate  
     G: Glycine      H: Histidine     I: Isoleucine      L: Leucine        K: Lysine         M: Methionine     F: Phenylalanine
     P: Proline      S: Serine        T: Threonine       W: Tryptophan     Y: Tyrosine       V: Valine
  */
  {  4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0},
  { -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3},
  { -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3},
  { -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3},
  {  0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1},
  { -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2},
  { -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2},
  {  0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3},
  { -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3},
  { -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3},
  { -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1},
  { -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2},
  { -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1},
  { -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1},
  { -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2},
  {  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2},
  {  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0},
  { -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3},
  { -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1},
  {  0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4}
}};


/* Alternative matrixes */

// Nucleotides
const std::array<std::array<float, 15>, 15> alternative_1_NTDegeneratedMatrix = {{
  { 1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
  { 0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
  { 0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
  { 0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
  { 0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
  { 0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
  { 0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0},
  { 0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0},
  { 0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0},
  { 0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0},
  { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0},
  { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0},
  { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0},
  { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0},
  { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1}
}};

/* BLOSUM62 Similarity Matrix extended with ambiguous amino acids */
const std::array<std::array<float, 25>, 25> extendedAAMatrix = {{
  /* A: Alanine      R: Arginine      N: Asparagine      D: Aspartate      C: Cysteine       Q: Glutamine      E: Glutamate  
     G: Glycine      H: Histidine     I: Isoleucine      L: Leucine        K: Lysine         M: Methionine     F: Phenylalanine
     P: Proline      S: Serine        T: Threonine       W: Tryptophan     Y: Tyrosine       V: Valine         B: Aspartate or Asparagine
     J: Leucine or Isoleucine         Z: Glutamine or Glutamate            X: Any            *: Stop codon
  */
  {  4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1, -1, -1, -4},
  { -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1, -2,  0, -1, -4},
  { -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  4, -3,  0, -1, -4},
  { -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4, -3,  1, -1, -4},
  {  0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -1, -3, -1, -4},
  { -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0, -2,  4, -1, -4},
  { -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1, -3,  4. -1, -4},
  {  0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -4, -2, -1, -4},
  { -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0, -3,  0, -1, -4},
  { -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3,  3, -3, -1, -4},
  { -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4,  3, -3, -1, -4},
  { -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0, -3,  1, -1, -4},
  { -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3,  2, -1, -1, -4},
  { -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3,  0, -3, -1, -4},
  { -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -3, -1, -1, -4},
  {  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0, -2,  0, -1, -4},
  {  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1, -1, -1, -4},
  { -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -2, -2, -1, -4},
  { -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -1, -2, -1, -4},
  {  0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3,  2, -2, -1, -4},
  { -2, -1,  4,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4, -3,  0, -1, -4},
  { -1, -2, -3, -3, -1, -2, -3, -4, -3,  3,  3, -3,  2,  0, -3, -2, -1, -2, -1,  2, -3,  3, -3, -1, -4},
  { -1,  0,  0,  1, -3,  4,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -2, -2, -2,  0, -3,  4, -1, -4},
  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -4},
  { -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1}
}};

#endif
#endif