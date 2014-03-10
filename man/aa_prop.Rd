\name{aa_prop}
\alias{aa_prop}
\docType{data}
\title{
Aminoacid properties
}
\description{
Properties of aminoacids. The original table contained 246 rows, but only 243 of them are described. I went deeper into this topic and the descriptions are missed for row: 244, 245 and 246. No big deal. 
}
\usage{data(aa_prop)}
\format{
  The format is:
 'data.frame':  246 obs. of  20 variables:
 $ A: num  0.61 0.77 0.37 0.357 52.6 ...
 $ R: num  0.6 0.72 0.84 0.529 109.1 ...
 $ N: num  0.06 0.55 0.97 0.463 75.7 ...
 $ D: num  0.46 0.65 0.97 0.511 68.4 ...
 $ C: num  1.07 0.65 0.84 0.346 68.3 ...
 $ Q: num  0 0.72 0.64 0.493 89.7 ...
 $ E: num  0.47 0.55 0.53 0.497 84.7 ...
 $ G: num  0.07 0.65 0.97 0.544 36.3 ...
 $ H: num  0.61 0.83 0.75 0.323 91.9 ...
 $ I: num  2.22 0.98 0.37 0.462 102 151 100 6.6 11 -2.26 ...
 $ L: num  1.53 0.83 0.53 0.365 102 145 108 20 15 -2.46 ...
 $ K: num  1.15 0.55 0.75 0.466 105.1 ...
 $ M: num  1.18 0.98 0.64 0.295 97.7 ...
 $ F: num  2.02 0.98 0.53 0.314 113.9 ...
 $ P: num  1.95 0.55 0.97 0.509 73.6 ...
 $ S: num  0.05 0.55 0.84 0.507 54.9 ...
 $ T: num  0.05 0.83 0.75 0.444 71.2 ...
 $ W: num  2.65 0.77 0.97 0.305 135.4 ...
 $ Y: num  1.88 0.83 0.84 0.42 116.2 ...
 $ V: num  1.32 0.98 0.37 0.386 85.1 ...
}
\details{
P1 Hydrophobicity index (Argos et al., 1982) Structural prediction of membrane-bound proteins E
ur. J. Biochem. 128, 565-575 (1982)

P2 Conformational parameter of beta-structure (Beghin-Dirkx, 1975) Une methode statistique simp
le de prediction des conformations proteiques Arch. Int. Physiol. Biochim. 83, 167-168 (1975)

P3 Conformational parameter of beta-turn (Beghin-Dirkx, 1975) Une methode statistique simple de
 prediction des conformations proteiques Arch. Int. Physiol. Biochim. 83, 167-168 (1975)

P4 Average flexibility indices (Bhaskaran-Ponnuswamy, 1988) Positional flexibilities of amino a
cid residues in globular proteins Int. J. Peptide Protein Res. 32, 241-255 (1988)

P5 Residue volume (Bigelow, 1967) On the average hydrophobicity of proteins and the relation be
tween it and protein structure J. Theor. Biol. 16, 187-211 (1967)

P6 Information value for accessibility; average fraction 35% (Biou et al., 1988) 

P7 Information value for accessibility; average fraction 23% (Biou et al., 1988) Secondary stru
cture prediction: combination of three different methods Protein Engineering 2, 185-191 (1988)

P8 Retention coefficient in TFA (Browne et al., 1982) The isolation of peptides by high-perform
ance liquid chromatography using predicted elution positions Anal. Biochem. 124, 201-208 (1982)

P9 Retention coefficient in HFBA (Browne et al., 1982) The isolation of peptides by high-perfor
mance liquid chromatography using predicted elution positions Anal. Biochem. 124, 201-208 (1982
)2)

P10 Transfer free energy to surface (Bull-Breese, 1974) Surface tension of amino acid solutions
: A hydrophobicity scale of the amino acid residues Arch. Biochem. Biophys. 161, 665-670 (1974)

P11 Apparent partial specific volume (Bull-Breese, 1974) Surface tension of amino acid solution
s: A hydrophobicity scale of the amino acid residues Arch. Biochem. Biophys. 161, 665-670 (1974
)4)

P12 Normalized frequency of alpha-helix (Burgess et al., 1974) Analysis of conformations of ami
no acid residues and prediction of backbone topography in proteins Isr. J. Chem. 12, 239-286 (1
974)

P13 Normalized frequency of extended structure (Burgess et al., 1974) Analysis of conformations
 of amino acid residues and prediction of backbone topography in proteins Isr. J. Chem. 12, 239
-286 (1974)

P14 Steric parameter (Charton, 1981) Protein folding and the genetic code: An alternative quant
itative model J. Theor. Biol. 91, 115-123 (1981)

P15 Polarizability parameter (Charton-Charton, 1982) The structural dependence of amino acid hy
drophobicity parameters J. Theor. Biol. 99, 629-644 (1982)

P16 The Chou-Fasman parameter of the coil conformation (Charton-Charton, 1983) 

P17 Average volume of buried residue (Chothia, 1975) Structural invariants in protein folding N
ature 254, 304-308 (1975)

P18  Residue accessible surface area in tripeptide (Chothia, 1976)	The nature of the accessible 
and buried surfaces in proteins	J. Mol. Biol. 105, 1-14 (1976)

P19 Residue accessible surface area in folded protein (Chothia, 1976) The nature of the accessi
ble and buried surfaces in proteins J. Mol. Biol. 105, 1-14 (1976)

P20 Proportion of residues 95% buried (Chothia, 1976) The nature of the accessible and buried s
urfaces in proteins J. Mol. Biol. 105, 1-14 (1976)

P21 Proportion of residues 100% buried (Chothia, 1976) The nature of the accessible and buried 
surfaces in proteins J. Mol. Biol. 105, 1-14 (1976)

P22 Normalized frequency of beta-turn (Chou-Fasman, 1978a) Empirical predictions of protein con
formation Ann. Rev. Biochem. 47, 251-276 (1978)

P23 Normalized frequency of alpha-helix (Chou-Fasman, 1978b) Prediction of the secondary struct
ure of proteins from their amino acid sequence Adv. Enzymol. 47, 45-148 (1978)

P24 Normalized frequency of beta-sheet (Chou-Fasman, 1978b) Prediction of the secondary structu
re of proteins from their amino acid sequence Adv. Enzymol. 47, 45-148 (1978)

P25 Normalized frequency of beta-turn (Chou-Fasman, 1978b) Prediction of the secondary structur
e of proteins from their amino acid sequence Adv. Enzymol. 47, 45-148 (1978)

P26	Normalized average hydrophobicity scales (Cid et al., 1992)	Hydrophobicity and structural c
lasses in proteins	Protein Engineering 5, 373-375 (1992)

P27	Partial specific volume (Cohn-Edsall, 1943)		"Protein, Amino Acid, and Peptides", Reinhold,
 New York (1943)

P28	Normalized frequency of beta-sheet (Crawford et al., 1973)	The reverse turn as a polypeptid
e conformation in globular proteins	Proc. Natl. Acad. Sci. USA 70, 538-542 (1973)

P29 Normalized frequency of turn (Crawford et al., 1973) The reverse turn as a polypeptide conf
ormation in globular proteins Proc. Natl. Acad. Sci. USA 70, 538-542 (1973)

P30 Size (Dawson, 1972)  In "The Biochemical Genetics of Man" (Brock, D.J.H. and Mayo, O.,

P31 Amino acid composition (Dayhoff et al., 1978a) Composition of proteins In "Atlas of Protein
 Sequence and Structure", Vol.5, Suppl.3 (Dayhoff,

P32 Relative mutability (Dayhoff et al., 1978b) A model of evolutionary change in proteins In "
Atlas of Protein Sequence and Structure", Vol.5, Suppl.3 (Dayhoff,

P33 Consensus normalized hydrophobicity scale (Eisenberg, 1984) Three-dimensional structure of 
membrane and surface proteins Ann. Rev. Biochem. 53, 595-623 (1984)

P34 Solvation free energy (Eisenberg-McLachlan, 1986) Solvation energy in protein folding and b
inding Nature 319, 199-203 (1986)

P35 Atom-based hydrophobic moment (Eisenberg-McLachlan, 1986) Solvation energy in protein foldi
ng and binding Nature 319, 199-203 (1986)

P36	Molecular weight (Fasman, 1976)		"Handbook of Biochemistry and Molecular Biology", 3rd ed.,
 Proteins -

P37	Hydrophobic parameter pi (Fauchere-Pliska, 1983)	Hydrophobic parameters pi of amino-acid si
de chains from the partitioning of N-acetyl-amino-acid amides	Eur. J. Med. Chem. 18, 369-375 (1
983)

P38	Graph shape index (Fauchere et al., 1988)	Amino acid side chain parameters for correlation 
studies in biology and pharmacology	Int. J. Peptide Protein Res. 32, 269-278 (1988)

P39 Normalized van der Waals volume (Fauchere et al., 1988) Amino acid side chain parameters fo
r correlation studies in biology and pharmacology Int. J. Peptide Protein Res. 32, 269-278 (198
8)

P40 STERIMOL length of the side chain (Fauchere et al., 1988) Amino acid side chain parameters 
for correlation studies in biology and pharmacology Int. J. Peptide Protein Res. 32, 269-278 (1
988)

P41 STERIMOL minimum width of the side chain (Fauchere et al., 1988) Amino acid side chain para
meters for correlation studies in biology and pharmacology Int. J. Peptide Protein Res. 32, 269
-278 (1988)

P42 STERIMOL maximum width of the side chain (Fauchere et al., 1988) Amino acid side chain para
meters for correlation studies in biology and pharmacology Int. J. Peptide Protein Res. 32, 269
-278 (1988)

P43 Number of hydrogen bond donors (Fauchere et al., 1988) Amino acid side chain parameters for
 correlation studies in biology and pharmacology Int. J. Peptide Protein Res. 32, 269-278 (1988
)8)

P44	Positive charge (Fauchere et al., 1988)	Amino acid side chain parameters for correlation st
udies in biology and pharmacology	Int. J. Peptide Protein Res. 32, 269-278 (1988)

P45	Negative charge (Fauchere et al., 1988)	Amino acid side chain parameters for correlation st
udies in biology and pharmacology	Int. J. Peptide Protein Res. 32, 269-278 (1988)

P46	Helix-coil equilibrium constant (Finkelstein-Ptitsyn, 1977)	Theory of protein molecule self
-organization. II. A comparison of calculated thermodynamic parameters of local secondary struc
tures with experiments	Biopolymers 16, 497-524 (1977)

P47	Partition coefficient (Garel et al., 1973)	Coefficients de partage d'aminoacides, nucleobas
es, nucleosides et nucleotides dans un systeme solvant salin	J. Chromatogr. 78, 381-391 (1973)

P48	Alpha-helix indices (Geisow-Roberts, 1980)	Amino acid preferences for secondary structure v
ary with protein class	Int. J. Biol. Macromol. 2, 387-389 (1980)

P49	Beta-strand indices (Geisow-Roberts, 1980)	Amino acid preferences for secondary structure v
ary with protein class	Int. J. Biol. Macromol. 2, 387-389 (1980)

P50	Hydrophobicity factor (Goldsack-Chalifoux, 1973)	Contribution of the free energy of mixing 
of hydrophobic side chains	J. Theor. Biol. 39, 645-651 (1973)

P51	Residue volume (Goldsack-Chalifoux, 1973)	Contribution of the free energy of mixing of hydr
ophobic side chains	J. Theor. Biol. 39, 645-651 (1973)

P52	Composition (Grantham, 1974)	Amino acid difference formula to help explain protein evolutio
n	Science 185, 862-864 (1974)

P53	Polarity (Grantham, 1974)	Amino acid difference formula to help explain protein evolution	S
cience 185, 862-864 (1974)

P54	Volume (Grantham, 1974)	Amino acid difference formula to help explain protein evolution	Sci
ence 185, 862-864 (1974)

P55 Partition energy (Guy, 1985) Amino acid side-chain partition energies and distribution of r
esidues in soluble proteins Biophys. J. 47, 61-70 (1985)

P56 Hydration number (Hopfinger, 1971), Cited by Charton-Charton (1982)  "Intermolecular Intera
ctions and Biomolecular Organizations", Wiley, New York (1977)

P57 Hydrophilicity value (Hopp-Woods, 1981) Prediction of protein antigenic determinants from a
mino acid sequecces Proc. Natl. Acad. Sci. USA 78, 3824-3828 (1981)

P58 Heat capacity (Hutchens, 1970) Heat capacities, absolute entropies, and entropies of format
ion of amino acids and related compounds In "Handbook of Biochemistry", 2nd ed. (Sober, H.A., e
d.), Chemical Rubber Co., Cleveland, Ohio, pp. B60-B61 (1970)

P59 Absolute entropy (Hutchens, 1970) Heat capacities, absolute entropies, and entropies of for
mation of amino acids and related compounds In "Handbook of Biochemistry", 2nd ed. (Sober, H.A.
, ed.), Chemical Rubber Co., Cleveland, Ohio, pp. B60-B61 (1970)

P60	Entropy of formation (Hutchens, 1970)	Heat capacities, absolute entropies, and entropies of
 formation of amino acids and related compounds	In "Handbook of Biochemistry", 2nd ed. (Sober, 
H.A., ed.), Chemical Rubber Co., Cleveland, Ohio, pp. B60-B61 (1970)

P61	Normalized relative frequency of alpha-helix (Isogai et al., 1980)	Characterization of mult
iple bends in proteins	Biopolymers 19, 1183-1210 (1980)

P62	Normalized relative frequency of extended structure (Isogai et al., 1980)	

P63 Normalized relative frequency of bend (Isogai et al., 1980) Characterization of multiple be
nds in proteins Biopolymers 19, 1183-1210 (1980)

P64 Normalized relative frequency of coil (Isogai et al., 1980) Characterization of multiple be
nds in proteins Biopolymers 19, 1183-1210 (1980)

P65 Average accessible surface area (Janin et al., 1978) Conformation of amino acid side-chains
 in proteins J. Mol. Biol. 125, 357-386 (1978)

P66 Percentage of buried residues (Janin et al., 1978) Conformation of amino acid side-chains i
n proteins J. Mol. Biol. 125, 357-386 (1978)

P67 Percentage of exposed residues (Janin et al., 1978) Conformation of amino acid side-chains 
in proteins J. Mol. Biol. 125, 357-386 (1978)

P68 Ratio of buried and accessible molar fractions (Janin, 1979) Surface and inside volumes in 
globular proteins Nature 277, 491-492 (1979)

P69 Transfer free energy (Janin, 1979) Surface and inside volumes in globular proteins Nature 2
77, 491-492 (1979)

P70	Hydrophobicity (Jones, 1975)	Amino acid properties and side-chain orientation in proteins: 
A cross correlation approach	J. Theor. Biol. 50, 167-183 (1975)

P71	Relative frequency of occurrence (Jones et al., 1992)	The rapid generation of mutation data
 matrices from protein sequences	CABIOS 8, 275-282 (1992)

P72	Relative mutability (Jones et al., 1992)	The rapid generation of mutation data matrices fro
m protein sequences	CABIOS 8, 275-282 (1992)

P73	Amino acid distribution (Jukes et al., 1975)	Amino acid composition of proteins: Selection 
against the genetic code	Science 189, 50-51 (1975)

P74	Sequence frequency (Jungck, 1978)	The genetic code as a periodic table	J. Mol. Evol. 11, 21
1-224 (1978)

P75	Average relative probability of helix (Kanehisa-Tsong, 1980)	Local hydrophobicity stabilize
s secondary structures in proteins	Biopolymers 19, 1617-1628 (1980)

P76	Average relative probability of beta-sheet (Kanehisa-Tsong, 1980)	Local hydrophobicity stab
ilizes secondary structures in proteins	Biopolymers 19, 1617-1628 (1980)

P77	Flexibility parameter for no rigid neighbors (Karplus-Schulz, 1985)	Prediction of chain fle
xibility in proteins	Naturwiss. 72, 212-213 (1985)

P78 Flexibility parameter for one rigid neighbor (Karplus-Schulz, 1985) Prediction of chain fle
xibility in proteins Naturwiss. 72, 212-213 (1985)

P79 Flexibility parameter for two rigid neighbors (Karplus-Schulz, 1985) Prediction of chain fl
exibility in proteins Naturwiss. 72, 212-213 (1985)

P80 Net charge (Klein et al., 1984) Prediction of protein function from sequence properties: Di
scriminant Biochim. Biophys. Acta 787, 221-226 (1984)

P81 Side chain interaction parameter (Krigbaum-Rubin, 1971) Local interactions as structure det
erminant for globular proteins Biochim. Biophys. Acta 229, 368-383 (1971)

P82 Side chain interaction parameter (Krigbaum-Komoriya, 1979) Local interactions as a structur
e determinant for protein molecules:II Biochim. Biophys. Acta 576, 204-228 (1979)

P83 Fraction of site occupied by water (Krigbaum-Komoriya, 1979) Local interactions as a struct
ure determinant for protein molecules: II Biochim. Biophys. Acta 576, 204-228 (1979)

P84 Side chain volume (Krigbaum-Komoriya, 1979) Local interactions as a structure determinant f
or protein molecules: II Biochim. Biophys. Acta 576, 204-228 (1979)

P85 Hydropathy index (Kyte-Doolittle, 1982) A simple method for displaying the hydropathic char
acter of a protein J. Mol. Biol. 157, 105-132 (1982)

P86 Transfer free energy, CHP/water (Lawson et al., 1984) A simple experimental model for hydro
phobic interactions in proteins J. Biol. Chem. 259, 2910-2912 (1984)

P87 Hydrophobic parameter (Levitt, 1976) A simplified representation of protein conformations f
or rapid simulation of protein folfing J. Mol. Biol. 104, 59-107 (1976)

P88 Distance between C-alpha and centroid of side chain (Levitt, 1976) A simplified representat
ion of protein conformations for rapid simulation of protein folfing J. Mol. Biol. 104, 59-107 
(1976)

P89 Side chain torsion angle phi(AAAR) (Levitt, 1976) A simplified representation of protein co
nformations for rapid simulation of protein folfing J. Mol. Biol. 104, 59-107 (1976)

P90 Radius of gyration of side chain (Levitt, 1976) A simplified representation of protein conf
ormations for rapid simulation of protein folfing J. Mol. Biol. 104, 59-107 (1976)

P91 van der Waals parameter R0 (Levitt, 1976) A simplified representation of protein conformati
ons for rapid simulation of protein folfing J. Mol. Biol. 104, 59-107 (1976)

P92 van der Waals parameter epsilon (Levitt, 1976) A simplified representation of protein confo
rmations for rapid simulation of protein folfing J. Mol. Biol. 104, 59-107 (1976)

P93	Normalized frequency of alpha-helix, with weights (Levitt, 1978)	Conformational preferences
 of amino acids in globular proteins	Biochemistry 17, 4277-4285 (1978)

P94	Normalized frequency of beta-sheet, with weights (Levitt, 1978)	Conformational preferences 
of amino acids in globular proteins	Biochemistry 17, 4277-4285 (1978)

P95	Normalized frequency of reverse turn, with weights (Levitt, 1978)	Conformational preference
s of amino acids in globular proteins	Biochemistry 17, 4277-4285 (1978)

P96	Normalized frequency of alpha-helix, unweighted (Levitt, 1978)	Conformational preferences o
f amino acids in globular proteins	Biochemistry 17, 4277-4285 (1978)

P97	Normalized frequency of beta-sheet, unweighted (Levitt, 1978)	Conformational preferences of
 amino acids in globular proteins	Biochemistry 17, 4277-4285 (1978)

P98	Normalized frequency of reverse turn, unweighted (Levitt, 1978)	Conformational preferences 
of amino acids in globular proteins	Biochemistry 17, 4277-4285 (1978)

P99 Conformational preference for all beta-strands (Lifson-Sander, 1979) Antiparallel and paral
lel beta-strands differ in amino acid residue preference Nature 282, 109-111 (1979)

P100 Conformational preference for parallel beta-strands (Lifson-Sander, 1979) 

P101 Conformational preference for antiparallel beta-strands (Lifson-Sander, 1979) 

P102 Average surrounding hydrophobicity (Manavalan-Ponnuswamy, 1978) Hydrophobic character of a
mino acid residues in globular proteins Nature 275, 673-674 (1978)

P103 Normalized frequency of alpha-helix (Maxfield-Scheraga, 1976) Status of empirical methods 
for the prediction of protein backbone topography Biochemistry 15, 5138-5153 (1976)

P104 Normalized frequency of extended structure (Maxfield-Scheraga, 1976) Status of empirical m
ethods for the prediction of protein backbone topography Biochemistry 15, 5138-5153 (1976)

P105 Normalized frequency of left-handed alpha-helix (Maxfield-Scheraga, 1976) 

P106 Refractivity (McMeekin et al., 1964), Cited by Jones (1975)  In "Amino Acids and Serum Pro
teins" (Stekol, J.A., ed.), American Chemical Society, Washington, D.C., p. 54 (1964)

P107 Retention coefficient in HPLC, pH7.4 (Meek, 1980) Prediction of peptide retention times in
 high-pressure liquid chromatography on the basis of amino acid composition Proc. Natl. Acad. S
ci. USA 77, 1632-1636 (1980)

P108 Retention coefficient in HPLC, pH2.1 (Meek, 1980) Prediction of peptide retention times in
 high-pressure liquid chromatography on the basis of amino acid composition Proc. Natl. Acad. S
ci. USA 77, 1632-1636 (1980)

P109 Retention coefficient in NaClO4 (Meek-Rossetti, 1981) Factors affecting retention and reso
lution of peptides in high-performance liquid chromatography J. Chromatogr. 211, 15-28 (1981)

P110 Retention coefficient in NaH2PO4 (Meek-Rossetti, 1981) Factors affecting retention and res
olution of peptides in high-performance liquid chromatography J. Chromatogr. 211, 15-28 (1981)

P111 Average reduced distance for C-alpha (Meirovitch et al., 1980) Empirical studies of hydrop
hobicity. 1. Effect of protein size on the hydrophobic behavior of amino acids Macromolecules 1
3, 1398-1405 (1980)

P112 Average reduced distance for side chain (Meirovitch et al., 1980) Empirical studies of hyd
rophobicity. 1. Effect of protein size on the hydrophobic behavior of amino acids Macromolecule
s 13, 1398-1405 (1980)

P113	Average side chain orientation angle (Meirovitch et al., 1980)	Empirical studies of hydrop
hobicity. 1. Effect of protein size on the hydrophobic behavior of amino acids	Macromolecules 1
3, 1398-1405 (1980)

P114	Effective partition energy (Miyazawa-Jernigan, 1985)	Estimation of effective interresidue 
contact energies from protein crystal structures: Quasi-chemical approximation	Macromolecules 1
8, 534-552 (1985)

P115	Normalized frequency of alpha-helix (Nagano, 1973)	Local analysis of the mechanism of prot
ein folding. I. Prediction of helices, loops, and beta-structures from primary structure	J. Mol
. Biol. 75, 401-420 (1973)

P116 Normalized frequency of beta-structure (Nagano, 1973) 

P117 Normalized frequency of coil (Nagano, 1973) Local analysis of the mechanism of protein fol
ding. I. Prediction of helices, loops, and beta-structures from primary structure J. Mol. Biol.
 75, 401-420 (1973)

P118 AA composition of total proteins (Nakashima et al., 1990) Distinct character in hydrophobi
city of amino acid composition of mitochondrial proteins PROTEINS 8, 173-178 (1990)

P119 Sof AA composition of total proteins (Nakashima et al., 1990) 

P120 AA composition of membrane proteins (Nakashima et al., 1990) Distinct character in hydroph
obicity of amino acid composition of mitochondrial proteins PROTEINS 8, 173-178 (1990)

P121 Normalized composition of membrane proteins (Nakashima et al., 1990) Distinct character in
 hydrophobicity of amino acid composition of mitochondrial proteins PROTEINS 8, 173-178 (1990)

P122 Ratio of average and computed composition (Nakashima et al., 1990) Distinct character in h
ydrophobicity of amino acid composition of mitochondrial proteins PROTEINS 8, 173-178 (1990)

P123 8 A contact number (Nishikawa-Ooi, 1980) Prediction of the surface-interior diagram of glo
bular proteins by an empirical method Int. J. Peptide Protein Res. 16, 19-32 (1980)

P124 14 A contact number (Nishikawa-Ooi, 1986) Radial locations of amino acid residues in a glo
bular protein: Correlation with the sequence J. Biochem. 100, 1043-1047 (1986)

P125 Average non-bonded energy per atom (Oobatake-Ooi, 1977) An analysis of non-bonded energy o
f proteins J. Theor. Biol. 67, 567-584 (1977)

P126 Short and medium range non-bonded energy per atom (Oobatake-Ooi, 1977) An analysis of non-
bonded energy of proteins J. Theor. Biol. 67, 567-584 (1977)

P127 Long range non-bonded energy per atom (Oobatake-Ooi, 1977) An analysis of non-bonded energ
y of proteins J. Theor. Biol. 67, 567-584 (1977)

P128 Average non-bonded energy per residue (Oobatake-Ooi, 1977) An analysis of non-bonded energ
y of proteins J. Theor. Biol. 67, 567-584 (1977)

P129 Short and medium range non-bonded energy per residue (Oobatake-Ooi, 1977) An analysis of n
on-bonded energy of proteins J. Theor. Biol. 67, 567-584 (1977)

P130	Optimized beta-structure-coil equilibrium constant (Oobatake et al., 1985)	

P131	Optimized propensity to form reverse turn (Oobatake et al., 1985)	Optimization of amino ac
id parameters for correspondence of sequence	Bull. Inst. Chem. Res., Kyoto Univ. 63, 82-94 (198
5)

P132	Optimized transfer energy parameter (Oobatake et al., 1985)	Optimization of amino acid par
ameters for correspondence of sequence to tertiary structures of proteuins	Bull. Inst. Chem. Re
s., Kyoto Univ. 63, 82-94 (1985)

P133	Optimized average non-bonded energy per atom (Oobatake et al., 1985)	Optimization of amino
 acid parameters for correspondence of sequence to tertiary structures of proteuins	Bull. Inst.
 Chem. Res., Kyoto Univ. 63, 82-94 (1985)

P134 Optimized side chain interaction parameter (Oobatake et al., 1985) Optimization of amino a
cid parameters for correspondence of sequence to tertiary structures of proteuins Bull. Inst. C
hem. Res., Kyoto Univ. 63, 82-94 (1985)

P135 Normalized frequency of alpha-helix from LG (Palau et al., 1981) Protein secondary structu
re Int. J. Peptide Protein Res. 19, 394-401 (1981)

P136 Normalized frequency of alpha-helix from CF (Palau et al., 1981) Protein secondary structu
re Int. J. Peptide Protein Res. 19, 394-401 (1981)

P137 Normalized frequency of beta-sheet from LG (Palau et al., 1981) Protein secondary structur
e Int. J. Peptide Protein Res. 19, 394-401 (1981)

P138 Normalized frequency of beta-sheet from CF (Palau et al., 1981) Protein secondary structur
e Int. J. Peptide Protein Res. 19, 394-401 (1981)

P139 Normalized frequency of turn from LG (Palau et al., 1981) Protein secondary structure Int.
 J. Peptide Protein Res. 19, 394-401 (1981)

P140 Normalized frequency of turn from CF (Palau et al., 1981) Protein secondary structure Int.
 J. Peptide Protein Res. 19, 394-401 (1981)

P141	HPLC parameter (Parker et al., 1986)	New hydrophilicity scale derived from high-performanc
e liquid chromatography peptide retention data: Correlation of predicted surface residues with 
antigencity and x-ray-derived accessible sites	Biochemistry 25, 5425-5432 (1986)

P142 Partition coefficient (Pliska et al., 1981) Partition coefficients of amino acids and hydr
ophobic parameters pi of their side-chains as measured by thin-layer chromatography J. Chromato
gr. 216, 79-92 (1981)

P143 Surrounding hydrophobicity in folded form (Ponnuswamy et al., 1980) Hydrophobic packing an
d spatial arrangement of amino acid residues in globular proteins Biochim. Biophys. Acta 623, 3
01-316 (1980)

P144 Average gain in surrounding hydrophobicity (Ponnuswamy et al., 1980) Hydrophobic packing a
nd spatial arrangement of amino acid residues in globular proteins Biochim. Biophys. Acta 623, 
301-316 (1980)

P145 Average gain ratio in surrounding hydrophobicity (Ponnuswamy et al., 1980) 

P146 Surrounding hydrophobicity in alpha-helix (Ponnuswamy et al., 1980) Hydrophobic packing an
d spatial arrangement of amino acid residues in globular proteins Biochim. Biophys. Acta 623, 3
01-316 (1980)

P147 Surrounding hydrophobicity in beta-sheet (Ponnuswamy et al., 1980) Hydrophobic packing and
 spatial arrangement of amino acid residues in globular proteins Biochim. Biophys. Acta 623, 30
1-316 (1980)

P148 Surrounding hydrophobicity in turn (Ponnuswamy et al., 1980) Hydrophobic packing and spati
al arrangement of amino acid residues in globular proteins Biochim. Biophys. Acta 623, 301-316 
(1980)

P149 Accessibility reduction ratio (Ponnuswamy et al., 1980) Hydrophobic packing and spatial ar
rangement of amino acid residues in globular proteins Biochim. Biophys. Acta 623, 301-316 (1980
)0)

P150 Average number of surrounding residues (Ponnuswamy et al., 1980) Hydrophobic packing and s
patial arrangement of amino acid residues in globular proteins Biochim. Biophys. Acta 623, 301-
316 (1980)

P151 Intercept in regression analysis (Prabhakaran-Ponnuswamy, 1982) Shape and surface features
 of globular proteins Macromolecules 15, 314-320 (1982)

P152 Slope in regression analysis x 1.0E1 (Prabhakaran-Ponnuswamy, 1982) Shape and surface feat
ures of globular proteins Macromolecules 15, 314-320 (1982)

P153 Correlation coefficient in regression analysis (Prabhakaran-Ponnuswamy, 1982) Shape and su
rface features of globular proteins Macromolecules 15, 314-320 (1982)

P154	Hydrophobicity (Prabhakaran, 1990)	The distribution of physical, chemical and conformation
al properties in signal and nascent peptides	Biochem. J. 269, 691-696 (1990)

P155	Relative frequency in alpha-helix (Prabhakaran, 1990)	The distribution of physical, chemic
al and conformational properties in signal and nascent peptides	Biochem. J. 269, 691-696 (1990)

P156	Relative frequency in beta-sheet (Prabhakaran, 1990)	The distribution of physical, chemica
l and conformational properties in signal and nascent peptides	Biochem. J. 269, 691-696 (1990)

P157	Relative frequency in reverse-turn (Prabhakaran, 1990)	The distribution of physical, chemi
cal and conformational properties in signal and nascent peptides	Biochem. J. 269, 691-696 (1990
)0)

P158	Helix-coil equilibrium constant (Ptitsyn-Finkelstein, 1983)	Theory of protein secondary st
ructure and algorithm of its prediction	Biopolymers 22, 15-25 (1983)

P159 Beta-coil equilibrium constant (Ptitsyn-Finkelstein, 1983) Theory of protein secondary str
ucture and algorithm of its prediction Biopolymers 22, 15-25 (1983)

P160 Average reduced distance for C-alpha (Rackovsky-Scheraga, 1977) Hydrophobicity, hydrophili
city, and the radial and orientational distributions of residues in native proteins Proc. Natl.
 Acad. Sci. USA 74, 5248-5251 (1977)

P161 Average reduced distance for side chain (Rackovsky-Scheraga, 1977) Hydrophobicity, hydroph
ilicity, and the radial and orientational distributions of residues in native proteins Proc. Na
tl. Acad. Sci. USA 74, 5248-5251 (1977)

P162 Side chain orientational preference (Rackovsky-Scheraga, 1977) Hydrophobicity, hydrophilic
ity, and the radial and orientational distributions of residues in native proteins Proc. Natl. 
Acad. Sci. USA 74, 5248-5251 (1977)

P163 Information measure for alpha-helix (Robson-Suzuki, 1976) Conformational properties of ami
no acid residues in globular proteins J. Mol. Biol. 107, 327-356 (1976)

P164	Information measure for extended (Robson-Suzuki, 1976)	Conformational properties of amino 
acid residues in globular proteins	J. Mol. Biol. 107, 327-356 (1976)

P165	Information measure for pleated-sheet (Robson-Suzuki, 1976)	Conformational properties of a
mino acid residues in globular proteins	J. Mol. Biol. 107, 327-356 (1976)

P166 Information measure for extended without H-bond (Robson-Suzuki, 1976) Conformational prope
rties of amino acid residues in globular proteins J. Mol. Biol. 107, 327-356 (1976)

P167 Information measure for turn (Robson-Suzuki, 1976) Conformational properties of amino acid
 residues in globular proteins J. Mol. Biol. 107, 327-356 (1976)

P168 Information measure for coil (Robson-Suzuki, 1976) Conformational properties of amino acid
 residues in globular proteins J. Mol. Biol. 107, 327-356 (1976)

P169 Information measure for loop (Robson-Suzuki, 1976) Conformational properties of amino acid
 residues in globular proteins J. Mol. Biol. 107, 327-356 (1976)

P170 Hydration free energy (Robson-Osguthorpe, 1979) Refined models for computer simulation of 
protein folding: Applications to the study of conserved secondary structure and flexible hinge 
points during the folding of pancreatic trypsin inhibitor J. Mol. Biol. 132, 19-51 (1979)

P171 Mean area buried on transfer (Rose et al., 1985) Hydrophobicity of amino acid residues in 
globular proteins Science 229, 834-838 (1985)

P172	Mean fractional area loss (Rose et al., 1985)	Hydrophobicity of amino acid residues in glo
bular proteins	Science 229, 834-838 (1985)

P173	Side chain hydropathy, uncorrected for solvation (Roseman, 1988)	Hydrophilicity of polar a
mino acid side-chains is markedly reduced by flanking peptide bonds	J. Mol. Biol. 200, 513-522 
(1988)

P174	Side chain hydropathy, corrected for solvation (Roseman, 1988)	Hydrophilicity of polar ami
no acid side-chains is markedly reduced by flanking peptide bonds	J. Mol. Biol. 200, 513-522 (1
988)

P175	Loss of Side chain hydropathy by helix formation (Roseman, 1988)	Hydrophilicity of polar a
mino acid side-chains is markedly reduced by flanking peptide bonds	J. Mol. Biol. 200, 513-522 
(1988)

P176	Transfer free energy (Simon, 1976), Cited by Charton-Charton (1982)		Quantum Biochemistry 
and Specific Interactions, Abacus Press,

P177	Principal component I (Sneath, 1966)	Relations between chemical structure and biological a
ctivity in peptides	J. Theor. Biol. 12, 157-195 (1966)

P178 Principal component II (Sneath, 1966) Relations between chemical structure and biological 
activity in peptides J. Theor. Biol. 12, 157-195 (1966)

P179 Principal component III (Sneath, 1966) Relations between chemical structure and biological
 activity in  peptides J. Theor. Biol. 12, 157-195 (1966)

P180 Principal component IV (Sneath, 1966) Relations between chemical structure and biological 
activity in peptides J. Theor. Biol. 12, 157-195 (1966)

P181 Zimm-Bragg parameter s at 20 C (Sueki et al., 1984) Helix-coil stability constants for the
 naturally occurring amino acids  in water. 22. Histidine parameters from random  poly{(hydroxy
butyl)glutamine-co-L-histidine} Macromolecules 17, 148-155 (1984)

P182 Zimm-Bragg parameter sigma x 1.0E4 (Sueki et al., 1984) Helix-coil stability constants for
 the naturally occurring amino acids  in water. 22. Histidine parameters from random  poly{(hyd
roxybutyl)glutamine-co-L-histidine} Macromolecules 17, 148-155 (1984)

P183 Optimal matching hydrophobicity (Sweet-Eisenberg, 1983) Correlation of sequence hydrophobi
cities measures similarity in three-dimensional protein structure J. Mol. Biol. 171, 479-488 (1
983)

P184 Normalized frequency of alpha-helix (Tanaka-Scheraga, 1977) Statistical mechanical treatme
nt of protein conformation. 5.  A  multiphasic model for specific-sequence copolymers of amino 
acids Macromolecules 10, 9-20 (1977)

P185 Normalized frequency of extended structure (Tanaka-Scheraga, 1977) Statistical mechanical 
treatment of protein conformation. 5.  A multiphasic model for specific-sequence copolymers of 
amino acids Macromolecules 10, 9-20 (1977)

P186 Normalized frequency of coil (Tanaka-Scheraga, 1977) Statistical mechanical treatment of p
rotein conformation. 5.  A multiphasic model for specific-sequence copolymers of amino acids Ma
cromolecules 10, 9-20 (1977)

P187 Relative population of conformational state A (Vasquez et al., 1983) Computed conformation
al states of the 20 naturally occurring amino acid residues and of the prototype residue alpha-
aminobutyric acid Macromolecules 16, 1043-1049 (1983)

P188 Relative population of conformational state C (Vasquez et al., 1983) Computed conformation
al states of the 20 naturally occurring amino acid residues and of the prototype residue alpha-
aminobutyric acid Macromolecules 16, 1043-1049 (1983)

P189 Relative population of conformational state E (Vasquez et al., 1983) Computed conformation
al states of the 20 naturally occurring amino acid residues and of the prototype residue alpha-
aminobutyric acid Macromolecules 16, 1043-1049 (1983)

P190 Electron-ion interaction potential (Veljkovic et al., 1985) Is it possible to analyze DNA 
and protein sequences by the method of digital signal processing? IEEE Trans. Biomed. Eng. 32, 
337-341 (1985)

P191 Transfer free energy to lipophilic phase (von Heijne-Blomberg, 1979) Trans-membrane transl
ocation of proteins: The direct transfer model Eur. J. Biochem. 97, 175-181 (1979)

P192 Average interactions per side chain atom (Warme-Morgan, 1978) A survey of amino acid side-
chain interactions in 21 proteins J. Mol. Biol. 118, 289-304 (1978)

P193 RF value in high salt chromatography (Weber-Lacey, 1978) Genetic code correlations: Amino 
acids and their anticodon nucleotides J. Mol. Evol. 11, 199-210 (1978)

P194 Propensity to be buried inside (Wertz-Scheraga, 1978) Influence of water on protein struct
ure. An analysis of the preferences of amino acid residues for the inside or outside and for sp
ecific conformations in a protein molecule Macromolecules 11, 9-15 (1978)

P195 Free energy change of epsilon(i) to epsilon(ex) (Wertz-Scheraga, 1978) Influence of water 
on protein structure. An analysis of the preferences of amino acid residues for the inside or o
utside and for specific conformations in a protein molecule Macromolecules 11, 9-15 (1978)

P196 Free energy change of alpha(Ri) to alpha(Rh) (Wertz-Scheraga, 1978) Influence of water on 
protein structure. An analysis of the preferences of amino acid residues for the inside or outs
ide and for specific conformations in a protein molecule Macromolecules 11, 9-15 (1978)

P197 Free energy change of epsilon(i) to alpha(Rh) (Wertz-Scheraga, 1978) Influence of water on
 protein structure. An analysis of the preferences of amino acid residues for the inside or out
side and for specific conformations in a protein molecule Macromolecules 11, 9-15 (1978)

P198 Polar requirement (Woese, 1973) Evolution of genetic code Naturwiss. 60, 447-459 (1973)

P199 Hydration potential (Wolfenden et al., 1981) Affinties of amino acid side chains for solve
nt water Biochemistry 20, 849-855 (1981)

P200 Hydrophobicity (Zimmerman et al., 1968) The characterization of amino acid sequences in pr
oteins by statistical methods J. Theor. Biol. 21, 170-201 (1968)

P201 Bulkiness (Zimmerman et al., 1968) The characterization of amino acid sequences in protein
s by statistical methods J. Theor. Biol. 21, 170-201 (1968)

P202	Polarity (Zimmerman et al., 1968)	The characterization of amino acid sequences in proteins
 by statistical methods	J. Theor. Biol. 21, 170-201 (1968)

P203	Isoelectric point (Zimmerman et al., 1968)	The characterization of amino acid sequences in
 proteins by statistical methods	J. Theor. Biol. 21, 170-201 (1968)

P204 RF rank (Zimmerman et al., 1968) The characterization of amino acid sequences in proteins 
by statistical methods J. Theor. Biol. 21, 170-201 (1968)

P205 Delta G values for the peptides extrapolated to 0 M urea (O'Neil and DeGrado, 1990) A ther
modynamic scale for the helix-forming tendencies of the commonly occurring amino acids. Science
 250, 646-651 (1990)

P206 Helix formation parameters (delta delta G) (O'Neil and DeGrado, 1990) A thermodynamic scal
e for the helix-forming tendencies of the commonly occurring amino acids. Science 250, 646-651 
(1990)

P207 Normalized flexibility parameters (B-values), average (Vihinen et al., 1994) Accuracy of p
rotein flexibility predictions. Proteins 19, 141-149 (1994)

P208 Normalized flexibility parameters (B-values) for each residue surrounded by 

P209 Normalized flexibility parameters (B-values) for each residue surrounded by 

P210 Normalized flexibility parameters (B-values) for each residue surrounded by 

P211 Free energy in alpha-helical conformation (Munoz and Serrano et al., 1994) Intrinsic secon
dary structure propensities of the amino acids, using statistical phi-psi matrices: comparison 
with experimental scales. Proteins 20, 301-311 (1994)

P212 Free energy in alpha-helical region (Munoz and Serrano et al., 1994) Intrinsic secondary s
tructure propensities of the amino acids, using statistical phi-psi matrices: comparison with e
xperimental scales. Proteins 20, 301-311 (1994)

P213 Free energy in beta-strand conformation (Munoz and Serrano et al., 1994) Intrinsic seconda
ry structure propensities of the amino acids, using statistical phi-psi matrices: comparison wi
th experimental scales. Proteins 20, 301-311 (1994)

P214 Free energy in beta-strand region (Munoz and Serrano et al., 1994) Intrinsic secondary str
ucture propensities of the amino acids, using statistical phi-psi matrices: comparison with exp
erimental scales. Proteins 20, 301-311 (1994)

P215 Free energy in beta-strand region (Munoz and Serrano et al., 1994) Intrinsic secondary str
ucture propensities of the amino acids, using statistical phi-psi matrices: comparison with exp
erimental scales. Proteins 20, 301-311 (1994)

P216 Free energies of transfer of AcWl-X-LL peptides from bilayer interface to 

P217 Number of codon(s) 

P218 Recognition factors 

P219 Hphob. / Eisenberg et al. 

P220 Hphob. / Abraham & Leo  

P221 Hphob. / Black 

P222 Hphob. / Bull & Breese 

P223 Hphob. / Miyazawa et al. 

P224 Hphob. / Roseman 

P225 Hphob. HPLC / Wilson & al 

P226 Normalized frequency of alpha region(Maxfield-Scheraga, 1976) 

P227 Hphob. HPLC pH3.4 / Cowan 

P228 Hphob. HPLC pH7.5 / Cowan 

P229 Hphob. Rf mobility / Aboderin 

P230 \% buried residues/Janin 

P231 \% accessible residues/Janin 

P232 alpha-helix / Deleage & Roux 

P233 beta-sheet / Deleage & Roux 

P234 beta-turn / Deleage & Roux 

P235 Coil / Deleage & Roux 

P236 A.A. composition/McCaldon 

P237 A.A. comp. in SWISS-PROT /Bairoch 

P238 E1 PCP New quantitative descriptors Venkatarajan, M.S., Braun W., 2001 New quantitative de
scriptors for amino acids based on multidimensional scaling of a large number of physical-chemi
cal properties J Mol Modeling 7,445-453 (2001) 

P239 E2 PCP New quantitative descriptors Venkatarajan, M.S., Braun W., 2001 New quantitative de
scriptors for amino acids based on multidimensional scaling of a large number of physical-chemi
cal properties J Mol Modeling 7,445-453 (2001)

P240 E3 PCP New quantitative descriptors Venkatarajan, M.S., Braun W., 2001 New quantitative de
scriptors for amino acids based on multidimensional scaling of a large number of physical-chemi
cal properties J Mol Modeling 7,445-453 (2001)

P241 E4 PCP New quantitative descriptors Venkatarajan, M.S., Braun W., 2001 New quantitative de
scriptors for amino acids based on multidimensional scaling of a large number of physical-chemi
cal properties J Mol Modeling 7,445-453 (2001)

P242 E5 PCP New quantitative descriptors Venkatarajan, M.S., Braun W., 2001 New quantitative de
scriptors for amino acids based on multidimensional scaling of a large number of physical-chemi
cal properties J Mol Modeling 7,445-453 (2001)

P243 Side-chain contribution to protein stablity delta delta G (Takano T & Yutani K) A new scal
e for side-chain contribution to protein stability based on the empirical stability analysis of
 mutant proteins Protein Engineering 14,525-528 (2001)

}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(aa_prop)
## how I created descriptions (tmp are copied original descriptions)
# writeLines(unlist(lapply(tmp, function(x) {
#   txt <- strsplit(x, "")[[1]]
#   inds <- c(seq(1, length(txt), by = 95), length(txt))
#   txts <- lapply(1:(length(inds) - 1), function(i)
#     paste0(txt[inds[i]:(inds[i + 1] - 1)], collapse = ""))
#   txts[length(txts)] <- paste0(txts[length(txts)], txt[length(txt)], "\n")
#   txts
# })), "tmp.txt")


}
\keyword{datasets}
