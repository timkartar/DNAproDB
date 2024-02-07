/* HBPLUS - Hydrogen Bond Calculator PDB input v 3.2 */
/* Copyright I K Mcdonald, D T Jones, D Naylor and J M Thornton 1993
             All Rights Reserved
   Copies are free to Academic Users - contact mcdonald@uk.ac.ucl.bsm.bioc.bsm

 The publication of research using the Software must reference "McDonald
 IK, Naylor DN, Jones DT & Thornton J M (1993), 'HBPLUS', Computer Program,
 Department of Biochemistry and Molecular Biology, University College
 London." or successor references as defined by the authors.

 Unless informed otherwise, you should be an academic user and have sent a
 signed confidentiality agreement to the authors (address at the end of the
 hbplus.h source code).  If you do not have a copy of HBPLUS and would like to
 receive one (free to academic users), please detach the confidentiality
 agreement from the end of this document, sign it, and send to the address
 given.  Please allow other people in your department to use your copy of
 HBPLUS, but do not allow them to make their own copy. */

/* Contents */
/* 1. Copyright **/
/* 2. Version Log **/
/* 3. Tokens as abbreviations **/

/* Prototypes */
short iswater(char * resid);
int string_truncate(char *string,int max_length);
/* <-- RAL 3 Jul 2012 */

/* 4. Tokens as flags and customisables **/

/* 5. Globals as flags and customisables **/

/* 6. Globals as name arrays **/

/* 8. Vector/Matrix Routines */

/* 8.5 Globals as arrays **/

/* 8.7Minor useful functions **/

/* 9. Hydrogen Position Calculation Subroutines */
/*10. Hydrogen Position Calculation Routine */
/*11. PDBOUT */

/*12. Connectivity Calculation Subroutines **/

/*13. Connectivity Calculation Routine **/

/*13.2 Domain calculation Routines */

/*13.5 Adding new residues to the selection **/

/*14. find_hb */

/*15. inpdb_file **/

/*16. main */

/****************************************************************************/
/* 2. Version Log */

/* Version 2.23    IM 4th January 1994
                   Separating into header files

   Version 2.25    IM 21st January 1994, further separation.

   Version 2.26    IM 24th January 1994, create array of secondary structure
                   segments.
*/
/****************************************************************************/
/* 3. Tokens as abbreviations */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
/* can't find
#include <stdarg.h>
 */
#include "hbplus.h"

/*#define BSM /* If this is defined, then the code compiles expecting to files in their locations in their home laboratory */

#define TF1 "%6.2f%6.2f" /* P.A. Keller 1993 */

#define TF "%8.3f%8.3f%8.3f"
#define VXYZ(v) v.x, v.y, v.z

#define unit unit_vector
#define to from_a_to_b

#define SHORTID(d) atom[d].chnid,atom[d].aanum,atom[d].resnam,atom[d].atmnam

#define SQR(X)  ((X)*(X))

/* A list of common PDB record types... */

#define HEADER             1
#define COMPND             2
#define SOURCE             3
#define AUTHOR             4
#define REVDAT             5
#define REMARK             6
#define SEQRES             7
#define FTNOTE             8
#define HET                9
#define FORMUL             10
#define HELIX              11
#define CRYST1             12
#define ORIGX1             13
#define ORIGX2             14
#define ORIGX3             15
#define SCALE1             16
#define SCALE2             17
#define SCALE3             18
#define ATOM               19
#define TER                20
#define HETATM             21
#define CONECT             22
#define ENDENT             23
#define JRNL               24
#define TURN               25
/* Bug-fix. RAL 1 May 1997 --> */                
#define ENDMDL             26
/* <-- RAL 1 May 1997 */                


/****************************************************************************/
/* 4. Tokens as flags and customisables */

#ifdef BSM
#define SSTHEAD "/idata/sst/p"
#endif

short TOTNAA = 69 ; /* total no of amino acids recognized by neighbour   */
#define MAXCNRECS 5000  /* max no of CONECT records allowed in a pdbfile     */
#define MAXCNCOLS    5  /* max no of columns of data in a CONECT record      */
#define IGNLB        0  /* ignore bonds between ligands and proteins?       */

/* Amendment. RAL 15 Nov 2007 --> */
//#define ARKAA       11  /* number of residues with ARKIVED donors & acceptors*/
#define ARKAA       16  /* no. of residues with ARKIVED donors & acceptors*/
/* <--- RAL 15 Nov 2007 */
/*MOVE*/

/* Bug-fix. RAL 1 May 1997 --> */                
/* #define NTOKENS            25 */
#define NTOKENS            26 /* tokens in pdb files */
/* <-- RAL 1 May 1997 */                

#define CHNBRKFLAG   4  /* mode for determination of chain breaks (1-4)      */
/*1 =) PEPBND   2=) CADBND  3=) unknown atoms 4=) all criteria */

#define   CADBND  5.0 /* max CA:CA distance , if such is criteria */
#define   PEPBND  2.5 /* max  N-C distance */



/****************************************************************************/
/* 5. Globals as flags and customisables */
int ssnum=0;

float   SSDIST = 3.0; /* internal linkage only */
#ifdef BSM /*If this is an internal copy*/
short inputsstflg=1;
#else
short inputsstflg=0;
#endif
/* =0 don't worry about sst, =1 do try to input them */
short cssflg=0;
short kshflg=0;
char only_chainid_lst[32]="\0";

/****************************************************************************/
/* 6. Globals as name arrays */

/* Residue name to allow conversion of a.a. name into numeric code */
char rnames[MAXNAA][4] =
{
 "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
 "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR",
 "AIB", "PHL", "SEC", "ALM", "MPR", "FRD", "LYM", "GLM", "PPH", "PGL", 
 "OLE", "ABA", "NLE", "B2V", "B2I", "B1F", "BNO", "B2A", "B2F", "IVA", 
 "LOV", "STA", "PVL", "CAL", "PHA", "DCI", "AHS", "CHS", "MSE", "ETA",
 "PCA", "ASX", "GLX", "UNK", "CYH", "CSS",
 "  C", "  A", "  U", "  G", "  T",
/* Amendment. RAL 15 Nov 2007 --> */
 " DC", " DA", " DU", " DG", " DT",
/* <--- RAL 15 Nov 2007 */
 "ATP", "COA", "FMN", "HEM", "MTX", "NAD",
 "ACE", "FOR"
};

/* HETATMS
 "HOH", "  C","  A","  T","  G","ATP","COA","FMN","HEM","MTX","NAD","NDP"
*/

/* synonym table
"CSH",cyh
"CYX",css
"PRZ",pro
"TRY",tyr
"ILU",ile
"1MA",__a
"5MC",__c
"OMC",__c
"1MG",__g
"2MG",__g
"M2G",__g
"7MG",__g
"OMG",__g
" YG",_yg
"  I"     ?inositine
" +U",__u
"H2U",__u
"5MU",??
"PSU",__u
*/

/*****************************************************************************/
/* Record names for decoding record types */
/* Bug-fix. RAL 1 May 1997 --> */                
/* const char           *tokstr[25] = */
const char           *tokstr[NTOKENS] =
/* <-- RAL 1 May 1997 */                
{
 "HEADER", "COMPND", "SOURCE", "AUTHOR", "REVDAT",
 "REMARK", "SEQRES", "FTNOTE", "HET",    "FORMUL",
 "HELIX",  "CRYST1", "ORIGX1", "ORIGX2", "ORIGX3",
 "SCALE1", "SCALE2", "SCALE3", "ATOM",   "TER",
 "HETATM", "CONECT", "END",    "JRNL",   "TURN",
/* Bug-fix. RAL 1 May 1997 --> */                
 "ENDMDL"
/* <-- RAL 1 May 1997 */                

};

/* Bonding Atom Dictionary - this is the list of atom names, arranged on a
   residue by residue basis in RESCODE order. 
*/
char necatm[MAXNAA][TOTNATM*4] =
{
      " N   CA  C   O   CB  OXT",
      " N   CA  C   O   CB  SG  OXT",
      " N   CA  C   O   CB  CG  OD1 OD2 OXT",
      " N   CA  C   O   CB  CG  CD  OE1 OE2 OXT",
      " N   CA  C   O   CB  CG  CD1 CD2 CE1 CE2 CZ  OXT",
      " N   CA  C   O   C   OXT",
      " N   CA  C   O   CB  CG  ND1 CD2 CE1 NE2 OXT",
      " N   CA  C   O   CB  CG1 CG2 CD1 OXT",
      " N   CA  C   O   CB  CG  CD  CE  NZ  OXT",
      " N   CA  C   O   CB  CG  CD1 CD2 OXT",
      " N   CA  C   O   CB  CG  SD  CE  OXT",
      " N   CA  C   O   CB  CG  OD1 ND2 OXT",
      " N   CA  C   O   CB  CG  CD  OXT",
      " N   CA  C   O   CB  CG  CD  OE1 NE2 OXT",
      " N   CA  C   O   CB  CG  CD  NE  CZ  NH1 NH2 OXT",
      " N   CA  C   O   CB  OG  OXT",
      " N   CA  C   O   CB  OG1 CG2 OXT",
      " N   CA  C   O   CB  CG1 CG2 OXT",
      " N   CA  C   O   CB  CG  CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2 OXT",/*Trp*/
      " N   CA  C   O   CB  CG  CD1 CD2 CE1 CE2 CZ  OH  OXT",/*Tyr*/
      " N   CA  C   O   CB1 CB2 OXT",
      " N   CA  C   O   CB  CG  CD1 CD2 CE1 CE2 CZ  OXT",
      " N   CA  C   O   CB SEG  OD1 OD2 OXT",
      " N   CA  C   O   CB  CM  OXT",
      " ??? CA  C   O   CB  SG ",
      " N   CA  C   ??? CB  CG  CD1 CD2 CE1 CE2 CZ  OXT",
      " N   CA  C   O   CB  CG  CD  CE  NZ  CM  OXT",
      " N   CA  C   O   CM  OXT",
      " N   CA  P   OP1 OP2 CB  CG  CD1 CD2 CE1 CE2 CZ  OXT",
      " N   CA  P   O1P O2P OXT",
      " ON  CA  C   O   CB  CG  CD1 CD2 OXT",
      " N   CA  C   O   CB  CG  OXT",
      " N   CA  C   O   CB  CG  CD  CE  OXT",
      " N   CA  B   O1  O2  CB  CG1 CG2 OXT",
      " N   CA  B   O1  O2  CB  CG1 CG2 CD1 OXT",
      " N   CA  B   O1  O2  CB  CG  CD1 CD2 CE1 CE2 CZ  OXT",
      " N   CA  B   O1  O2  CB  CG  CD  CE  OXT",
      " N   CA  B   O1  O2  CB  OXT",
      " N   CA  B   O1  O2  CB  CG  CD1 CD2 CE1 CE2 CZ  OXT",
      " ??? CA  C   O   CB  CG1 CG2 OXT",
      " N   CA  C   O   CB  CG1 CG2 CD1 CD2 C1G C1B C1A CS  OS  CT  OXT",
      " N   CA  C   O   CB  CG  CD1 CD2 CH  OH  CM  OXT",
      " ON  CA  C   O   CB ",
      " N   CA  C   O   CB  CG  CD1 CD2 CE1 CE2 CZ  CH  OH  CM  CA2 CB2 CG2 CD3 CD4 OXT",
      " N   CA  C   O   CB  CG  CD1 CD2 CE1 CE2 CZ  OXT",
      " N   CA  ??? ??? CB  CG1 CG2 CD1",
      " N   CA  C   O   CB  CG  CD1 CD2 CE1 CE2 CZ  CH  OH  CM  N1  CB2 CG2 CD3 CD4 OXT",
      " N   CA  C   O   CB  CG  CD1 CD2 CE1 CE2 CZ  CH  OH  CM  OXT",
      " N   CA  C   O   CB  CG SED  CE  OXT",
      " N   CA  ??? O   CB ",
      " N   CA  C   O   CB  CG  CD  OE  OXT",
      " N   CA  C   O   CB  CG  AD1 AD2 OXT",/*asx*/
      " N   CA  C   O   CB  CG  CD  AE1 AE2 OXT", /*glx*/
      " N   CA  C   O   CB  OXT", /*unk - do edit*/
      " N   CA  C   O   CB  SG  OXT",/*cyh*/
      " N   CA  C   O   CB  SG  OXT",/*css*/
/* Amendment. RAL 3 Mar 2011 --> */
      " ??? OP3 P   OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N1  C2  O2  N3  C4  N4  C5  C6  OXT",  /*  C*/
      " ??? OP3 P   OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N9  C8  N7  C5  C6  N6  N1  C2  N3  C4  OXT",/*  A*/
      " ??? OP3 P   OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N1  C2  O2  N3  C4  O4  C5  C6  OXT",  /*  U*/
      " ??? OP3 P   OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N9  C8  N7  C5  C6  O6  N1  C2  N2  N3  C4  OXT",/*  G*/
      " ??? OP3 P   OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N1  C2  O2  N3  C4  O4  C5  C5M C6  OXT",  /*  T*/
//      " ??? O3P P   O1P O2P O5* C5* C4* O4* C3* O3* C2* O2* C1* N1  C2  O2  N3  C4  N4  C5  C6  OXT",  /*  C*/
//      " ??? O3P P   O1P O2P O5* C5* C4* O4* C3* O3* C2* O2* C1* N9  C8  N7  C5  C6  N6  N1  C2  N3  C4  OXT",/*  A*/
//      " ??? O3P P   O1P O2P O5* C5* C4* O4* C3* O3* C2* O2* C1* N1  C2  O2  N3  C4  O4  C5  C6  OXT",  /*  U*/
//      " ??? O3P P   O1P O2P O5* C5* C4* O4* C3* O3* C2* O2* C1* N9  C8  N7  C5  C6  O6  N1  C2  N2  N3  C4  OXT",/*  G*/
//      " ??? O3P P   O1P O2P O5* C5* C4* O4* C3* O3* C2* O2* C1* N1  C2  O2  N3  C4  O4  C5  C5M C6  OXT",  /*  T*/
/* <--- RAL RAL 3 Mar 2011 */
/* Amendment. RAL 15 Nov 2007 --> */
      " ??? OP3 P   OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N1  C2  O2  N3  C4  N4  C5  C6  OXT",  /*  New C */
      " ??? OP3 P   OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N9  C8  N7  C5  C6  N6  N1  C2  N3  C4  OXT",/*  New A */
      " ??? OP3 P   OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N1  C2  O2  N3  C4  O4  C5  C6  OXT",  /*  New U */
      " ??? OP3 P   OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N9  C8  N7  C5  C6  O6  N1  C2  N2  N3  C4  OXT",/*  New G */
      " ??? OP3 P   OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N1  C2  O2  N3  C4  O4  C5  C7  C6  OXT",  /*  T*/
/* <--- RAL 15 Nov 2007 */
/*       00  01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36  37  38  39  40  41  42  43  44  45  46  47  48  49  50  51*/

      " ??? PG  O1G O2G O3G PB  O1B O2B O3B PA  O1A O2A O3A O5* C5* C4* O4* C3* O3* C2* O2* C1* N9  C8  N7  C5  C6  N6  N1  C2  N3  C4  OXT",/*ATP*/
      " ???AO6 AP2 AO4 AO5 AO3 AP1 AO1 AO2 AO5*AC5*AC4*AO4*AC3*AO3*AP3*AO7 AO8 AO9 AC2*AO2*AC1*AN9 AC8 AN7 AC5 AC6 AN6 AN1 AC2 AN3 AC4 PS1 PC2 PC3 PN4 PC5 PO5 PC6 PC7 PN8 PC9 PO9 PC10PO10PC11PC12PC13PC14", /*COA*/
      " ??? N1  C2  O2  N3  C4  O4  C4A N5  C5A C6  C7  C7M C8  C8M C9  C9A N10 C10 C1* C2* O2* C3* O3* C4* O4* C5* O5* P   OP1 OP2 OP3 OXT", /*FMN*/
      " ???FE   CHA CHB CHC CHD N A C1A C2A C3A C4A CMA CAA CBA CGA O1A O2A N B C1B C2B C3B C4B CMB CAB CBB N C C1C C2C C3C C4C CMC CAC CBC N D C1D C2D C3D C4D CMD CAD CBD CGD O1D O2D", /*HEM*/
      " ??? N1  C2  NA2 N3  C4  NA4 C4A N5  C6  C7  N8  C8A C9  N10 CM  C11 C12 C13 C14 C15 C16 C   O   N   CA  CT  O1  O2  CB  CG  CD  OE1 OE2", /*MTX*/
      " ???AP  AO1 AO2 AO5*AC5*AC4*AO4*AC3*AO3*AC2*AO2*AC1*AN9 AC8 AN7 AC5 AC6 AN6 AN1 AC2 AN3 AC4  O3 NP  NO1 NO2 NO5*NC5*NC4*NO4*NC3*NO3*NC2*NO2*NC1*NN1 NC2 NNC3 NC7 NO7 NN7 NC4 NC5 NC6 ", /*NAD*/

      " ??? CH3 C   O   OXT", /*ACE *official MISC*/
      " ??? ??? C   O   OXT" /*FOR *official MISC*/
};
/*    " ??? CM1 C   O   CM2" ACN */
/*    "     CH3 C   O   O1 " ACT */

/* Number of permissible bonds for hydrogen bond acceptors */
short accepts[MAXNAA][TOTNATM] =
{
 {0, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*ala*/
 {0, 0, 0, 2, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*cys*/
 {0, 0, 0, 2, 0, 0, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 0, 0, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*gly*/
 {0, 0, 0, 2, 0, 0, 1, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*his*/
 {0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 0, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {2, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0},
 {2, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 2},
 {0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 2},
 {0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 0, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 2, 0, 0, 0, 2, 2, 2,-0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, /*glx*/
 {0, 0, 0, 2, 0,-0,-0,-0,-0,-0,-0,-0,-0,-0,-0,-0,-0,-0,-0,-0}, /*unk*/
 {0, 0, 0, 2, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*cys*/
 {0, 0, 0, 2, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*css*/
/*0,01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19*/
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*C*/
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*A*/
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*U*/
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*G*/
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*T*/
/* Amendment. RAL 15 Nov 2007 --> */
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*New C*/
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*New A*/
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*New U*/
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*New G*/
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*New T*/
/* <--- RAL 15 Nov 2007 */
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*ATP*/
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /*CoA*/
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, /*FMN*/
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, /*Hem*/
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, /*Mtx*/
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, /*NAD*/
 {0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*ACE*/
 {0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}/*FOR*/


};
/* listing them . . .*/

struct index
{
    short n;
    char s[128];
};

struct index accepts_archive[ARKAA] =
{
/* Amendment. RAL 3 Mar 2011 --> */
    {__C," OP1 OP2 O5' O4' O3' O2' O2  N3  OXT"},
    {__A," OP1 OP2 O5' O4' O3' O2' N1  N3  N7  OXT"},
    {__U," OP1 OP2 O5' O4' O3' O2' O2  O4  OXT"},
    {__G," OP1 OP2 O5' O4' O3' O2' N3  O6  N7  OXT"},
    {__T," OP1 OP2 O5' O4' O3' O2' O2  O4  OXT"},
//    {__C," O1P O2P O5* O4* O3* O2* O2  N3  OXT"},
//    {__A," O1P O2P O5* O4* O3* O2* N1  N3  N7  OXT"},
//    {__U," O1P O2P O5* O4* O3* O2* O2  O4  OXT"},
//    {__G," O1P O2P O5* O4* O3* O2* N3  O6  N7  OXT"},
//    {__T," O1P O2P O5* O4* O3* O2* O2  O4  OXT"},
/* <--- RAL RAL 3 Mar 2011 */

/* Amendment. RAL 15 Nov 2007 --> */
    {_DC," OP1 OP2 O5' O4' O3' O2' O2  N3  OXT"},
    {_DA," OP1 OP2 O5' O4' O3' O2' N1  N3  N7  OXT"},
    {_DU," OP1 OP2 O5' O4' O3' O2' O2  O4  OXT"},
    {_DG," OP1 OP2 O5' O4' O3' O2' N3  O6  N7  OXT"},
    {_DT," OP1 OP2 O5' O4' O3' O2' O2  O4  OXT"},
/* <--- RAL 15 Nov 2007 */

    {ATP," O1G O2G O3G O1B O2B O3B O1A O2A O3A O5* O4* O3* O2* N1  N3  N7  OXT"},
    {CoA,"AO6 AO4 AO5 AO3 AO1 AO2 AO5*AO4*AO3*AO7 AO8 AO9 AO2*AN7 AN1 AN3 PS1 PO5 PO9 PO10"},
    {FMN," N1  O2  O4  N5  O2* O3* O4* O5* OP1 OP2 OP3 OXT"},
    {Hem," O1D O2D O1A O2A"},
    {Mtx," N1  N3  N5  O   O1  O2  OE1 OE2"},
    {NAD,"AO1 AO2 AO5*AO4*AO3*AO2*AN1 AN3 AN7 NO1 NO2 NO5*NO4*NO3*NO2*NO7 "}
};


/* {0, 0, 2, 2, 2, 0, 0, 2, 0, 2, 0, 2, 0, 0, 0, 2, 1, 0, 0, 0, 0, 2}, /*c*/
/* {0, 0, 2, 2, 2, 0, 0, 2, 0, 2, 0, 2, 0, 0, 0, 2, 0, 0, 0, 2, 0, 2, 0, 2}, /*a*/
/* {0, 0, 2, 2, 2, 0, 0, 2, 0, 2, 0, 2, 0, 0, 0, 2, 2, 0, 2, 0, 0, 2}, /*u*/
/* {0, 0, 2, 2, 2, 0, 0, 2, 0, 2, 0, 2, 0, 0, 0, 2, 0, 0, 2, 0, 0, 0, 1, 0, 2}, /*g*/
/* {0, 0, 2, 2. 2, 0, 0, 2, 0, 2, 0, 2, 0, 0, 0, 2, 0, 0, 2, 0, 0, 0, 2}, /*t*/
/* {0, 0, 2, 2, 2, 0, 2, 2, 2, 0, 2, 2, 2, 0, 0, 2, 0, 2, 0, 2, /*atp*/
/*  0, 0, 0, 2, 0, 0, 0, 2, 0, 2, 0, 2},*/

/* Number of permissible bonds for hydrogen bond donors */
short donors[MAXNAA][TOTNATM] =
{
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*cys*/
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*gly*/
 {1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*his*/
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*lys*/
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*met*/
 {1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 1, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*ser*/
 {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*thr*/
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {1, 0, 0, 0, 0, 0, 0, 2, 2}, /*GLX*/
 {1, 0, 0, 0, 0, 0}, /*UNK*/
 {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*cyh*/
 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*css*/
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*C*/
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*A*/
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*U*/
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*G*/
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*T*/
/* Amendment. RAL 15 Nov 2007 --> */
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*New C*/
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*New A*/
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*New U*/
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*New G*/
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*New T*/
/* <--- RAL 15 Nov 2007 */
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*ATP*/
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*CoA*/
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*FMN*/
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*Hem*/
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*Mtx*/
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},/*NAD*/
 {0, 0, 0, 0, 0}, /*ACE*/
 {0, 0, 0, 0, 0} /*FOR*/
};
/* A ADDITIONAL TABLE FOR THE REALLY LONG AND DIFFICULT ONES */
struct index donors_archive[ARKAA] =
{
/* Amendment. RAL 3 Mar 2011 --> */
  {__C," O2' O3' N4  N4 "},
  {__A," O2' O3' N6  N6 "},
  {__U," O2' O3' N3 "},
  {__G," O2' O3' N1  N2  N2 "},
  {__T," O2' O3' N3 "},
//  {__C," O2* O3* N1  N4  N4 "},
//  {__A," O2* O3* N6  N6 "},
//  {__U," O2* O3* N3 "},
//  {__G," O2* O3* N1  N2  N2 "},
//  {__T," O2* O3* N3 "},
/* <--- RAL RAL 3 Mar 2011 */

/* Amendment. RAL 15 Nov 2007 --> */
  {_DC," O2' O3' N1  N4  N4 "},
  {_DA," O2' O3' N6  N6 "},
  {_DU," O2' O3' N3 "},
  {_DG," O2' O3' N1  N2  N2 "},
  {_DT," O2' O3' N3 "},
/* <--- RAL 15 Nov 2007 */

  {ATP," O2* O3* N6  N6 "},
  {CoA,"AO2*AN6 AN6 PS1 PN4 PN8 PO10"},
  {FMN," N3  O2* O3* O4*"},
  {Hem,""},
  {Mtx," NA2 NA2 NA4 NA4 N8  N  "},
  {NAD,"AO3*AO2*AN6 AN6 NO3*NO2*NN7 NN7 "}
};
            
/**************************************************************************/
/* 8.5 Globals as global variables */

FILE      *ifp=NULL; /* *.new file ? *//*inpdb only*/
int       bonds[MAXCNRECS][MAXCNCOLS] /* 0.9m */; /*copy of CONECT records
                                       used to generate icon[][]*/
/****************************************************************************/
/* 8.7 Minor useful functions */

/***************************************************************************/

/*****************************************************************************/

void getcoord(float *x,float * y,float * z,float * occ, float *b, char * chain,int * n,int * aacode,char * resnam,char * atmnam,int * atmnum,char * buf)
{
    int       i, itemp, debug=0 /*shuts up the debugging chatter here only*/;
    char tempbuf[80];
    
    /* Extracts the following info from the line of input held in BUF[] :
       x,y,z     - the XYZ coordinates
       chain     - the one letter chain code (if any)
       n         - the amino acid resdiue number
       resnam    - the 3-letter residue name
       atmnam    - the 4-letter atom name
       aacode    - an index into RNAMES for this RESNAM 
                   (0 to STDAA-1 if a standard AA)
                   (STDAA to TOTNAA-1 if a non-standard AA)
                   (= TOTNAA if not recognized)
       atmnum    - the atom id number
    */
    if (debug==2) printf("At line %d",__LINE__);
    if (debug) printf("At start of get-co-ord");
    
    strncpy(atmnam, buf+12, 4);
    if (debug==2) printf("At line %d",__LINE__);
    atmnam[4] = '\0';
/*    if (atmnam[2] == ' ')
        atmnam[3] = ' ';/*Why why why ?*/
    if (debug==2) printf("At line %d",__LINE__);
    
    sscanf(buf+6, "%d", atmnum);                 /* atom id number */
    if (debug==2) printf("At line %d",__LINE__);
    sscanf(buf+17, "%c%c%c", resnam, resnam+1, resnam+2);
    resnam[3] = '\0';
    if (debug==2) printf("At line %d",__LINE__);
    *chain = buf[21];
    if (*chain == ' ') *chain = '-';
    
    if (debug==2) printf("At line %d",__LINE__);
    sscanf(buf+22, "%4d", n);
    if (debug==2) printf("At line %d",__LINE__);
    sscanf(buf+27, "%f%f%f", x, y, z);
    /*new in 2.11, inspiration ---pak--- */
    *b= -99;
    *occ= -99;
    strncpy(tempbuf,buf+54,6);
    tempbuf[6]='\0';
    sscanf(tempbuf,"%6f",occ);
    strncpy(tempbuf,buf+60,6);
    sscanf(tempbuf,"%6f",b);
    /*end --- pak --- */
    if (debug==2) printf("At line %d",__LINE__);
    itemp = TOTNAA;
    for (i = 0; i < TOTNAA; i++)
    {
        if (debug) printf("At line %d",__LINE__);
        
        if (!strncmp(rnames[i], resnam, 3))
        {
            itemp = i;
            if (debug) printf("At line %d",__LINE__);
            break;
        }
    }
    *aacode = itemp;
    if (debug) printf("At line %d",__LINE__);
}

void initialise_arrays(void)
{
    int i,atmtyp;
    char *j,tmpstr[5]="\0\0\0\0\0";
    
    /*set up accepts & donors from the archive*/
    for(i=0;i<ARKAA;i++)
    {
        for(j=donors_archive[i].s;*j;j+=4)
        {
            strncpy(tmpstr,j,4);
            atmtyp=scanatmnam(tmpstr,donors_archive[i].n);
            if (atmtyp < TOTNATM)
                donors[donors_archive[i].n][atmtyp]++;
            else
                printf("BUG: %.3s%.4s in donors_archive not recognised\n",
                       rnames[donors_archive[i].n],tmpstr);
        }
        for(j=accepts_archive[i].s;*j;j+=4)
        {
            strncpy(tmpstr,j,4);
            atmtyp=scanatmnam(tmpstr,accepts_archive[i].n);
            if (atmtyp < TOTNATM)
                accepts[accepts_archive[i].n][atmtyp]++;
            else
                printf("BUG: %3.3s%4.4s in accepts_archive not recognised\n",
                       rnames[accepts_archive[i].n],tmpstr);
        }
    }
}


void check_icon(){
    int i,j;
    char buf1[20],buf2[20];
    float length;
    
    for(i=0;i<natoms;i++)
        for(j=0;j<MAXCON;j++)
            if (j< atom[i].ncon)
                if ( (length=vector_length( to(atom[i].p, atom[ icon[i][j] ].p ) ) ) > 3.5 )
                    printf("WARNING: %s - %s bond is %5.3fA long.\n",atomid(i,buf1),atomid(icon[i][j],buf2),length);
                else
                    ;
            else
                if ( icon[i][j]>-1 )
                    printf("BUG: Discounted connectivity records for %s.\n", atomid(i,buf1));
    return;
}

char simplify_sstruc(char sstruc)
{
    if( sstruc == 'S' || sstruc == 'T' || sstruc == 't' )
        return(' ');
    if( sstruc == 'B' )
        return(' ');
    if( islower(sstruc) )
        return(toupper(sstruc));
    return(sstruc);
}

void set_pdbsstruc(struct pdbsstruc * domain, char * buf)
{
    domain -> type = simplify_sstruc( buf[25] );
    if ( domain -> aggr > ' ' )
        domain -> aggr = buf[28];
}

/****************************************************************************/
/*12. Connectivity Calculation Subroutines*/

/* NBONDS contains the number of bonds contained by each of the STDAA amino
   acids, arranged in RNAMES/RESCODE order.
   Each subsequent array (e.g.alabonds) contains the names of the two atoms
   that make up each of the bonds for that particular residue.
*/

int   nbonds[STDAA]={5,6,8,9,12,4,11,8,9,8,8,8,8,9,11,6,7,7,16,13};

char resbonds[MAXNAA][TOTNATM*9]=
{
" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",

" N   CA : CA  C  : C   O  : C   OXT: CA  CB : CB  SG :",/*cys*/

/* It shouldn't matter that the CYS residues are divided into CSS and CYH,
because this should take place after the connectivity records are done. */

" N   CA : CA  C  : C   O  : C   OXT: CA  CB : CB  CG : CG  OD1: CG  OD2:",

" N   CA : CA  C  : C   O  : C   OXT: CA  CB : CB  CG : CG  CD : CD  OE1: CD  OE2:",

" N   CA : CA  C  : C   O  : C   OXT: CA  CB : CB  CG : CG  CD1: CG  CD2: CD1 CE1: CD2 CE2: CE1 CZ : CE2 CZ :",/*phe*/

" N   CA : CA  C  : C   O  : C   OXT:",/*gly*/

" N   CA : CA  C  : C   O  : C   OXT: CA  CB : CB  CG : CG  ND1: CG  CD2: ND1 CE1: CD2 NE2: CE1 NE2:",

" N   CA : CA  C  : C   O  : C   OXT: CA  CB : CB  CG1: CB  CG2: CG1 CD1:", /*ile*/

" N   CA : CA  C  : C   O  : C   OXT: CA  CB : CB  CG : CG  CD : CD  CE : CE  NZ :",/*lys*/

" N   CA : CA  C  : C   O  : C   OXT: CA  CB : CB  CG : CG  CD1: CG  CD2:",/*leu*/

" N   CA : CA  C  : C   O  : C   OXT: CA  CB : CB  CG : CG  SD : SD  CE :", /*met*/

" N   CA : CA  C  : C   O  : C   OXT: CA  CB : CB  CG : CG  OD1: CG  ND2:",

" N   CA : CA  C  : C   O  : C   OXT: CA  CB : CB  CG : CG  CD : CD  N  :" ,

" N   CA : CA  C  : C   O  : C   OXT: CA  CB : CB  CG : CG  CD : CD  OE1: CD  NE2:" ,/*gln*/

" N   CA : CA  C  : C   O  : C   OXT: CA  CB : CB  CG : CG  CD : CD  NE : NE  CZ : CZ  NH1: CZ  NH2:" ,/*arg*/

" N   CA : CA  C  : C   O  : C   OXT: CA  CB : CB  OG :" ,/*ser*/

" N   CA : CA  C  : C   O  : C   OXT: CA  CB : CB  OG1: CB  CG2:" ,/*thr*/

" N   CA : CA  C  : C   O  : C   OXT: CA  CB : CB  CG1: CB  CG2:" ,/*val*/

" N   CA : CA  C  : C   O  : C   OXT: CA  CB : CB  CG : CG  CD1: CG  CD2: CD1 NE1: NE1 CE2\
: CD2 CE3: CD2 CE2: CE3 CZ3: CE2 CZ2: CZ2 CH2: CZ3 CH2:" ,/*trp*/

" N   CA : CA  C  : C   O  : C   OXT: CA  CB : CB  CG : CG  CD1: CG  CD2: CD1 CE1: CD2 CE2\
: CE1 CZ : CE2 CZ : CZ  OH :",  /*tyr*/

" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",
" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",
" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",
" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",
" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",
" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",
" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",
" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",
" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",
" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",/*29*/

" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",
" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",
" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",
" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",
" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",
" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",
" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",
" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",
" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",
" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",/*39*/

" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",
" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",
" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",
" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",
" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",
" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",
" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",
" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",
" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",
" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",/*49*/

" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",
" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",
" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",
" N   CA : CA  C  : C   O  : C   OXT: CA  CB :",
" N   CA : CA  C  : C   O  : C   OXT: CA  CB : CB  SG :",/*cyh*/
" N   CA : CA  C  : C   O  : C   OXT: CA  CB : CB  SG :",/*css*/

/*note the nucleotide backbone is included separately*/
/* Amendment. RAL 3 Mar 2011 --> */
" N1  C1': N1  C2 : C2  O2 : C2  N3 : N3  C4 : C4  N4 : C4  C5 : C5  C6 : C6  N1 :", /*__C*/
" N9  C1': N9  C8 : C8  N7 : N7  C5 : C5  C6 : C6  N6 : C6  N1 : N1  C2 : C2  N3 : N3  C4 \
: C4  C5 : C4  N9 :",/*__A*/
" N1  C1': N1  C2 : C2  O2 : C2  N3 : N3  C4 : C4  O4 : C4  C5 : C5  C6 : C6  N1 :",/*__U*/
" N9  C1': N9  C8 : C8  N7 : N7  C5 : C5  C6 : C6  O6 : C6  N1 : N1  C2 : C2  N2 : C2  N3 \
: N3  C4 : C4  C5 : C4  N9 :",/*__G*/
" N1  C1': N1  C2 : C2  O2 : C2  N3 : N3  C4 : C4  O4 : C4  C5 : C5  C5M: C5  C6 : C6  N1 :",/*__T*/
//" N1  C1*: N1  C2 : C2  O2 : C2  N3 : N3  C4 : C4  N4 : C4  C5 : C5  C6 : C6  N1 :", /*__C*/
//" N9  C1*: N9  C8 : C8  N7 : N7  C5 : C5  C6 : C6  N6 : C6  N1 : N1  C2 : C2  N3 : N3  C4 \
//: C4  C5 : C4  N9 :",/*__A*/
//" N1  C1*: N1  C2 : C2  O2 : C2  N3 : N3  C4 : C4  O4 : C4  C5 : C5  C6 : C6  N1 :",/*__U*/
//" N9  C1*: N9  C8 : C8  N7 : N7  C5 : C5  C6 : C6  O6 : C6  N1 : N1  C2 : C2  N2 : C2  N3 \
//: N3  C4 : C4  C5 : C4  N9 :",/*__G*/
//" N1  C1*: N1  C2 : C2  O2 : C2  N3 : N3  C4 : C4  O4 : C4  C5 : C5  C5M: C5  C6 : C6  N1 :",/*__T*/
/* <--- RAL RAL 3 Mar 2011 */
/* Amendment. RAL 15 Nov 2007 --> */
" N1  C1': N1  C2 : C2  O2 : C2  N3 : N3  C4 : C4  N4 : C4  C5 : C5  C6 : C6  N1 :", /*_DC*/
" N9  C1': N9  C8 : C8  N7 : N7  C5 : C5  C6 : C6  N6 : C6  N1 : N1  C2 : C2  N3 : N3  C4 \
: C4  C5 : C4  N9 :",/*_DA*/
" N1  C1': N1  C2 : C2  O2 : C2  N3 : N3  C4 : C4  O4 : C4  C5 : C5  C6 : C6  N1 :",/*_DU*/
" N9  C1': N9  C8 : C8  N7 : N7  C5 : C5  C6 : C6  O6 : C6  N1 : N1  C2 : C2  N2 : C2  N3 \
: N3  C4 : C4  C5 : C4  N9 :",/*_DG*/
" N1  C1': N1  C2 : C2  O2 : C2  N3 : N3  C4 : C4  O4 : C4  C5 : C5  C7 : C5  C6 : C6  N1 :",/*_DT*/
/* <--- RAL 15 Nov 2007 */

/* other bonds
ATP*/

" PG  O1G: PG  O2G: PG  O3G: PG  O3B: O3B PB : O1B PB : O2B PB \
: O3A PB : PA  O3A: O1A PA : O2A PA : PA  O5*: O5* C5*: C5* C4*\
: C4* C3*: C3* O3*: C3* C2*: C2* O2*: C2* C1*: C1* O4*: O4* C4*\
: C1* N9 : N9  C8 : C8  N7 : N7  C5 : C5  C4 : C4  N9 : C5  C6 \
: C6  N6 : C6  N1 : N1  C2 : C2  N3 : N3  C4 :",

/*COA*/

"PS1 PC2 :PC2 PC3 :PC3 PN4 :PN4 PC5 :PC5 PO5 :PC5 PC6 :PC6 PC7 \
:PC7 PC8 :PC8 PC9 :PC9 PO9 :PC9 PC10:PC10PC11:PC10PO10:PC11PC12\
:PC11PC13:PC11PC14\
:PC12AO6 \
:AO6 AP2 :AP2 AO4 :AP2 AO5 :AP2 AO3 :AO3 AP1 :AP1 AO1 :AP1 AO2 \
:AP1 AP5*\
:AP5*AC5*:AC5*AC4*:AC4*AC3*:AC3*AO3*:AO3*AP3*:AP3*AO8 :AP3*AO9 \
:AP3*AO7 :AC3*AC2*:AC2* O2*:AC2*AC1*:AC1*AO4*:AO4*AC4*\
:AC1*AN9 \
:AN9 AC8 :AC8 AN7 :AN7 AC5 :AC5 AC4 :AC4 AN9 :AC4 AN3 : AN3 AC2\
:AC2 AN1 :AN1 AC6 :AC6 AN6 :AC6 AC5 :",

/*FMN*/

" O2  C2 : C2  N3 : N3  C4 : C4  O4 : C4  C4A: C4A C10: C10 N1 \
: N1  C2 : C10 N10: N10 C9A: C9A C5A: C5A N5 : N5  C4A: C9A C9 \
: C9  C8 : C8  C7 : C7  C6 : C6  C5A: C8  C8M: C7  C7M\
: N10 C1*\
: C1* C2*: C2* O2*: C2* C3*: C3* O3*: C3* C4*: C4* O4*: C4* C5*\
: C5* O5*: O5* P  : P   OP1: P   OP2: P   OP3:",

/*HEM*/
" C1A C2A: C2A C3A: C3A C4A: C4A N A: N A C1A: N AFE  : C1A CHA\
: C3A CMA: C2A CAA: CAA CBA: CBA CGA: CGA O1A: CGA O2A\
: C4A CHB\
: C1D C2D: C2D C3D: C3D C4D: C4D N D: N D C1D: N DFE  : C1D CHD\
: C2D CMD: C3D CAD: CAD CBD: CBD CGD: CGD O1D: CGD O2D\
: C4D CHA\
: C1C C2C: C2C C3C: C3C C4C: C4C N C: N C C1C: N CFE  : C1C CHC\
: C2C CMC: C3C CAC: CAC CBC\
: C4C CHD\
: C1B C2B: C2B C3B: C3B C4B: C4B N B: N B C1B: N BFE  : C1B CHB\
: C2B CMB: C3B CAB: CAB CBB\
: C4B CHC:",

/*MTX*/
" N1  C2 : C2  N3 : N3  C4 : C4  C4A: C4A C8A: C8A N1 : C4A N5 \
: N5  C6 : C6  C7 : C7  N8 : N8  C8A: C2  NA2: C4  NA4: C2  NA2\
: C4  NA4: C6  C9 : C9  N10: N10 CM : N10 C11\
: C11 C12: C12 C13: C13 C14: C14 C15: C15 C16: C16 C11\
: C14 C  : C   N  : C   O  : N   CA : CA  CT : CA  CB : CT  O1 \
: CT  O2 : CB  CG : CG  CD : CD  OE1: CD  OE2:",

/*NAD*/
"AP  AO1 :AP  AO2 :AP   O3 :AP  AO5*:AO5*AC5*:AC5*AC4*:AC4*AO4*\
:AO4*AC1*:AC1*AC2*:AC2*AC3*:AC3*AC4*:AC3*AO3*:AC2*AO2*:AC1*AN9 \
:AN9 AC8 :AC8 AN7 :AN7 AC5 :AC5 AC4 :AC4 AN9 :AC5 AC6 :AC6 AN1 \
:AN1 AC2 :AC2 AN3 :AN3 AC4 :AC6 AN6 \
:NP   O3 \
:NP  NO1 :NP  NO2 :NP  NO5*:NO5*NC5*:NC5*NC4*:NC4*NC3*:NC3*NC2*\
:NC2*NC1*:NC1*NO4*:NO4*NC5*:NC1*NN1 :NN1 NC2 :NC2 NC3 :NC3 NC4 \
:NC4 NC5 :NC5 NC6 :NC6 NN1 :NC3 NC7 :NC7 NO7 :NC7 NN7 :",

" CH3 C  : C  O   : C   OXT:", /*ACE*/
" C   O  : C  OXT :" /*FOR*/

};



/*******************************************************************************/

int  alreadybonded(int atom1,int atom2)


/* returns 1(true) if atoms 1+2 are already bonded (according to ICON), else 0*/

{
    int  i;

    for (i=0; i<MAXCON; i++)
    {
        if (icon[atom1][i]==atom2 || icon[atom2][i]==atom1)
            {
                return(1);
            }
        
    }
    return(0);
}

/******************************************************************************/

int nearlybonded(int atom1,int atom2)
/* returns 1 (true) if atoms 1+2 are 1-3 or 1-4 atoms (according to ICON), or
   else 0
*/

{
     int  i, j, itemp;
     int k,l; /* used when 1-4 bonds are nearlybonded*/
     itemp = debug;
     
/* check for 1-3 pairs first */

     for (i=0; i<atom[atom1].ncon; i++)
     {
         j = icon[atom1][i];
         if (j == -1)
             continue;
         if (alreadybonded(j,atom2))
         {
             if (debug == 1)
             {
                 fprintf(dbgfp, "Atoms %5d and %5d are 1-3 neighbours\n",
                         atom1+1, atom2+1);
                 printf("Atoms %d (%d %s) and %d (%d %s) are 1-3 neighbours through %d (%d %s)\n",
                        atom1,atom[atom1].aanum,atom[atom1].atmnam,
                        atom2,atom[atom2].aanum,atom[atom2].atmnam,
                        j,atom[j].aanum,atom[j].atmnam);
                 
             }

             debug = itemp;
             return(1);
         } /* end if */
     } /* end i loop */
/* this is a line put in by IM to stop nearlybonded throwing out 1-4 pairs */
/*    return(0); */
     
     /* now check for 1-4 pairs */
     for (i=0; i<atom[atom1].ncon; i++)
     {
         k = icon[atom1][i];
         if (k == -1)
             continue;
         for (j=0; j<atom[atom2].ncon; j++)
         {
             l = icon[atom2][j];
             if (l == -1)
                 continue;
             if (alreadybonded(k,l))
             {
                 if (debug == 1)
                 {
                     fprintf(dbgfp, "Atoms %5d and %5d are 1-4 neighbours\n",
                             atom1+1, atom2+1);
                     printf("Atoms %d (%d %s) and %d (%d %s) are 1-4 neighbours through %d (%d %s) and %d (%d %s)\n",
                        atom1,atom[atom1].aanum,atom[atom1].atmnam,
                        atom2,atom[atom2].aanum,atom[atom2].atmnam,
                        k,atom[k].aanum,atom[k].atmnam,
                        l,atom[l].aanum,atom[l].atmnam);
                 }
                 debug = itemp;
                 return(1);
             } /* end if */
         }  /* end of j loop */
     } /* end of i loop */
     debug = itemp;
     return(0);
     
 } /* end of procedure */

/******************************************************************************/
int bondedwithin(int atom1, int atom2, int nbonds)
/* returns 1 (true) if atoms 1+2 are 1-(nbonds+1) atoms are nearer
   (according to ICON) */
{
    int i, j, debug=0;
    if (nbonds<2)
        return( alreadybonded(atom1,atom2));
    else
    {
        for (i=0; i<atom[atom1].ncon; i++)
        {
            j = icon[atom1][i];
            if (j== -1)
                continue;
            if (bondedwithin(j,atom2,nbonds-1))
            {
                if (debug == 1)
                {
                    fprintf(dbgfp, "Atoms %5d and %5d are 1-3 neighbours\n",
                            atom1+1, atom2+1);
                    printf("Atoms %d (%d %s) and %d (%d %s) are 1-3 neighbours through %d (%d %s)\n",
                           atom1,atom[atom1].aanum,atom[atom1].atmnam,
                           atom2,atom[atom2].aanum,atom[atom2].atmnam,
                           j,atom[j].aanum,atom[j].atmnam);
                    
                }
                
                return(1);
            }/* end if bondedwithin*/
        }/*end i loop*/
        return(0);
    }/*end of routine*/
    /*printf("BUG: END OF BONDEDWITHIN REACHED\n");*/
    /*This was removed because it showed up on compilers and worried
      people :) v3.06 */
}

    
        
/******************************************************************************/

int  atomindex(int atomid)
{
    int  i;

    for (i=0; i<natoms; i++)
        if (atom[i].atmnum == atomid)
            return(i);
    return(-1);
}

/******************************************************************************/

void printicon(char * string) /* a debugging tool, I think */
{
    int i, j, k;
    char buf[20];
    
    fprintf(dbgfp, "\n%s\n", string);
    for (i=0; i<natoms; i++)
    {
        fprintf(dbgfp, "%s:", atomid(i,buf));
        for (j=0; j<MAXCON; j++)
        {
            k = icon[i][j];
            if (k < 0)
                fprintf(dbgfp, "             :");
            else
                fprintf(dbgfp, "%s:", atomid(k,buf)+5);
        }
        fprintf(dbgfp, "\n");
    }
}

/******************************************************************************/

void printbonds(void) /* debugging tool */
{
    int i, j;
        
    fprintf(dbgfp, "\nContents of BONDS array:\n");
    for (i=0; i<nconrecs; i++)
    {
        fprintf(dbgfp, "  CONECT %5d : ", bonds[i][0]);
        for (j=1; j<MAXCNCOLS; j++)
        {
            fprintf(dbgfp, "%5d", bonds[i][j]);
        }
        fprintf(dbgfp, "\n");
    }
}

/******************************************************************************/

int readconect(char * fname) /*formerly readpdb*/
/* file fname CONECT records => bonds array, which is a simple copy of CONECT */
{
    FILE  *pdbfp=NULL;
    char  line[160], keywd[7], stemp[6];
    int   i, iatom, offset, numfnd, itemp, token, dbgtemp;
/*note jatom[10] and pdbfn[128] not used.  I have deleted them * -IM */

/* Bug-fix. RAL 1 May 1997 --> */                
/* Bug-fix. RAL 4 Jul 2012 --> */
//    char     keystring[8], message[30];
    char     keystring[8], message[MESSAGE_LEN];
/* <-- RAL 4 Jul 2012 */
    int highest_atom_number;
    int still_reading;

    highest_atom_number = 0;
    still_reading = TRUE;
/* <-- RAL 1 May 1997 */                

    nconrecs = 0;
    dbgtemp  = debug;
    debug = 0;
    pdbfp = fopen(fname,"r");
    if (!pdbfp)
    {
        printf("NOTE: Failed to open file %s for CONECTs\n", fname);
        return (0);
    }
    else
        printf("Reading PDB file \"%s\" for CONECTs . . .\n",fname);
    
    while (!feof(pdbfp))
    {
        if (!fgets(line, 160, pdbfp))
            break;
        
/* Bug-fix. RAL 1 May 1997 --> */                
/*        sscanf (line, "%s", keywd); */
        strncpy(keystring,line,6);
        keystring[6] = ' ';
        keystring[7] = '\0';
        sscanf(keystring, "%s", keywd);
/* <-- RAL 1 May 1997 */                
        if (!keywd[0])
            break;
        token = 0;
        for (i=1; i<=NTOKENS; i++)
            if (!strcmp(keywd, tokstr[i-1]))
                token = i;
        switch (token)
        {

/* Bug-fix. RAL 1 May 1997 --> */                
        case ATOM:
        case HETATM:

          /* If this is an ATOM or HETATM record, and haven't yet hit
             an ENDMDL record, then store the current atom-number as
             the highest encountered so far */
          if (still_reading == TRUE)
            {
              strncpy(keystring,line+6,6);
              keystring[6] = '\0';
              highest_atom_number = atoi(keystring);
            }
          break;

        case ENDMDL:

          /* If this is an ENDMDL record, aren't interested in the
             remainded of the structure */
          if (still_reading == TRUE)
            printf("ENDMDL record encountered after atom number %d\n",
                   highest_atom_number);
          still_reading = FALSE;
          break;
/* <-- RAL 1 May 1997 */                
        case CONECT:
            
            if (nconrecs>=MAXCNRECS)
              {
                sprintf(message,"Too many CONECT records: %d\0",nconrecs);
                fail(message);
              }
            for (i=0; i<MAXCNCOLS; i++)
                bonds[nconrecs][i] = 0;
            offset = 6;
            strncpy (stemp, line+offset, 5);
            stemp[5] ='\0';
            numfnd = sscanf (stemp, "%d", &iatom);
            if (numfnd != 1)
            {
                printf("WARNING: Cannot read first common of CONECT record below\n");
                printf("%s\n",line);
                break;
            }
            
/* Bug-fix. RAL 1 May 1997 --> */                
          /* If the first CONECT atom points to an atom higher than the
             highest read in, then not interested in this CONECT record */
          if (iatom <= highest_atom_number)
            {
/* <-- RAL 1 May 1997 */                
              bonds[nconrecs][0] = iatom;
              for (i=1; i<MAXCNCOLS; i++)
                {
                  offset += 5;
                  strncpy (stemp, line+offset, 5);
                  stemp[5] = '\0';
                numfnd = sscanf(stemp, "%d", &itemp);
                  if (numfnd == 1)
                    bonds[nconrecs][i] = itemp;
                }
              nconrecs++;
/* Bug-fix. RAL 1 May 1997 --> */                
            }
/* <-- RAL 1 May 1997 */                
            break;
        default:
            break;
        }          /* end of switch statement */
    }              /* end of while statement */
    fclose (pdbfp);
    printf ("PDB file contained %d CONECT records \n", nconrecs);
    if (debug == 1)
        printbonds();
    debug = dbgtemp;
    return(1);
    
}                  /* end of procedure */
/*****************************************************************************/

void load_ststan(void) /*set up startup and stopat arrays*/
{
    int     i,kount, stres, curres;
    char    stchn, curchn, stinscode, curinscode;
    int     debug=0;
    
    /* Load start and stop atom numbers for each residue */
    
    


    stchn     = atom[0].chnid;
    stres     = atom[0].aanum;
    stinscode = atom[0].inscode;
    kount     = 0;
    residue[kount].strtat = 0; /* assume the first atom belongs to the first residue!*/

    for (i=0; i < natoms; i++)
    {
        curchn     = atom[i].chnid;
        curres     = atom[i].aanum;
        curinscode = atom[i].inscode;
        if (curchn     != stchn ||
            curres     != stres ||
            curinscode != stinscode)
        {
            kount++;
/*            printf("NOTE: Found %c%4d%c%3s\n", atom[i].chnid, atom[i].aanum, atom[i].inscode, atom[i].resnam );*/
            
            residue[kount-1].stopat = i-1;
            if (debug == 1)
            {
            fprintf(dbgfp,
                   "Residue %4d starts at atom %5d and stops at atom %5d \n",
                    reskount+1, residue[reskount].strtat+1, residue[reskount].stopat+1);
            }

            residue[kount].strtat = i;
            stchn     = curchn;
            stres     = curres;
            stinscode = curinscode;
        }                              /* end if */
    }                                  /* end of i loop */
    residue[kount].stopat = natoms-1;
    if (debug == 1)
    {
        fprintf(dbgfp,"Residue %4d starts at atom %5d and stops at atom %5d \n",
                       reskount+1, residue[reskount].strtat+1, residue[reskount].stopat+1);
    } /* end if */
    if ( kount+1 != reskount )
        printf("BUG: load_ststan finds only %d residues out of %d.\n",kount+1,reskount);
    
} /* end of procedure load_ststan */

/******************************************************************************/
/* Procedure to load into the icon[][] array all the bonds belonging to
   residue ires (numbered from 0 to (reskount-1)). The residue is of the type
  "restyp", numbered from 0 to (STDAA -1). Only applies to standard amino
   acids! */

int apply_bond_to_residue(char * tempstr, int ires)
{
    int iatom,jatom, nconi, nconj;
    char *p1, *p2;
    if (tempstr[8])
        printf("BUG: bond string %s wrong length in apply_bond.\n",tempstr); 

    for (iatom = residue[ires].strtat; iatom<residue[ires].stopat; iatom++)
    {
        p1 = instr(tempstr, atom[iatom].atmnam);
        if (p1== NULL)
            continue;  /* atom "iatom" is NOT involved in bond "ibond": */
        for (jatom = iatom + 1; jatom<=residue[ires].stopat; jatom++)
        {
            char atomidbuf[20];
            
            p2 = instr(tempstr, atom[jatom].atmnam);
            if (p2 == NULL)
                continue; /* atom "jatom" is NOT involved in bond "ibond" */
            if (p1 == p2)
                continue; /* they are the same atom - assumedly with different
                             alternate location indicators */
            
            nconi = atom[iatom].ncon;
            nconj = atom[jatom].ncon;
            if (nconi >= MAXCON)
            {
                printf("WARNING: Excess bonds to %s are ignored.\n", atomid(iatom,atomidbuf));
                break; /* exit jatom loop, try next iatom */
            }
            if (nconj >= MAXCON)
            {
                printf("WARNING: Excess bonds to %s are ignored.\n", atomid(jatom,atomidbuf));
                continue; /* try next jatom */
            }
            if (alreadybonded(iatom,jatom))
                continue;
            icon[iatom][nconi] = jatom;
            icon[jatom][nconj] = iatom;
            atom[iatom].ncon++;
            atom[jatom].ncon++;
            if (debug == 1)
                fprintf(dbgfp, " found bond between atoms %5d and %5d\n",
                        iatom+1, jatom+1);
            return(1);
            
            
        }  /* end of jatom loop */
    }  /* end of iatom loop */
    return (0);
    
}


void load_resbonds(int restyp,int ires)
{
    int  itemp ;
    char tempstr[9];
    char * p;
/* Amendment. RAL 25 May 1997 --> */
/*    char bondlist[1024]="\0"; */
    char bondlist[2048]="\0";    
/* <-- RAL 25 May 1997 */
    int debug=0;

    strcat(bondlist,resbonds[restyp]);
    
    switch (restyp)
    {
    case __C:
    case __A:
    case __U:
    case __G:
    case __T:
/* Amendment. RAL 15 Nov 2007 --> */
    case _DC:
    case _DA:
    case _DU:
    case _DG:
    case _DT:
/* <--- RAL 15 Nov 2007 */
/* Amendment. RAL 3 Mar 2011 --> */
        strcat(bondlist, " P   OP1: P   OP2: P   OP3: P   O5': O5' C5': C5' C4': C4' O4': C4' C3': C3' O3': C3' C2': C2' O2': C2' C1': C1' O4':"); /*phosphoribose nucleotide backbone*/
//        strcat(bondlist, " P   O1P: P   O2P: P   O3P: P   O5*: O5* C5*: C5* C4*: C4* O4*: C4* C3*: C3* O3*: C3* C2*: C2* O2*: C2* C1*: C1* O4*:"); /*phosphoribose nucleotide backbone*/
/* <--- RAL RAL 3 Mar 2011 */
        break;
    default:
        break;
    }
    
        
    if (debug == 1)
    {
        itemp= residue[ires].strtat;
        fprintf(dbgfp, "Loading residue bonds for residue %5d (%3s)\n",
                ires+1, atom[itemp].resnam);
    }
    
    /*    kount =nbonds[restyp];
    for (ibond=0; ibond<kount; ibond++)
    {
        switch (restyp)
        {

        default:
            printf("\nBUG: Unforseen restyp for residue %d in load_resbonds.\n", ires);
            printf("Please mail mcdonald@uk.ac.ucl.biochemistry.bsm\n");
            break;
            
        }   /* end of switch block */
            
    for(p=bondlist;*p;p+=9)
    {
        strncpy(tempstr,p,9);

        tempstr[8] = '\0' ; /* this may be superfluous, but does no harm*/
        if (!apply_bond_to_residue(tempstr,ires) )
            /*printf("WARNING: Covalent Bond %s could not be built.\n",tempstr)*/;
    }  /* end of ibond loop */
    
}  /* end of procedure */

/****************************************************************************/
/*13. Connectivity Calculation Routine*/
/*****************************************************************************/
/* Initializes and loads the ICON[][] array.  Firstly it loads in all bonds
specified by CONECT records as stored in the BONDS[] array.  Then it loads in
all intra-residue bonds that can be determined from a knowledge of the relevant
atom names and residue type.  Lastly it adds in inter-residue (peptide) bonds
based on a distance cutoff criterion between the C atom of one residue and the
N atom of the next. */

void load_icon(void)
{
    int    i, j, bondatom1, bondatom2, iatom1, iatom2, ncon1, ncon2,
           ires, restyp, jres, ityp, jtyp, itemp, jtemp, katom;
    float  dist;
    char   atomidbuf[128];
    
    int debug   = 0;
    
    for (i=0; i<natoms; i++)
    {
        atom[i].ncon = 0;
        for (j=0; j<MAXCON; j++)
            icon[i][j] = -1;

    }


    /* Begin by loading the contents of the BONDS array into ICON */
    printf("%d CONRECS used.\n", nconrecs);
    
    
    if (debug == 1)
        printicon("  ICON array before phase 1:");
    for (i=0; i<nconrecs; i++)
    {
            
        bondatom1 = bonds[i][0];
        if (bondatom1 == 0)
            continue;   /* increment i - try next row of BONDS */
        iatom1 = atomindex(bondatom1);
        if (iatom1 < 0)
        {
            printf("PDBFILE ERROR: Could not find atom with id number %5d in ", bondatom1);
            printf("atom[].atmnum array.\n");
            printf("PDBFILE ERROR: atom %5d was defined by CONECT record %4d \n",
                    bondatom1, i+1);
            continue;   /* increment i - try next row of BONDS */
        }

         
        
        for (j=1; j<MAXCNCOLS; j++)
        {
            
            bondatom2 = bonds[i][j];
            if (bondatom2 == 0)
                continue;              /* next j */
            iatom2 = atomindex(bondatom2);
            
            
            if (iatom2 < 0)
            {
                printf("Could not find atom with id number %5d in ", bondatom2);
                printf("atom[].atmnum array.\n");
                printf("  atom %5d was defined by CONECT record %4d\n",
                        bondatom2, i+1);
                continue;   /* this must be continue not break */
            }
            
            
            if (iatom2 == iatom1)
            {
                printf("CONECT record %4d: ATOM2=ATOM1! \n", i+1);
                continue;     /* next j */
            }
            
            
            if (IGNLB == 1)
            {
                if (atom[iatom1].hetflg && (atom[iatom1].aacode == TOTNAA)&&
                                           (atom[iatom2].aacode < TOTNAA))
                {
                /*  printf("Ignoring CONECT bond from record %d", i+1);
                    printf(" - atom ids are %d and %d\n", bondatom1, bondatom2); */
                    continue; /* try next column j */
                }
                if (atom[iatom2].hetflg && (atom[iatom2].aacode == TOTNAA) &&
                                           (atom[iatom1].aacode < TOTNAA))
                {
                /*  printf("Ignoring CONECT bond from record %d", i+1);
                    printf(" - atom ids are %d and %d\n", bondatom2, bondatom1); */
                    continue;  /* try next column j */
                }
            }
            
            
            if (atom[iatom1].ncon >= MAXCON)
            {
                printf("Warning: too many connectivities for atom id %5d\n",
                        bondatom1);
                break;    /* goto end j loop, next conect record */
            }
            
            
            if (atom[iatom2].ncon >= MAXCON)
            {
                printf("Warning: too many connectivities for atom id %5d\n",
                        bondatom2);
                continue;   /* try next j column */
            }
            if (alreadybonded(iatom1,iatom2))
                continue;
            
            
            ncon1 = atom[iatom1].ncon;
            ncon2 = atom[iatom2].ncon;
            if ( vector_length( to (atom[iatom1].p , atom[iatom2].p ) ) > 5.0 )
                printf("BUG: Long Covalent bond in PDB file\n");
            
                
            
            icon[iatom1][ncon1] = iatom2;
            icon[iatom2][ncon2] = iatom1;
            
            atom[iatom1].ncon++;
            atom[iatom2].ncon++;
        }    /* end of j loop */
    }        /* end of i loop */
    if (debug == 1)
        printicon("  ICON array after phase 1:\n");

    /* Explicit connectivities from CONECT records now loaded into ICON[][]. We
       are now ready to load implicit connectivities based on atom and residue
       names! (for the standard residues only at present) */


    for (ires=0; ires<reskount; ires++)
    {
        itemp = residue[ires].strtat;
        restyp = atom[itemp].aacode;
        if (restyp < TOTNAA)
            load_resbonds(restyp, ires);
    }
    if (debug == 1)
        printicon("  ICON array after phase 2:");


    /* Now load into ICON[][] all peptide bonds, from C of 1 residue to N of the
       next provided they are (a) from the same chain, and (b) close enough
       together! */

    for (ires=0; ires < reskount-1; ires++)
    {
        jres = ires + 1;
        itemp = residue[ires].strtat;
        jtemp = residue[jres].strtat;
        ityp = atom[itemp].aacode;
        jtyp = atom[jtemp].aacode;
        iatom1 = -1;
        iatom2 = -1;
        
        if ( ityp == TOTNAA || jtyp == TOTNAA )
            continue;
        
        if ( !strncmp(necatm[ityp]," N  ",4) && !strncmp(necatm[jtyp]," N  ",4) )
        {

            for (i=residue[ires].strtat; i<=residue[ires].stopat; i++)
            {
                if (atom[i].atmtyp == 2)    /* C atom */
                {
                    iatom1 = i;
                    break;
                }                    /* end if */
            }                        /* end of i loop */
            if (iatom1 == -1)
            {
                printf ("Cannot find atom C for residue %5d ", ires+1);
                katom = residue[ires].strtat;
                printf("(%3s %c %04d%c)\n", atom[katom].resnam,
                                           atom[katom].chnid,
                                           atom[katom].aanum,
                                           atom[katom].inscode);
                continue;   /* try next ires */
            }

            for (j=residue[jres].strtat; j<=residue[jres].stopat; j++)
            {
                if (atom[j].atmtyp == 0)  /* N atom */
                {
                    iatom2 = j;
                    break;
                }                   /* end if */
            }                       /* end of j loop */
            if (iatom2 == -1)
            {
                printf ("Cannot find atom N for residue %5d ", jres+1);
                katom = residue[jres].strtat;
                printf("(%3s %c %04d%c)\n", atom[katom].resnam,
                                           atom[katom].chnid,
                                           atom[katom].aanum,
                                           atom[katom].inscode);
                continue;   /* try next ires */
            }
        }
/*        if (debug)
        {
            printf("Residue %d is ",ires);
            if (isnucleotide(rnames[ityp])) printf("not ");
            printf("a nucleotide, ");

            printf("%d is ",jres);
            if (isnucleotide(rnames[jtyp])) printf("not ");
            printf("one.\n");
            printf("(%d, %d)\n",isnucleotide(rnames[ityp]),isnucleotide(rnames[jtyp]));
            
        }*/
        
        if ( isnucleotide(atom[residue[ires].strtat].resnam) && isnucleotide(atom[residue[jres].strtat].resnam) )
        {
            if (debug) printf("Testing %d and %d for nucleotide bond.\n",ires,jres);
            
            for (i=residue[ires].strtat; i<=residue[ires].stopat; i++)
            {
                if (atom[i].atmtyp == 10) /* O3* */
                {
                    iatom1 = i;
                    break;
                }
            }
            if (iatom1== -1)
            {
                printf ("Cannot find atom O3* for residue of %s.\n",atomid(residue[ires].strtat,atomidbuf));
                continue;
            }
            for (j=residue[jres].strtat; j<=residue[jres].stopat; j++)
            {
                if (atom[j].atmtyp == 2) /* P  */
                {
                    iatom2 = j;
                    break;
                }
            }
            if (iatom2== -1)
            {
                printf("Cannot find atom P for residues of %s.\n",atomid(residue[jres].strtat,atomidbuf));
                continue;
            }
        } 
        if (iatom1 > -1 && iatom2 > -1)
        {
            ncon1 = atom[iatom1].ncon;
            ncon2 = atom[iatom2].ncon;
            if (atom[iatom1].chnid != atom[iatom2].chnid)
            {
                if (debug) printf("Nucleotide-Nucleotide Bond Rejected bec different chains.\n");
                
                continue;  /* try next ires */
            }
            dist = SQR(atom[iatom1].p.x - atom[iatom2].p.x) +
                SQR(atom[iatom1].p.y - atom[iatom2].p.y) +
                    SQR(atom[iatom1].p.z - atom[iatom2].p.z);
            
            if (dist > SQR(PEPBND))
            {
                if (debug) printf("N-N backbone bond rejected because of distance.\n");
                
                continue;  /* try next ires */
            }
            if (alreadybonded(iatom1, iatom2))
                continue;  /* try next ires */
            
            if (ncon1 >= MAXCON)
            {
                char buf[20];
                printf("WARNING: too many bonds with %s, peptide bond ignored.\n",atomid(iatom1,buf));
                
                continue;   /* try next ires */
            }
            if (ncon2 >= MAXCON)
            {
                char buf[20];
                printf("WARNING: too many bonds with %s, peptide bond ignored.\n",atomid(iatom2,buf));
                
                continue;   /* try next ires */
            }
            icon[iatom1][ncon1] = iatom2;
            icon[iatom2][ncon2] = iatom1;
            atom[iatom1].ncon++;
            atom[iatom2].ncon++;
        }   /* end if a couple of atoms have been found */
    }   /* end of ires loop */

    if (debug == 1)
        printicon("  ICON array after phase 3:");

}   /* end of procedure */

/******************************************************************************/
void  find_brakes(void)

/* Procedure to identify chain breaks. Any that are found are stored in the
   BRAKE[][] array such that BRAKE[i][0] contains the lower residue number
   and BRAKE[i][1] the higher residue number of the two residues between
   which there is chain break I. There are 3 criteria for what constitutes
   a chain break, under the control of the CHNBRKFLAG constant:
      1 ==> chain breaks defined by CA-CA distance greater than CADBND;
      2 ==> chain breaks defined by  C-N  distance greater than PEPBND;
      3 ==> chain breaks defined by CA-C-N-CA atom(s) missing;
      4 ==> use all criteria!
*/
{
    float  dist;
    int    ires, jres, iatom, jatom, dbgtemp;
/* Amendment. RAL 23 Feb 2004 --> */        
    static int first = TRUE;
/* <-- RAL 23 Feb 2004 */                

    dbgtemp = debug;
    debug   =  0;
    nbrakes = 0;

    for (ires = 0; ires<(reskount-1); ires++)
    {
        jres = ires + 1;
        if (CHNBRKFLAG == 1 || CHNBRKFLAG == 4)
        {
            dist = 999.9;
            if (residue[ires].ca && residue[jres].ca)
            {
                dist = length_squared( to(*(residue[ires].ca), *(residue[jres].ca)));
            }                
            if (dist>SQR(CADBND))
            {
/* Amendment. RAL 24 Jun 1997 --> */                
/*                if (nbrakes >= MAXBRKS)
                    fail("Too many chain breaks."); */
/* Amendment. RAL 23 Feb 2004 --> */                
/*                if (nbrakes >= MAXBRKS)
                  {
                    printf("*** Too many chain breaks:  %d",nbrakes);
                    fail("Too many chain breaks.");
                    } */
                if (nbrakes < MAXBRKS - 1)
                  {
/* <-- RAL 23 Feb 2004 */                
/* <-- RAL 24 Jun 1997 */                
                    brake[nbrakes][0] = ires;
                    brake[nbrakes][1] = jres;
                    nbrakes = nbrakes + 1;
/* Amendment. RAL 23 Feb 2004 --> */                
                  }
                else if (first == TRUE)
                  {
                    printf("*** Warning. Too many chain breaks:  %d",
                           nbrakes);
                    first = FALSE;
                  }
/* <-- RAL 23 Feb 2004 */                
                if (debug == 1)
                {
                    fprintf(dbgfp,"Chnbrk detected between residues ");
                    fprintf(dbgfp, "%5d and %5d ", ires+1, jres+1);
                    iatom = residue[ires].strtat;
                    jatom = residue[jres].strtat;
                    fprintf(dbgfp,"(%3s %c%04d%c and %3s %c%04d%c)\n",
                            atom[iatom].resnam, atom[iatom].chnid,
                            atom[iatom].aanum, atom[iatom].inscode,
                            atom[jatom].resnam, atom[jatom].chnid,
                            atom[jatom].aanum, atom[jatom].inscode);
                }
                continue;  /* look at next residue pair */
            }   /* end if */
        }   /* end if */
        if (CHNBRKFLAG == 1)
            continue;

        if (CHNBRKFLAG == 2 || CHNBRKFLAG == 4)
        {
            dist = 999.9;                        

            if (residue[ires].c && residue[jres].n)
            {
                dist = length_squared( to(*(residue[ires].c),*(residue[jres].n)));
            }
            if (dist > SQR(PEPBND))
            {
/* Amendment. RAL 24 Jun 1997 --> */                
/*                if (nbrakes >= MAXBRKS)
                    fail("Too many chain breaks."); */
/* Amendment. RAL 23 Feb 2004 --> */                
/*                if (nbrakes >= MAXBRKS)
                  {
                    printf("*** Too many chain breaks:  %d",nbrakes);
                    fail("Too many chain breaks.");
                  } */
                if (nbrakes < MAXBRKS - 1)
                  {
/* <-- RAL 23 Feb 2004 */                
/* <-- RAL 24 Jun 1997 */                
                    brake[nbrakes][0] = ires;
                    brake[nbrakes][1] = jres;
                    nbrakes = nbrakes + 1;
/* Amendment. RAL 23 Feb 2004 --> */                
                  }
                else if (first == TRUE)
                  {
                    printf("*** Warning. Too many chain breaks:  %d",
                           nbrakes);
                    first = FALSE;
                  }
/* <-- RAL 23 Feb 2004 */                
                if (debug == 1)
                {
                    fprintf (dbgfp,"Chnbrk detected between residues ");
                    fprintf (dbgfp, "%5d and %5d ", ires+1, jres+1);
                    iatom = residue[ires].strtat;
                    jatom = residue[jres].strtat;
                    fprintf(dbgfp,"(%3s %c%04d%c and %3s %c%04d%c)\n",
                            atom[iatom].resnam, atom[iatom].chnid,
                            atom[iatom].aanum, atom[iatom].inscode,
                            atom[jatom].resnam, atom[jatom].chnid,
                            atom[jatom].aanum, atom[jatom].inscode);
                }
                continue;  /* look at next residue pair */
            }   /* end if */
        }   /* end if */
        if (CHNBRKFLAG == 2)
            continue;

        if (CHNBRKFLAG == 3 || CHNBRKFLAG == 4)
        {
            if (!residue[ires].ca || !residue[ires].c ||
                 !residue[jres].n || !residue[jres].ca)
            {
/* Amendment. RAL 24 Jun 1997 --> */                
/*                if (nbrakes >= MAXBRKS)
                    fail("Too many chain breaks."); */
/* Amendment. RAL 23 Feb 2004 --> */                
/*                if (nbrakes >= MAXBRKS)
                  {
                    printf("*** Too many chain breaks:  %d",nbrakes);
                    fail("Too many chain breaks.");
                  } */
                if (nbrakes < MAXBRKS - 1)
                  {
/* <-- RAL 23 Feb 2004 */                
/* <-- RAL 24 Jun 1997 */                
                    brake[nbrakes][0] = ires;
                    brake[nbrakes][1] = jres;
                    nbrakes = nbrakes + 1;
/* Amendment. RAL 23 Feb 2004 --> */                
                  }
                else if (first == TRUE)
                  {
                    printf("*** Warning. Too many chain breaks:  %d",
                           nbrakes);
                    first = FALSE;
                  }
/* <-- RAL 23 Feb 2004 */                
                if (debug == 1)
                {
                    fprintf(dbgfp,"Chnbrk detected between residues ");
                    fprintf(dbgfp,"%5d and %5d ", ires+1, jres+1);
                    iatom = residue[ires].strtat;
                    jatom = residue[jres].strtat;
                    fprintf(dbgfp,"(%3s %c%04d%c and %3s %c%04d%c)\n",
                            atom[iatom].resnam, atom[iatom].chnid,
                            atom[iatom].aanum, atom[iatom].inscode,
                            atom[jatom].resnam, atom[jatom].chnid,
                            atom[jatom].aanum, atom[jatom].inscode);
                }
                continue; /* next ires-jres pair */
            } /* end if */
        }     /* end if */
    }         /* end of ires loop */
    debug = dbgtemp;
}             /* end procedure */

int cagap(int residue1,int residue2)
{
    int i;
    
    if (!nbrakes)
        return(abs(residue1-residue2));
    if(residue1>residue2)
    {
        i=residue1;
        residue1=residue2;
        residue2=i;
    }
    
    for(i=0;i<nbrakes;i++)
    {
        if(brake[i][0]>=residue1 && brake[i][1]<=residue2)
            return(-1);
    }
    return (abs(residue1-residue2));
}

/***************************************************************************/
/*13.5 Adding new residues to the selection                                */
/***************************************************************************/

int add_residue_type( char * residue, char * similar_residue )
{
    int similar_resnum=TOTNAA, resnum=TOTNAA, i;

    if ( MAXNAA == TOTNAA )
    {
        printf("BUG: \"%s\" is one residue type too many.\n",residue);
        return(0);
    }
    
    if (strlen(residue) != 3)
    {
        printf("ERROR: \"%s\" must be 3 characters in add_residue_type.\n",residue);
        return(0);
    }
    if ( * similar_residue )
    {
        if (strlen(similar_residue) != 3)
        {
            printf("ERROR: \"%s\" must be 3 characters in add_residue_type.\n",similar_residue);
            return(0);
        }
        similar_resnum= scanresnam( similar_residue);
        if (similar_resnum==TOTNAA)
        {
            printf("ERROR: \"%s\" not recognised in add_residue_type.\n",similar_residue);
            return(0);
        }
    }

/* The bugchecking has been passed - 
   we are now dealing with a "real" residue */
    
/* Change rnames, TOTNAA . . . */
    
    strncpy(rnames[TOTNAA],residue,3);
    TOTNAA++;
    
    if ( !*similar_residue )
        return(1);
    
/* Set the resbonds, donors, accepts, and necatm arrays */
    resnum = TOTNAA-1; /*3.02*/
    
    strcpy( necatm[resnum], necatm[similar_resnum]);
    strcpy( resbonds[resnum], resbonds[similar_resnum]);
    for (i=0;i<TOTNATM;i++)
    {
        donors[resnum][i]=donors[similar_resnum][i];
        accepts[resnum][i]=accepts[similar_resnum][i];
    }
    return(1);

}

int add_atoms( char * residue, char * atoms)
{
    int resnum;
    if ( (resnum=scanresnam(residue)) == TOTNAA)
    {
        printf("ERROR: Residue %s not recognised in add_atoms.\n", residue);
        return(0);
    }
    if ( strlen( atoms ) % 4 )
    {
        printf("ERROR: Following atom list for %s does not divide into 4 character atoms.\n", residue);
        printf("ERROR: %s\n", atoms);
        return(0);
    }
    if ( strlen( necatm[resnum] ) + strlen( atoms ) > TOTNATM*4 )
    {
        printf("BUG: Too many atoms in necatm for residue %s.\n", residue);
        return(0);
    }
    else
    {
        strcat( necatm[resnum], atoms );
        return(1);
    }
}

int add_bonds( char * residue, char * new_bonds)
{
    int resnum;
    char *i;
    
    if ( (resnum=scanresnam(residue)) == TOTNAA)
    {
        printf("ERROR: Residue %s not recognised in add_bonds.\n", residue);
        return(0);
    }
    if ( strlen( new_bonds ) % 9 )
    {
        printf("ERROR: Following bonds list for %s does not divide into 9 character bonds.\n", residue);
        printf("ERROR: %s\n", new_bonds);
        return(0);
    }
    for(i=resbonds[resnum];*i;i+=9)
        if ( *(i+8) != ':' )
        {
            printf("ERROR: Colon separator missing for %s in add_bonds.\n",residue);
            printf("ERROR: %s\n", new_bonds);
            return(0);
        }

    if ( strlen( resbonds[resnum] ) + strlen( new_bonds ) > TOTNATM*9 )
    {
        printf("ERROR: Too many bonds in resbonds for residue %s.\n", residue);
        return(0);
    }
    else
    {
        strcat( resbonds[resnum], new_bonds );
        return(1);
    }
}

int add_donacc ( char * residue, char * atoms, short number, short donflg )
{
    int resnum,atmnum;
    char * i;
    
    if ( (resnum=scanresnam(residue))==TOTNAA)
    {
        printf("BUG: Residue %s not recognised in add_donacc.\n", residue);
        return(0);
    }
    if ( strlen( atoms ) % 4 )
    {
        printf("BUG: Following atom list for %s does not divide into 4 character atoms.\n", residue);
        printf("BUG: %s\n", atoms );
        return(0);
    }
    for(i=atoms;*i;i+=4)
    {
        atmnum=scanatmnam(i,resnum);
        if (atmnum==TOTNATM)
            printf("BUG: Atom %.3s%.4s not recognised in add_donacc.\n", residue, i);
        else
            if (donflg)
                donors[resnum][atmnum]=number;
            else
                accepts[resnum][atmnum]=number;
    }
    return(1);
            
}

int add_donor(char * residue, char * atoms, short number)
{
    return(add_donacc(residue,atoms,number,1));
}
    
int add_acceptor( char * residue, char * atoms, short number )
{
    return(add_donacc(residue,atoms,number,0));
}

void supplement_arrays()
{
    add_residue_type("SO4","");
    add_atoms("SO4"," S   O1  O2  O3  O4 ");
    add_bonds("SO4"," S   O1 : S   O2 : S   O3 : S   O4 :");
    add_acceptor("SO4", " O1  O2  O3  O4 ",2);

    add_residue_type("FMD","");
    add_atoms("FMD"," N   C   O  ");
    add_bonds("FMD"," N   C  : C   O  :");
    add_donor("FMD"," N  ",2);
    add_acceptor("FMD"," O  ",2);
    
    add_residue_type("FML","");
    add_atoms("FML"," C   O  ");
    add_bonds("FML"," C   O  :");
    add_acceptor("FML"," O  ",2);
    

/*ADD NEW RESIDUE TYPES BELOW THIS LINE=====================================*/


}


/****************************************************************************/

void find_ss(void)
{
    int i,j;
    float d;
    
    /* Check for disulphide bridges */

    printf("Checking for disulphide bridges . . . \n");
    for (i = 0; i < natoms; i++)
        for (j = i + 1; j < natoms; j++)
            if (atom[i].aacode == Cys && atom[i].atmtyp == 5 && 
                atom[j].aacode == Cys && atom[j].atmtyp == 5 &&
                fabs(atom[i].p.x - atom[j].p.x) < SSDIST && 
                fabs(atom[i].p.y - atom[j].p.y) < SSDIST && 
                fabs(atom[i].p.z - atom[j].p.z) < SSDIST)
            {
                d = SQR(atom[i].p.x - atom[j].p.x) + 
                    SQR(atom[i].p.y - atom[j].p.y) + 
                    SQR(atom[i].p.z - atom[j].p.z);
                if (d <= SQR(SSDIST))
                    {
                        char buf1[20],buf2[20];
                        
                        atom[i].ssflg = atom[j].ssflg = TRUE;
                        atom[i].nh = atom[j].nh = 0;
                        /* this should have hbout from outputting the
                           SH hydrogens in their default positions */
                        ssnum++;
                        printf("Disulphide bond between %s and %s.\n",atomid(i,buf1),atomid(j,buf2));
                    }
            }
    printf("%d disulphide bonds found.\n\n",ssnum);
}

int scanatmnam(char * atmnam,int res)
{
    char *p;
    char tmpstr[5]="\0\0\0\0\0";
    
    strncpy(tmpstr,atmnam,4);
    
    /* return TOTNATM for failure */
    /* otherwise the atomnumber*/

    p = instr(necatm[res],tmpstr);
    if (p!=NULL)
        return ((p-necatm[res])/4);
    else
        return (TOTNATM);
}


int scanresatmnam(char inatmnam[8],int * res, int * atm)
{
/* identify the residue and atom numbers of res/atm id.
   Return "0" for success,
   Return "1" for resid only,
   Return "2" for neither identified.*/
    int i,j,debug=0;
    char tmpatmnam[8]="\0\0\0\0\0\0\0\0";
    strncpy(tmpatmnam,inatmnam,7);
    
    
    for (i = 0; i < TOTNAA; i++)
    {
        if (!strncmp(rnames[i], tmpatmnam, 3))
        {
            break;
        }
    }
    *res = i;
    
    if (i==TOTNAA)
    {
        *res=TOTNAA;
        *atm=TOTNATM;
        if (debug) printf("Residue Unrecognised\n");
        
        return(2);
    }
    
    j=scanatmnam(tmpatmnam+3,i);
        
    if (j==TOTNATM)
    {
        *res=i;
        *atm=TOTNATM;
        if (debug) printf("Atom Unrecognised\n");
        
        return(1);
    }
    else
    {
        *res=i;
        *atm=j;
        return(0);
        /* success !! */
    }
}

int scanresnam(char * inresnam)
{
    /* identify the residue numbers of resid. */

    int i;
    char tmpresnam[4]="\0\0\0\0";
    strncpy(tmpresnam,inresnam,3);
    
    
    for (i = 0; i < TOTNAA; i++)
    {
        if (!strncmp(rnames[i], tmpresnam, 3))
        {
            return(i);
            
        }
    }
    return(i);

}

void print_atoms(char * string) /*debugging only*/
{
    int i;
    char buf[20];
    
    printf("%s\n",string);
    for(i=0;i<natoms;i++)
    {
        printf("%6d %s\n",i,atomid(i,buf));
    }
}



/****************************************************************************/
/*15. inpdb_file */

short inpdb_file(char * fname, char * inpdbfn)
/* returns flag of whether structure usable, takes clean and unclean pdb */
{
    int       i,j, token, namino /*=].aanum*/, aac/*=.aacode*/;
    int       oldresnum, atmnum;
    float     occ, b; /*---*pak*---*/
    
    short     caonly; /*flag to determine if structure usable*/
    char     *p; /*pointer to position in necatm*/
    char      atmnam[5]="\0\0\0\0\0", resnam[4]="\0\0\0\0";
    char chain, strucsum, inscode,
    oldchain, oldinscode, space, altcode, lurn[4]="???\0";
              /*lurn is a store for any unrecognised residues*/
    float     x, y, z;
    char      sstfn[128]="\0", sstbuf[160]="\0";
    char      pdbfn[128]="\0" /* If inpdbfn absent find a pdb file*/;

    char     keyword[7],buf[128];
/* Bug-fix. RAL 1 May 1997 --> */                
    char     keystring[8];

    int      keep_reading;
/* <-- RAL 1 May 1997 */                
    
    FILE     *ifp=NULL, *sstfp=NULL;
    int      hkount = 0/*where to put next atom's h*/;
    int      debug=0;
    int      firstofchain /* flag for whether or not to put NH3 in */;
    int      ownh=0 /*flag for if the record has any of its own Hs - only
                     used so that if such a message is sent to the user,
                     it is only sent once*/;
    int      pdbinflg=0; /* has the pdb file been gleaned for CONECT records ? */
    short int first_ATOM_flg=1; /* so the first sstrucseg starts at the first residue, and any N-terminyl additions don't mess up the sstrucsegs.*/
    
    int last_aares = 0;
    char last_sstruc = 0, sstrucseg_aggr = ' ', last_chnbrk_chr;
/* Bug-fix. RAL 4 Jul 2012 --> */
//    char message[30];
    char message[MESSAGE_LEN];
/* <-- RAL 4 Jul 2012 */

/* Bug-fix. RAL 1 May 1997 --> */                
    keep_reading = TRUE;
/* <-- RAL 1 May 1997 */                


    printf("Processing file \"%s\" . . .\n", fname);
    printf("Opening \"%s\" for protein co-ordinates. . .\n", fname);
    ifp = fopen(fname, "r");
    if (!ifp)
    {
        printf("Failed to open specified input file %s!\n",fname);
        return(0);
    }

    printf("\n");

    natoms     = 0;
    space      = ' ';
    caonly     = TRUE;
    chain      = '?';
    oldchain   = '?';
    oldresnum  = -999;
    oldinscode = '?';
    reskount   = 0;
    aareskount = 0;
    
    for (i=0;i<MAXNATM;i++)
    {
        h_atm[i].p=h_atm[i].a=vector_of(-99.9,-99.9,-99.9);
        h_atm[i].typ=0;
        atom[i].nh=0;
        atom[i].h_ptr= -99; /*I don't know how effective this would be */
        atom[i].occ= -99.9;
        atom[i].b= -99.9;
        atom[i].aanum= -99;
        atom[i].aacode= -99;
        atom[i].atmtyp= -99;
        atom[i].caindex= -99;
        atom[i].atmnum= -99;
        atom[i].acc = -99.9;
        atom[i].ndonhb=0;
        atom[i].nacchb=0;
        atom[i].don_p=NULL;
        atom[i].acc_p=NULL;
        
        /* char, hetflg, ssflg and ncon are not set to -99 just in case . . */
        /* h_atm[i].typ set to 0, because that represents no H */
    }
    
    

    if (debug) printf("Entering the loop\n");
/* Bug-fix. RAL 1 May 1997 --> */                
/*    while (!feof(ifp)) */
    while (!feof(ifp) && keep_reading == TRUE)
/* <-- RAL 1 May 1997 */                
    {
        int debug=0;/*setting debug for within the loop*/
        
/* Amendment. RAL 3 Jul 2012 --> */
//        if (!fgets(buf, 160, ifp))
//            break;
        if (!fgets(buf, LINE_LEN, ifp))
            break;
        string_truncate(buf,LINE_LEN);
/* <-- RAL 3 Jul 2012 */

/* Bug-fix. RAL 1 May 1997 --> */                
/*        sscanf(buf, "%s", keyword);  */                /* Read the record name */
        strncpy(keystring,buf,6);
        keystring[6] = ' ';
        keystring[7] = '\0';
        sscanf(keystring, "%s", keyword);
/* <-- RAL 1 May 1997 */                
        if (debug) printf("Keyword %s", keyword);
        
        if (!keyword[0])           /* Odd - there isn't a record name! Exit. */
            break;

        token = 0;
        for (i = 1; i <= NTOKENS; i++)                 /* Decode record type */
            if (!strcmp(keyword, tokstr[i - 1]))
                token = i;

        if (token == ENDENT)
            break;
        
        if (debug) printf("  Token %d, %s about to be processed\n", token, keyword);
        switch (token)
        {


        case HEADER:
            if (buf[9] != space)
            {
                printf("HEADER continuation line identified and ignored!\n");
                break;
            }
            strncpy(brcode, buf+62, 4);
            brcode[4] = '\0';

            if (*inpdbfn)
            {
                printf("Trying to open named file %s for CONECT records . . .\n",inpdbfn);
                
                pdbinflg = readconect(inpdbfn); /* use any GIVEN unclean file*/
            }
            
            if (!pdbinflg) /*find local unclean pdb file */
            {
                printf("Looking for uncleaned pdb format file for CONECT records . . .");
                

                strcpy(pdbfn, "p");
                strcat(pdbfn, brcode);
                strcat(pdbfn, ".pdb");
                for (p = pdbfn; *p; p++)
                    if (isupper(*p))
                        *p = tolower(*p);
                pdbinflg = readconect(pdbfn);
            }
#ifdef BSM
            if (!pdbinflg)
            {
                strcpy(pdbfn, "/home/bsm/mcdonald/s/data/set/p");
                strcat(pdbfn, brcode);
                strcat(pdbfn, ".pdb");
                for (p = pdbfn; *p; p++)
                    if (isupper(*p))
                        *p = tolower(*p);
                pdbinflg = readconect(pdbfn);
            }

            if (!pdbinflg)
            {
                strcpy(pdbfn, "/home/bsm/mcdonald/s/data/pdb/pdb");
                strcat(pdbfn, brcode);
                strcat(pdbfn, ".ent");
                for (p = pdbfn; *p; p++)
                    if (isupper(*p))
                        *p = tolower(*p);
                pdbinflg = readconect(pdbfn);
            }

            if (!pdbinflg)
            {
                strcpy(pdbfn, "/data/pdb/prerelease/pdb");
                strcat(pdbfn, brcode);
                strcat(pdbfn, ".ent");
                for (p = pdbfn; *p; p++)
                    if (isupper(*p))
                        *p = tolower(*p);
                pdbinflg = readconect(pdbfn);
            }

            if (!pdbinflg)
            {
                strcpy(pdbfn, "/data/pdb_release/jan94/p");
                strcat(pdbfn, brcode);
                strcat(pdbfn, ".pdb");
                for (p = pdbfn; *p; p++)
                    if (isupper(*p))
                        *p = tolower(*p);
                pdbinflg = readconect(pdbfn);
            }

            if (!pdbinflg)
            {
                strcpy(pdbfn, "/data/pdb/p");
                strcat(pdbfn, brcode);
                strcat(pdbfn, ".pdb");
                for (p = pdbfn; *p; p++)
                    if (isupper(*p))
                        *p = tolower(*p);
                pdbinflg = readconect(pdbfn);
            }

#else /*matches ifdef BSM */
            if (!pdbinflg) /* last conect resort is to use only given file */
            {
                printf("Looking at named input file for CONECT records . . .\n");
                pdbinflg = readconect(fname);
            }
#endif            
            if (debug==1) printf("finished readconect function\n");

            
            /* now we've finished conect, look for sst */
            if (inputsstflg) 
            {
                /* start of bit to get *.sst */
#ifdef BSM
                strcpy(sstfn, SSTHEAD );   
                strcat(sstfn, brcode);
                strcat(sstfn, ".sst");
                if (debug==1) printf("About to set sstfn to lower case\n");
                
                for (p=sstfn; *p; p++)
                    if (isupper(*p))
                        *p = tolower(*p);
                if (debug==1) printf ("Just set sstfn to lower case\n");
                printf("Attempting to open standard SST file %s\n", sstfn);
                sstfp = fopen(sstfn, "r");
#endif
                if (!sstfp)
                {
#ifdef BSM
                    printf("Failed to open standard SST file %s\n", sstfn);
#endif
                    strcpy(sstfn, "p");   
                    strcat(sstfn, brcode);
                    strcat(sstfn, ".sst\0");
                    if (debug==1) printf("About to set sstfn to lower case\n");
                    
                    for (p=sstfn; *p; p++)
                        if (isupper(*p))
                            *p = tolower(*p);
                    if (debug==1) printf ("Just set sstfn to lower case\n");
                    sstfp = fopen(sstfn, "r");
                    if (!sstfp)
                    {
                        printf("Failed to open current directory SST file %s\n", sstfn);
                        /* beep point ? */
                        /* before version 1.0w onwards, runs that failed to
                           find and sst file were cancelled. */
                    }
                    
                }
                printf("Opened SST file %s\n", sstfn);
                /* note it would not crash the computer if sstfp = 0 */
            }
            else /* matches if inputsstflg */
                sstfp=NULL;
            
            for (i=0; i<7 && sstfp ; i++)
                fgets(sstbuf,160,sstfp);
            break;

        case ENDMDL:

            /* If this is an ENDMDL record signifying the end of an NMR
               structure, then want to stop reading here */
            keep_reading = FALSE;

            break;
            
        /* Start of the "input a file" routine ----------------------------*/
        case HETATM:
        case ATOM:

            if (debug==1){
                printf("%s\n",buf);                
                debug=2;
            }
            
            if (natoms >= MAXNATM)
              {
                sprintf(message,"Too many atoms (%d)! Increase MAXNATM.\0",
                        natoms);
                fail(message);
              }
            
/* get the co-ords, the inscode, and the altcode */
            altcode = buf[16];
            if (altcode == space)
                altcode = '-';


/* remove alternates was commented out in version 2.0c on the grounds that
we are looking for potential hydrogen bonds, and these includes those which
involve alternative locations.  Of course, if a "new" file is in use, this
should not be a problem. */
/*   if (altcode != space && altcode != 'A') continue; */

            inscode = buf[26];
            if (inscode == space)
                inscode = '-';
                        
            getcoord(&x, &y, &z, &occ, &b, &chain, &namino, &aac, resnam, atmnam, &atmnum,buf);
            if (*only_chainid_lst && !strchr(only_chainid_lst,chain))
                break;
            
            if (aac == TOTNAA && !iswater(resnam)) 
                /* unrecognised ATOM residue */
            {
                if (strncmp(resnam, lurn, 3))
                {
                    printf("WARNING: Residue %3s is not recognized by HBPLUS \n", resnam);
                    strncpy(lurn, resnam, 3);
                }
                
                /*break;*/
                /* ignore unrecognized residues */
                /* this was changed in 2.0c to regard unrecognised residues as HETs */
                /* and again in 2.08 to ignore intatm/hetatm difference in PDB file when dealing with it */                  
            }
            /* Either this atom is a from
               a recognized residue. Check first whether or not it is the
               first record of a new residue. If it is, update the 
               residue counters, and initialize the CAXYZ coordinate array.
               If it is an ATOM residue as well then read in its SSSUM from 
               the SST file. */

            firstofchain= oldchain!=chain;
                /*IM- so that NH3s can be put on multi-chain proteins*/
            if (  chain != oldchain   || 
                 namino != oldresnum  ||
                inscode != oldinscode) /* if newresidue */
            {
                /*v2.29*/       if(oldchain != chain) num_chains++;
                oldchain   = chain;             
                oldresnum  = namino;
                oldinscode = inscode;
                /*debug line*/
/*                printf("NOTE: Found %c%04d%c%3s\n",chain,namino,inscode,resnam);*/
                reskount   ++; /* update counts */
                    
/* Bug-fix. RAL 3 Jul 2012 --> */
//                fprintf(dbgfp, "Residue %5d is %3s %c %4d%c%c\n",
//                                reskount, resnam, chain, namino, inscode,
//                                strucsum);
                fprintf(dbgfp, "Residue %5d is %3s %c %4d%c\n",
                                reskount, resnam, chain, namino, inscode);
/* <-- RAL 3 Jul 2012 */
                if ((reskount-1) >= MAXNRES)
                  {
                    sprintf(message,
                            "Too many residues (%d)! Increase MAXNRES.\0",
                            reskount-1);
                    fail(message);
                  }

                /* set all the arrays to blank */
                residue[reskount-1].ca = NULL;
                residue[reskount-1].c = NULL;
                residue[reskount-1].n = NULL;
                residue[reskount-1].n = NULL;
                
            }

            p = NULL;
            
/* A quick summary of if . . . thens
            if ((aac < TOTNAA)) /* note ATOM=) acc < TOTNAA 
            IF RESIDUE RECOGNISED    
                if ((p = instr(necatm[aac], atmnam)) != NULL)
                    IF ATOM RECOGNISED ACT NORMALLY
                else
                    if (atmnam[1]=='H'||atmnam[1]=='D')/* then it's H 
                        IF HYDROGEN TRY TO PLACE IT
                    elsed
                        GIVE UP ON THAT LINE

            else  /* atom is from a HETATM record of not-identified class
                if (!(strncmp(resnam,"DO",2) && strncmp(resnam,"HO",2)))
                IT IS WATER
                    if (atmnam[1]=='O')
                        IT IS A WATER OXYGEN
                    else
                        ASSUME IT IS HYDROGEN
                
            IF TOTALLY UNRECOGNISED RESIDUE REGISTER IT ANYWAY    
*/
            if ((aac < TOTNAA)) /*because unrecs -> HETATM */
            /* if recognised ATOM or HETATM */
            {
                if (debug)
                    printf("\tResidue type recognised\n");
                
                if ((p = instr(necatm[aac], atmnam)) != NULL)
                /* if a recognised heavy atom */
                {
                    if (debug)
                        printf("\tAtom type recognised\n");
                    
                    atom[natoms].caindex =reskount -1;
                    atom[natoms].atmtyp  = (p-necatm[aac])/4;
                    /* I hope that that use of that array works! */
                    atom[natoms].nh = donors[aac][atom[natoms].atmtyp];
                    if (firstofchain && atom[natoms].atmtyp == 0)
                        atom[natoms].nh = 3 /*NH3 terminus*/;
                    atom[natoms].h_ptr = hkount;
                    hkount += atom[natoms].nh;
                        

                    if ((instr(" N   CA  C   O  ", atmnam)) != NULL)/* if amino acid backbone*/
                    {
                        
                        
                        if (atom[natoms].atmtyp == 0)
                            residue[reskount-1].n= &(atom[natoms].p);
                        
                        
                        if (atom[natoms].atmtyp == 1) /*set up ca*/
                            residue[reskount-1].ca= &(atom[natoms].p);
                        else
                            caonly = FALSE;
                        
                        if (atom[natoms].atmtyp == 2)
                            residue[reskount-1].c = &(atom[natoms].p);
                        
                        if (atom[natoms].atmtyp == 3)
                            residue[reskount-1].o = &(atom[natoms].p);
                    }/* matches 'if it recognised residue and recognised as mc */
                   
                }/* matches 'if it is a recognised heavy atom' */
                else
                /* if it is not a recognised heavy atom */
                {
                    if (debug)
                        printf("\tAtom not recognised\n");
                    
/* Amendment. RAL 11 Dec 2013 --> */
                    // if (atmnam[1]=='H'||atmnam[1]=='D')/* then it's H */
                    if (atmnam[1] == 'H' || atmnam[1] == 'D' ||
                        (atmnam[0] == 'H' && atmnam[3] != ' '))
/* <-- RAL 11 Dec 2013 */
                    /* it is not recognised - and it is H */
                    {
                        int i;
                        int      hvy_atm /* what the H is mounted on */;
                        if (!ownh || debug) {
                            printf("\tHydrogens Recognized in Brookhaven File\n");
                            ownh=1;
                        }
                        
/*                        if (debug) printf ("%d && %d\n",(int)(strncmp(atmnam+1,"H  ",3)),(int)( strncmp(atmnam+1,"D  ",3)) );*/
                        if (strncmp(atmnam+1,"H  ",3) && strncmp(atmnam+1,"D  ",3))
                            hvy_atm=find_atom2(reskount-1,chain,atmnam,natoms-1); /* this would place " H  " on " O  " */
                        else
                            hvy_atm=find_atom4(chain,namino,inscode,resnam," N  ",natoms-1);
                        /* This is to compensate for the way NHs are shown */
                        
                        if (hvy_atm == -99 )
                        {
                            printf("WARNING: Heavy atom for pdb file Hydrogen not found.\n");
                            printf("%s\n",buf);
                        }
                        else
                        {
                            i=0;
                            if (debug) printf("Hydrogen Assigned to /%c%04d%c%.3s%.4s\n",atom[hvy_atm].chnid,atom[hvy_atm].aanum,atom[hvy_atm].inscode,atom[hvy_atm].resnam,atom[hvy_atm].atmnam);
                            /*                            if (debug) printf("hvy_atm %d atom[hvy_atm].nh %d atom[hvy_atm].h_ptr %d\n",hvy_atm,atom[hvy_atm].nh,atom[hvy_atm].h_ptr);*/
                            
                            while (i<atom[hvy_atm].nh && h_atm[atom[hvy_atm].h_ptr+i].typ)
                                
                            {
                                i++;
                                /*printf("Hydrogen %d type %d\n",atom[hvy_atm].h_ptr,h_atm[atom[hvy_atm].h_ptr].typ);*/
                            }
                            
                            if (i<atom[hvy_atm].nh)
                            {
                                h_atm[atom[hvy_atm].h_ptr+i].p=vector_of(x,y,z);
                                h_atm[atom[hvy_atm].h_ptr+i].typ=fixed;
                            }
                        }
                        
                        break;
                    }/* matches 'if it is a H' */
                    else
                    {
                        printf("WARNING: Unrecognized atom name %4s from residue %3s \n",
                            atmnam, resnam);
                        /*formerly break;  ignore unrecognized atoms 2.24*/
                        atom[natoms].atmtyp  = -99;
                        atom[natoms].caindex = -99;
                        
                    }/* matches 'if it isn't an H' */
                    
                }/* matches 'if it isn't in necatm heavy atoms' */
            }
            else  /* atom is from a HETATM record of not-identified class */
            {
                if (debug)
                printf("\tUnidentified HETATM residue '%s'\n",resnam);
                
                if ( iswater(resnam) )
                {
                    if(debug)
                        printf("\tWater\n");
                    
                    if (atmnam[1]=='O')
                    {
/*                        printf ("*Found %s for ca# %d of %c \n",atmnam,namino,chain);*/
                        atom[natoms].nh=2;
                        atom[natoms].h_ptr=hkount;
                        hkount+=2;
                    }/* water O block ends*/
                    else
                    {
                        int i=0,hvy_atm;
                        if (debug) printf ("*Found H for ca# %d of %c \n",namino,chain);
                        hvy_atm=find_atom3(namino,chain," O  ",natoms-1);
                        if (hvy_atm>=0)
                        {
                            while(i<2 && h_atm[atom[hvy_atm].h_ptr+i].typ)
                                i++;
                            if (i<2)
                            {
                                if (debug) printf ("*Placed H %d for ca# %d of %c \n",i,namino,chain);
                                h_atm[atom[hvy_atm].h_ptr+i].p=vector_of(x,y,z);
                                h_atm[atom[hvy_atm].h_ptr+i].typ=fixed;
                            }
                            break;
                        }/* matching if-they-find-heavy-atom */
                        else
                        {
                            printf("WARNING: Heavy atom for pdb file Hydrogen not found.\n");
                            printf("%s\n",buf);
                        }
                    }/* water H or D block ends*/
                }/* water block ends */
                else /* HETATM, not identified, even as water */
                    ;
                /* just note the values, don't do anything */
                atom[natoms].atmtyp  = -99;/* also set to -99 earlier */
                atom[natoms].caindex = -99;
            }


            atom[natoms].atmnum   = atmnum;
            atom[natoms].p.x        = x;
            atom[natoms].p.y        = y;
            atom[natoms].p.z        = z;
            atom[natoms].occ      = occ;
            atom[natoms].b        = b;
            atom[natoms].aanum    = namino;
            atom[natoms].aacode   = aac;
/*            atom[natoms].strucsum = strucsum; done later 2.24*/
            atom[natoms].chnid    = (chain == ' ') ? '-' : chain;
            atom[natoms].altcode  = (altcode == ' ') ? '-' : altcode;
            atom[natoms].inscode  = inscode;
            atom[natoms].hetflg   = (token == HETATM);
            atom[natoms].ssflg    = FALSE;
            strcpy(atom[natoms].resnam, resnam);
            strcpy(atom[natoms].atmnam, atmnam);
            
            natoms++;
            break;


        default:                   /* Ignore all other types in this version */
            break;
        }   /* end switch */
        if (hkount>=MAXNATM-3)
          {
            sprintf(message,"Too many hydrogens: %d\0",hkount);
            fail(message);
          }
        debug=0;
    }


    printf("%d atoms selected from %d residues.\n\n", 
            natoms, reskount);
    if (natoms/* && !caonly*/)
/*        Removed !caonly for v3.14 after a basic check 
        that this wouldn't cause problems.*/
    {
/*        print_atoms("After reading PDB file\n");*/
                                    
        load_ststan();

        for(i=0;i<reskount;i++)/*set hetatm/sstruc of residue[], then atom[]*/
        {
            /*int numsstruc = 0;
            char last_sstruc = 0;*/
            
            residue[i].hetatm = (residue[i].ca==NULL);
            
            if (!residue[i].hetatm && sstfp)
            {
/*                printf("i=%d",i);*/
                
                while(fgets(sstbuf, 160, sstfp) && strcmp( atm2molaa( residue[i].strtat ), strnremspc(sstbuf+6,6) ) )
                {
                    printf("BUG: Line in sstfile %s differs from line in pdb input %s.\n",strnremspc(sstbuf+6,8), atm2molaa(residue[i].strtat) );
                    printf("SST: %.30s\n",sstbuf);
                    printf("PDB: %s\n",atomid(residue[i].strtat,buf));
                    
                }
/*                printf("%s matched with %s.\n",strnremspc(sstbuf+6,8), atm2molaa(residue[i].strtat) );*/
                
                
                residue[i].strucsum = sstbuf[25];
                if (first_ATOM_flg) /*starting up the first residue*/
                {
                    printf("Starting with %s.\n",atomid(residue[i].strtat,buf));
                    
                    sstrucseg[numsstruc].start = i;
                    sstrucseg[numsstruc].type = simplify_sstruc( sstbuf[25] );
                    first_ATOM_flg=0;
                    last_sstruc=sstbuf[28];
                    /*this stops the *normal* domain boundary commands from being triggered*/                                
                }
                
                if ((simplify_sstruc(last_sstruc) != simplify_sstruc(sstbuf[25]) && i)
                    || (islower(last_sstruc) && islower(sstbuf[25]))  
                    || (sstbuf[28]>' ' && sstrucseg_aggr>' ' && sstbuf[28] != sstrucseg_aggr) 
                    || (last_chnbrk_chr=='!' && sstbuf[5]=='!') )/* new sstrucseg */
                {
                    sstrucseg[numsstruc].stop = last_aares;
                    sstrucseg[numsstruc].aggr = sstrucseg_aggr;
                    
/*                    printf("Ending domain %d at #%d %s.\n",numsstruc,last_aares,atm2molaa( residue[i-1].strtat ));*/
              
                    numsstruc++;
/*                    printf("Beginning domain %d at #%d %s.\n",numsstruc,i,atm2molaa( residue[i].strtat ));*/
                    
                    sstrucseg[numsstruc].start = i;
                    sstrucseg[numsstruc].type = simplify_sstruc( sstbuf[25] );
                    if(sstrucseg[numsstruc].type=='T')
                        printf("BUG: T is sstrucseg[%d].typ.\n",numsstruc); 
                    sstrucseg_aggr=' ';
                }
                residue[i].sstrucnum = numsstruc;
                last_sstruc = sstbuf[25];
                last_chnbrk_chr=sstbuf[5];
                
                last_aares = i;
                if (sstbuf[28]>' ' && sstrucseg_aggr<=' ')
                    sstrucseg_aggr=sstbuf[28];
                
            }
            else
            {
                residue[i].strucsum = ' ';
/*                printf("#%d HETATM is %s\n",i,atm2molaa(residue[i].strtat));*/
                
            }
            
                    
            for(j=residue[i].strtat;j<=residue[i].stopat;j++)
            {
                if (j==0 && i > 0) printf("BUG: Residue %d/%d runs from %d to %d.\n",i+1,reskount,residue[i].strtat,residue[i].stopat);
                atom[j].hetflg= residue[i].hetatm;
                atom[j].strucsum = residue[i].strucsum;
                atom[j].residue_p = &(residue[i]);
            }
        }
        numsstruc++;
        sstrucseg[numsstruc-1].stop = last_aares;
        sstrucseg[numsstruc-1].aggr = sstrucseg_aggr;
        
        printf("Found %d elements of secondary structure and %d chains.\n", numsstruc, num_chains);
        
/*        for(i=0;i<numsstruc;i++)
        {
            printf("%5d %s+ %c%c",i,atomid( residue[sstrucseg[i].start].strtat,buf),sstrucseg[i].type,sstrucseg[i].aggr);
            printf("\n");
        }    */




        /* mc[] is an array of pointers to bits of the residue array. */
        /* only real amino acids get into mc[], not all residue */
        for(i=0;i<reskount;i++)
        {
            /*if ca is found, ie HETATM not set*/
            /* then increase aakount and add in the co-ordinates */
            
            if( !residue[i].hetatm )
            {
                residue[i].mc_p= &mc[aareskount];
                
                for(j=residue[i].strtat;j<=residue[i].stopat;j++)
                {
                    atom[j].mc_p= &mc[aareskount];
                    
            /*    if (sstfp) 
                          {
                          fgets(sstbuf, 160, sstfp);
                          atom[j].strucsum = sstbuf[25];
                          }
                          */  
                    switch(atom[j].atmtyp)
                    {
                    case 0:
                        mc[aareskount].n=j;
                        break;
                    case 1:
                        mc[aareskount].ca=j;
                        break;
                    case 2:
                        mc[aareskount].c=j;
                        break;
                    case 3:
                        mc[aareskount].o=j;
                        break;
                    default:
                        break;
                    }/*switch*/
                }/* match for j*/
                mc[aareskount].residue_p = &(residue[i]);
                aareskount++;
            }
            
            else /*not recognised aa*/
            {
                for(j=residue[i].strtat;j<=residue[i].stopat;j++)
                {
                    atom[j].caindex= -99;
                    atom[j].strucsum = space;
                }/* match for j */                
            }/*match if recognised aa*/
        }/*match for(i . . . )*/
        
        load_icon();
        find_brakes();
        if (debug) printf("Back in inpdb_file\n");
        find_ss();
        if (caonly)
/* Bug-fix. RAL 1 May 1997 --> */                
/*            printf("WARNING: C-Alpha Only File\nWARNING: Will continue to examne HETATMs\n"); */
/* Bug-fix. RAL 16 May 1999 --> */                
          {
/* <-- RAL 16 May 1999 */                
            printf("WARNING: C-Alpha Only or Nucleic Acid only File\n");
            printf("WARNING: Will continue to examine HETATMs\n");
/* Bug-fix. RAL 16 May 1999 --> */                
          }
/* <-- RAL 16 May 1999 */                
/* <-- RAL 1 May 1997 */                
        
/* Bug-fix. RAL 17 Aug 1999 --> */
        if (sstfp != NULL)
/* <-- RAL 17 Aug 1999 */
          fclose(sstfp);   
/* Bug-fix. RAL 17 Aug 1999 --> */
        if (ifp != NULL)
/* <-- RAL 17 Aug 1999 */
          fclose(ifp);
        
    }
    return(natoms /*&& !caonly*/);
    /* !caonly removed 3.13 to enable the program to work with non-Protein structures that are placed in a PDB format file */
}


