/* HBPLUS - Hydrogen Bond Calculator Header File v 3.2 */
/* Copyright I K Mcdonald, D T Jones, D Naylor and J M Thornton 1993
             All Rights Reserved
   Copies are free to Academic Users - contact mcdonald@uk.ac.ucl.bsm.bioc.bsm

 The publication of research using the Software must reference "McDonald
 IK, Naylor DN, Jones DT & Thornton J M (1993), 'HBPLUS', Computer Program,
 Department of Biochemistry and Molecular Biology, University College
 London." or successor references as defined by the authors.

 Unless informed otherwise, you should be an academic user and have sent a
 signed confidentiality agreement to the authors (address at the end of
 this file).  If you do not have a copy of HBPLUS and would like to receive
 one (free to academic users), please detach the confidentiality agreement
 from the end of this document, sign it, and send to the address given.
 Please allow other people in your department to use your copy of HBPLUS,
 but do not allow them to make their own copy. */

/* Contents */
/* 1. Copyright - all*/
/* 2. Version Log - all*/
/* 3. Tokens as abbreviations */
/* 4. Tokens as flags and customisables */
/* 5. Globals as flags and customisables */
/* 6. Globals as name arrays */
/* 8. Vector/Matrix Routines */
/* 8.5 Globals as arrays */
/* 8.7Minor useful functions */
/* 9. Hydrogen Position Calculation Subroutines */
/*10. Hydrogen Position Calculation Routine */
/*11. PDBOUT */
/*12. Connectivity Calculation Subroutines */
/*13. Connectivity Calculation Routine */
/*13.2 Domain calculation Routines */
/*13.3 Other BSMU results-generating code*/
/*13.5 Adding new residues to the selection */
/*14. find_hb */
/*15. inpdb_file */
/*16. main */

/****************************************************************************/
/* 2. Version Log */

#define VERSION "3.2"

/* Version 0.9     by David T. Jones, August 1990 */

/* Version 0.9a    modified by DN, 9th August 1991 */

/* Version 0.9b    modified by DN, 28th August 1991 */

/* Version 0.9c    modified by DN, 28th August 1991 
                   to accommodate insertion codes   */

/* Version 0.9d    modified by DN, 25th October 1991
                   to calculate gap properly        */

/* Version 0.9e    modified by DN, 30th October 1991
                   to treat water oxygens as both donors and acceptors */

/* Version 0.9f    modified by DN, 11th March 1992
                   to add non-standard amino acid */

/* Version 0.9g    modified by D. Naylor, 8th April 1992
                   to take account of non-standard amino-acids, 
                   missing atoms (non-standard chain breaks), and
                   ligand molecules. */

/* Version 1.0     by Ian McDonald 22nd April 1992
                   to find the position of hydrogens where possible,
  HBPLUS.C     	   and limit hydrogen bonds to those bonds where the angle
		   at the hydrogen is over a particular minimum, and the
		   distance between acceptor and Hydrogen is under a set limit.
		   */

/* Version 1.0h    DNs 13th April 1992 correction of "a few typos" */

/* Version 1.0i    DN 14th April 1992
                   fixed bug that caused program to ignore last residue */

/* Version 1.0j2   IM 20th June 1992
                   -o command line argument: output a pdb file */

/* Version 1.0h2   IM 28th July 1992
                   include HIS NE1 as donor */

/* Version 1.0i2   IM 6th August 1992
                   adjust angles around NH */

/* Version 1.0j3   IM 29th August 1992
                   add angle criteria at the acceptor */

/* Version 1.0k    IM 30th August 1992
                   ensure that 1-4 bonds are allowed.
		   allow CYS.ss to act as an acceptor
		   throw out CYS.ss SH Hydrogens when the SS bridges are
		   found in find_hb */

/* Version 1.0l    IM 31st August 1992
                   throw out all bonds where the hydrogen cannot be positioned
		   but the donor-acceptor distance is more than the allowed
		   hydrogen-acceptor distance plus one Donor-H bond length.
		   (set at one Angstrom) */

/* Version 1.0m    IM 19th September 1992
                   allow the -c option which refers to CYS SG atoms as either
		   CSS SG or CYH SG depending on whether they are Cystines or
		   Cysteines */

/* Version 1.0n    IM  9th October 1992
                   Check /data/pdb/prerelease/pdb????.ent as well as
		   /data/pdb/p????.pdb */

/* Version 1.0p    IM  7th November 1992
                   Tighten up the positioning of NHs on atoms with insertion
		   codes and the listing of pdbout (ie including said Hydrogen
		   positions) files. */

/* Version 1.0q    IM 23rd December 1992 (<- Hard Worker, eh ?)
                   Redo the lines to read hydrogens from files
		   change oracle to idata
		   change angle at OH of Tyr from 120 to 110 */

/* Version 1.0r    IM 27th December 1992
                   Remove un-needed debugging line */

/* Version 1.0s    IM 28th December 1992 redo OH Tyr angle */

/* Version 1.0t    IM  7th April 1993 fix dha,haaa,daaa separately */

/* Version 1.0u    IM 7th May 1993
                   The smaller angle at the acceptor is used, not the larger
		   Command line argument (-x) allow for H-Bonds with "wrong" 
		   atoms of Asn, Gln and His. */

/* Version 1.0v    IM 21st May 1993
                   Command line argument (-X) calculates only H-Bonds with
		   swappable atoms of Asn, Gln and His
		   Atom identifiers with a space in the third (/4) position are
		   now used complete. */

/* Version 1.0w    IM 2nd June 1993
                   longoutflg governs output format */

/* Version 2.0     IM 3rd June 1993
                   -Ii options control whether sstflag is important
		   options may be added together on the same line.
		   filenames (new,pdb) may be included in commandline 
		   arguments

   Version 2.01    IM 14th June 1993
                   Remove extraneous comman from REMARKS in *.h files

   Version 2.02    IM 13th July 1993
                   -X options also allows Prolines to accept

   Version 2.03    IM 16th July 1993
                   Return the -X option to what it was before

   Version 2.04    IM 22nd July 1933
                   Rationalise the *.hb2 and *.h headers
		   Remove hbdebug.dat after running copies that were compiled 
		   by compilers that did not have BSM predefined (ie
		   somewhere other than on the machine on which it was 
		   developed.)

   Version 2.05    IM 24th July 1993
                   add the -e and -E options to add extra donors (-E)
		   and acceptors (-e).

   Version 2.06    IM 28th July 1993
                   change of output format and order of functions

   Version 2.07    IM 28th July 1993
                   stop alternate locations from being counted as covalently 
		   bonded

   Version 2.10    IM  1st August 1993 

                   HBPLUS is three this month !  

                   Buoyed on by fact that someone actually wants to use
                   HBPLUS, I make several adjustments.

                   NEW Burke/Petsko interactions, and the final enabling of
                   the -R options.  

		   NEW Recogniton of C, A, U, G, T, Glx, Unk, ATP, FMN,
		   Methotrexate, CoA, Heme and NAD
		   
		   And replacing the posH table with deductions on the
		   basis of atom type, nearby connectivities, and the
		   geometry of the immediate antecedent.
 
   Version 2.11    P.A.Keller, (applied by IM) 1st August

                   Adding occupancies and B-values to the *.h file.
		   Thankyou to PAK for making the alteration available and
		   letting me "make any use of these alterations".

   Version 2.12    IM 8th August
                   Minor bugfix - getting the topological hydrogen placement
		   routine to work with the HQN exchange.

   Version 2.13    IM 11th August
                   Add a set of routines to make it easy to add new residues.

   Version 2.14    IM 26th August
                   Bugfix to -B option routine

   Version 2.15    IM  5th September
                   Bugfix to procinfilename

   Version 2.16    IM  7th October
                   Add "You must cite this" to output
		   Add "Nn" option (find neighbours instead of H-bonds).
		   Make a few flags "short" instead of "int"

   Version 2.17    IM 10th October
                   Make the "unusual" CA-CA distance dependent on "debug".

   Version 2.18    IM 13th October
                   Add an option to change how many covalent bonds count as
		   "nearly bonded".

   Version 2.19    IM 20th October
                   Tell the BSM version to look in ~/data/pdb for 
		   pdb files
		   
   Version 2.20    IM 29th November
                   Incorporate -f -T -U and -M atoms, allowing command-line 
		   options to include the structures of atoms as yet unknown 
		   to HBPLUS (ready for alpha testing?).

   Version 2.21    IM 4th December
                   bond-type (MSH/MSH thing) is "H" if .hetflg is on
		   previously it was only "H" if the atom was unrecog.
		   Debug -P (print list of donors/acceptors)
		   
   Version 2.22    IM 8th December
                   Preparing the Domains Routine - Introduce aareskount

   Version 2.24    IM 16th January 1994
                   Major division, into several files
                   HETATM/ATOM division now according to whether the CA
		   of the residue is located.
		   
		   An unrecognised residue no longer implies a chain break.

   Version 2.25    IM 21st January 1994
                   Further division.  The source files are now hbplus.h and
		   hbp_{gen,inpdb,findh,hhb,dom,main}.c .  hbp_dom.c is still 
		   left out of the public release.
		   
   Version 2.26    IM 24th January 1994
                   Add routines to hbp_inpdb.c and arrays and structures to
		   hbplus.h/hbp_gen.c to input secondary structure segments
		   from sstruc files.
		   
   Version 2.27    IM 30th January 1994
                   Assign domains to secondary structure segments.

   Version 2.28    IM 3rd February 1994
                   Take the first character as the element name in find_h if
		   the second character is a digit.

   Version 2.29    Add chain count and domplus chain option.

   Version 3.0     IM 15th June
                   
                   (additions and editings will from hereon be noted only
		   in this header file)
		   Include the hbp_qnh algorithm.
		   Note the atom#s of all H-bond partners to and from each 
		   atom.
		   Release
	   
	   3.01    IM Allow MSW options to be included

   Version 3.02    IM 25th July
                   Fix major bug in -U option
		   
   Version 3.03    IM  8th August
                   Program in .hbplusrc option
		   Rewrite add-donacc code to comply with manual
		   Fix add-bonds

   Version 3.04    IM 10th August
                   Increase TOTNATM to 100, from 60.

   Version 3.05    IM 26th August
                   Add *O3 as a donor to ATP at SLM's request. (hbp_inpdb)

   Version 3.06    IM  5th October
                   Remove ="???\0" for certain string initalisations in
                   hbp_gen.  Also add spaces to '*j++=*i;' to remove
		   ambiguity.

		   Change Bug warning for donor atoms not recognised as N,
		   O or S.  Remove a warning that a line should not be
		   reached (it merely caused compiler warnings about a line
		   that could not be reached).  
		   
   Version 3.07    IM 8th October 

                   Rationalise Bug Warnings for unusual numbers of covalent
		   bonds to donor atoms in find_h, remove NOTE: angle in
		   isplanar, add newline at start of output for each
		   protein.

   Version 3.08    IM 9th October
                      
                   Adjust hbp_qnh so that histidine side-chains covalently
		   bonding to other residues (theoretically metal ions)
		   register as hydrogen bonded under a '-X' Asn/Gln/His
		   analysis.  

   Version 3.09    IM 18th October
                   
                   Update location of default PDB files.

   Version 3.10    IM 7th November
                   
                   Install foolproofing and error-checking into the routine
                   to print out hydrogen bonds, now that it is separate
                   from the main hydrogen bond finding routine.

   Version 3.11    IM 9th November

                   Place spaces in '=&' and '=-' combinations.
		   (These combinations were confusing some compilers)

   Version 3.12    IM 11th November

                   Change the BSM spaces in line with the new /data structure
		   at /pdb.

   Version 3.13    IM 17th November

                   Add '-i' OMLINSERT option for Susan
		   Change ???, to suit
		   Allow program to continue for CA only structures.

   Version 3.14    IM 31st December
                   
                   Allow program to process CONECT records for structures
		   that register as 'CA only', as some of these structures
		   are actually 'DNA only' and quite 'H bond
		   calculation-able'.  For Andrew Martin

   Version 3.15    IM 1st January 1995
   
                   Introduce basic error catching for 'angle' function when
		   small floating point errors produce an error in
		   sprint_hb.  Error still unexplained, but disappears when
		   the call of 'angle' is moved to a different line to the
		   sprintf statement.

   Version 3.2     RAL 15th November 2007
   
                   Amendments to recognize new nucleotide definitions
                   following the wwPDB's remediation exercise.

*/
/****************************************************************************/
/* 3. Tokens as abbreviations */

#include <stdio.h>

/*#define BSM /* If this is defined, then the code compiles expecting to files in their locations in their home laboratory */

#define NO 0
#define YES 1

#define TF1 "%6.2f%6.2f" /* P.A. Keller 1993 */

#define TF "%8.3f%8.3f%8.3f"

#define VXYZ(v) v.x, v.y, v.z

#define SHORTID(d) atom[d].chnid,atom[d].aanum,atom[d].resnam,atom[d].atmnam

#define cos30 0.8660254

enum RESCODE
{
    Ala, Cys, Asp, Glu, Phe, Gly, His, Ile, Lys, Leu,
    Met, Asn, Pro, Gln, Arg, Ser, Thr, Val, Trp, Tyr,
    Aib, Phl, Sec, Alm, Mpr, Frd, Lym, Glm, Pph, Pgl,
    Ole, Aba, Nle, B2v, B2i, Blf, Bno, B2a, B2f, Iva,
    Lov, Sta, Pvl, Cal, Pha, Dci, Ahs, Chs, Mse, Eta,
    Pca, Asx, Glx, Unk, Cyh, Css,
    __C, __A, __U, __G, __T,
/* Amendment. RAL 15 Nov 2007 --> */
    _DC, _DA, _DU, _DG, _DT,
/* <--- RAL 15 Nov 2007 */
    ATP, CoA, FMN, Hem, Mtx, NAD,
    Ace, For
};

#define max(x,y) ((x)>(y) ? (x) : (y))
#define min(x,y) ((x)<(y) ? (x) : (y))

#define FALSE              0
#define TRUE               1

#define SQR(X)  ((X)*(X))

/****************************************************************************/
/* 4. Tokens as flags and customisables */

#define STDAA       20  /* no of standard amino acids recognized by hbplus   */
#define NONSTDAA    32  /* no of non-standard amino acids recognized by hb   */
#define MAXNAA     100  /* total no of rows in the amino acid tables */
extern short TOTNAA ; /*=69 7.1 * total no of amino acids recognized by neighbour   */
/* Amendment. RAL 25 May 1997 --> */
#define TOTNATM    200  /* total no of columns necatm/donor and other tables */
/* <-- RAL 25 May 1997 */
#define MAXCON       8  /* max no of bonds allowed for any 1 atom            */
#define MAXNHB_PER_ATM 40 /* to be stored in the hb_prtnr array */
#define MAXNHB_ATOMS  250000 /* Atoms that donate and accept */
/* #ifdef __MSDOS__
#define MAXNATM  1000  num atoms in atom[] 
#define MAXNRES   150  num residues in ?[] 
#define MAXBRKS   200  num chain breaks in 
#define MAXNRES_DOMAIN 0  num residues in k
#else */
/* Amendment. RAL 23 Jun 2008 --> */
/* #define MAXNATM 25000 */
#define MAXNATM 100000
/* <-- Amendment. RAL 23 Jun 2008 */
/* Amendment. RAL 12 Apr 2002 --> */
/* #define MAXNRES 1200  */
/* Amendment. RAL 23 Jun 2008 --> */
/* #define MAXNRES 6000 */
#define MAXNRES 24000
/* <-- Amendment. RAL 23 Jun 2008 */
/* <-- RAL 12 Apr 2002 */
/* Amendment. RAL 23 Feb 2004 --> */
/* #define MAXBRKS  2000 */ /*num chain brakes in array */
#define MAXBRKS  5000 /*num chain brakes in array */
/* <-- RAL 23 Feb 2004 */
#define MAXNRES_DOMAIN 200 /* num residues in kshb[] */
/* Amendment. RAL 9 Mar 2009 --> */
#define FILENAME_LEN   1024
/* <-- RAL 9 Mar 2009 */
/* Amendment. RAL 3 Jul 2012 --> */
#define LINE_LEN        512
#define MESSAGE_LEN     256
/* <-- RAL 3 Jul 2012 */

/****************************************************************************/
/* 5. Globals as flags and customisables */
extern int check;
extern int ssnum;
extern float sinMIN_DHA ;
extern float cosMIN_DHA ; /* this means the dot product of the two vectors
			would come above zero if it below this minimum*/
extern double mindha  ;	/*angle */

extern float sinMIN_DAAA ;
extern float cosMIN_DAAA ; /* this means the dot product of the two vectors
			would come above zero if it below this minimum*/
extern double mindaaa  ;	/*angle */

extern float sinMIN_HAAA ;
extern float cosMIN_HAAA ; /* this means the dot product of the two vectors
			would come above zero if it below this minimum*/
extern double minhaaa;	/*angle */


extern float sinMAX_DAAX ;
extern float cosMAX_DAAX ;
extern double maxdaax ;

extern float sinMAX_HAAX ;
extern float cosMAX_HAAX ;
extern double maxhaax ;

extern float MAX_HA ;

extern float CAWARN ;

extern short pdboutflg;
extern short cssflg ; /* How to represent CYS in output - CYS or CYH */
extern short exchangeflg; /* =0 HIS, ASN, GLN hydrogen bond normally */
                   /* =1 HIS, ASN, GLN hydrogen bond as if swappable */
                   /* =2 Only the bonds with swappable atoms are counted
		    and the output is given*/
extern short nnbflg; /* =0 to output hhb/hb2 file */
                /* =1 to output nnb/nb2 file */
/*extern short domainflg; cancelled in the great split :) */
extern short kshflg; /* =1 to use Kabsch/Sander positions of Hydrogens */
extern short asaflg; /* =1 to input an asa or buratm file */

extern short longoutflg; 
extern short inputsstflg; /* =0 don't worry about sst, =1 do try to input them */
extern char only_chainid_lst[32];

extern float HBDIST ;
extern float SSDIST;

extern short numcovbonds ; /* number of covalent bonds that count as "nearly bonded" */

extern short bsmoptflg; /* v3.01 */
/* =0 no special bsm options, 1 = output MSW files */
extern short OMLINSERT;
/* v3.13 =0 normal =1 replace '-' with 'h' for hetatm insertion codes*/


/****************************************************************************/
/* 6. Globals as name arrays */

/* Residue name to allow conversion of a.a. name into numeric code */
extern char rnames[MAXNAA][4] ;

/*****************************************************************************/
/* Record names for decoding record types */
/* Bug-fix. RAL 1 May 1997 --> */		
/* extern const char           *tokstr[25] ; */
extern const char           *tokstr[26] ;
/* <-- RAL 1 May 1997 */		

/* Bonding Atom Dictionary - this is the list of atom names, arranged on a
   residue by residue basis in RESCODE order. 
*/
extern char necatm[MAXNAA][TOTNATM*4] ;

extern short accepts[MAXNAA][TOTNATM] ;

extern short donors[MAXNAA][TOTNATM] ;

/****************************************************************************/
/* 8. Vector/Matrix Routines */

/* these are the vector arithmetic routines I wrote on 31/1/91 */

/* firstly a vector structure is important */
struct vect
{
    float x,y,z;
};

struct loci
{
    short typ;
    struct vect p,a;
};

struct mat3x3 /* a three by three matrix */
{
    struct vect r1,r2,r3;
};
/* note that I premultiply a column vector by a transformation matrix, 
   NOT postmultiply a row vector. */

/*perpendicular, unit_vector, vector_plus, float_times_vect*/
struct vect vector_array(float * array);
struct vect vector_of(float x, float y, float z);
struct vect perpendicular(struct vect a, struct vect b, struct vect c);
struct vect vector_plus(struct vect a, struct vect b);
struct vect float_times_vect(float f, struct vect v);
struct vect unit_vector(struct vect);
struct vect cross_product(struct vect a, struct vect b);
struct vect intersect(struct vect normal, struct vect point, struct vect first, struct vect last);
struct vect onto_sphere(struct vect point, struct vect centre, float radius);
struct vect from_a_to_b(struct vect a, struct vect b);
#define to from_a_to_b
#define unit unit_vector
float length_squared(struct vect a);
double double_length_squared(struct vect a);
double dot_product(struct vect a, struct vect b);
float vector_length(struct vect a);
double double_vector_length(struct vect a);
float angle(struct vect a, struct vect b, struct vect c);
struct vect mlt_mtx_vct(struct mat3x3 pre, struct vect post);
struct mat3x3 axs_rot(struct vect x_axis, struct vect y_axis, struct vect z_axis);


/**************************************************************************/
/* 8.5 Globals as global variables */

extern struct hb_prtnr
{
    short unsigned int n[MAXNHB_PER_ATM]; 
} hb_prtnr[MAXNHB_ATOMS];

/* Each structure is a list of HB partners that can be pointed to by either 
   don_p or acc_p in a atom[] structure.

See also pdbatm.ndonhb, pdbatm.nacchb, pdbatm.don_p pdbatm.acc_p */

extern struct pdbatm 
{
    struct vect   p; /* position */
    float         occ, b; /* pak --- occupancy and B factor */
    float         acc; /* accessibility -99.9, 0.0-1.0, 2.0 */
    int           aanum ;/*residue count from pdb file:namino*/
    int           aacode/*01-*/, atmtyp/*as recognised from necatm, def -99*/;
    int           caindex /*HETs increment this count but have it set to -99 */; 
    int           atmnum /*from pdb*/, ncon/*number of connectivity records*/, h_ptr, nh;
    int           ndonhb,nacchb;
    char          altcode, inscode, chnid /*- is default for alt and chn*/, strucsum, 
                  resnam[4], atmnam[5] /*obvious*/;
    short         hetflg,ssflg;/*TRUE/FALSE*/
    struct main_chain * mc_p;
    struct pdbres * residue_p;
    struct hb_prtnr * don_p, * acc_p;

/* The hetflg is set . . . */
}   atom[MAXNATM];

/* h_ptr is the element number in the hydrogen array of the first hydrogen
   connected to that hydrogen - nh is the number of Hydrogens.*/

extern struct loci h_atm[MAXNATM]; /*array of loci of all hydrogen */

extern struct pdbres
{
    struct vect *n, *ca, *c, *o;
    int strtat,stopat;
    short hetatm;
    char strucsum;
    int sstrucnum;
    struct main_chain * mc_p;
} residue[MAXNRES];

/* main-chain amino acid co-ordinates */

extern struct main_chain{
    int c,o,ca,n;
    struct pdbres *residue_p;  /* pointer to place in residue[]*/
} mc[MAXNRES];

/* the atom array ordinates of the main_chain heavy atoms */

extern struct pdbsstruc 
{
    int start, stop; /* in CAINDEX values */
    char type; /* in simplified KS+ [HEGP ] */
    short aggr; /* eg which beta sheet ? */
} sstrucseg[MAXNRES];

extern short   brake[MAXBRKS][2];
/* allof these are used by in order to find chain-breaks */

/* used to help calculate connectivities */

/*char    *pdbfn[160], csdfn[160], logfn[160], keyword[40], buf[160];
  not logical as global variables, csdfn and logfn are not even mentined*/
extern char    brcode[5] /*, inpdbfn[160]*/;/*headerise*/
extern char    hbdn[160]; /* only used in main and procfile. HB directory */

extern FILE      *ofp ;/* HB2 file declared - defined hbp_hhb.c */
extern FILE      *dbgfp; /* new in 0.9m */

extern int       natoms, debug, nbrakes, nconrecs, reskount /*headerise*/;
extern int       numsstruc;
extern int       num_chains; /*v2.29*/
extern int       aareskount ; /* new v2.2, so that find_domains only deal with real aas. */
extern int       icon[MAXNATM][MAXCON]; /* 0.9m */ /*the REAL connectivity*/
/*extern int       bonds[MAXCNRECS][MAXCNCOLS] ;  0.9m */ /*copy of CONECT records - internal to hbp_inpdb.c*/

/****************************************************************************/
/****************************************************************************/
/* 8.7 Minor useful functions */

/***************************************************************************/
char * atomid(int a, char * string);
short iswater(char * resid);
short isnucleotide(char * resid);
short isplanar(int i);
void fail(char * errstr);

FILE * fmfopen(char* fn, char * mode); /* fopen with autodestruct */

char * instr(char *ct, char *cs);
int scanatmnam(char * atmnam,int res);
int scanresatmnam(char inatmnam[8],int * res, int * atm);
int scanresnam(char * inresnam);
void print_atoms(char * string) /*debugging only*/;

int find_atom(int caindex, char chnid, int atmtyp, int atomnum);
int find_atom2(int caindex, char chnid,char * atmnam, int atomnum);
int find_atom3(int aanum, char chnid,char * atmnam, int atomnum);
int find_atom4(char chnid, int aanum, char inscode,char * resnam, char * atmnam, int startindex);

char * atm2molaa(int i);
char * strnremspc(char * string, int n);

int find_gap(int caindex_i, int caindex_j);

void enable_exchange(void);
void enable_aromatic(void);
void disable_aromatic(void);

struct iif /* transferred here from hbp_dom.c 2.24 */
{
    int a,b;
    float f;
};
void parsefn(char * fname, char * rootnam); 

int inasa_file(char * asafname);
void printf_numhb(void);
void chkqnh(void);
void supplement_arrays(void);
void initialise_arrays(void);


/**********************************************************************/
/** 9.  Hydrogen Position Calculation Subroutines
/****
/*** HLOCI is used in interpreting hydrogen positions */
enum HLOCI 
{
    fixed=1, alternatives, circle 
};
/* 10. Hydrogen Position Calculation Routine */
void find_h(void);



/*11. PDBOUT */
void pdbout(char * outfn,char * pdbname, short int filter_flg,short int *filter);

/*12. Connectivity Calculation Routines*/
int alreadybonded(int atom1,int atom2);
int bondedwithin(int atom1, int atom2, int nbonds);

/* 13.2 Domain Calculation Routines */
#ifdef DOM
extern float PEPBNDVAL /*= -0.50*/;
extern float VDVVAL  /*-5.00*/;
extern float LTLDOMBON /*= -2.00*/;
extern float HBVAL  /*1.00*/;
void find_domains();
#endif

/*13.3 Other code for generating results for BSMU research*/
#ifdef BSM
void do_msw(char * basename);
#endif

/*13.5 Adding new residues to the selection*/
int add_residue_type( char * residue, char * similar_residue );
int add_atoms( char * residue, char * atoms);
int add_bonds( char * residue, char * new_bonds);

/*14. find_hb */
void find_hb (char * fname);
void sprintf_hb(char * o_str, int i, int j, struct vect nearest, float d, float ha, float daaa_ang, float haaa_ang);
void fprintf_hb(char * o_str, int i, int j);
int calc_hb(char * o_str, int i, int j);



/*15. inpdb_file*/
short inpdb_file(char * fname, char * inpdbfn);
int add_donacc ( char * residue, char * atoms, short number, short donflg );
int  alreadybonded(int atom1,int atom2);



/* The Confidentiality Agreement (for up-to-date copy see the bsm ftp site)
Correspondence to:
Ian McDonald (PG)
Biological Structure and Modelling Unit
Department of Biochemistry and Molecular Biology
University College London
Gower Street
LONDON WC1E 6BT
UK / EC

Email mcdonald@uk.ac.ucl.bsm


                    HBPLUS - Hydrogen Bond Calculation
                    ----------------------------------          

			CONFIDENTIALITY AGREEMENT
			-------------------------



In regard to the HBPLUS program, specified in Appendix 1 herewith (the
Software) supplied to us, the copyright and other intellectual property
rights to which belong to the authors, we

    __________________________________________________________________

undertake to the authors that we shall be bound by the following terms and
conditions:-

1. We will receive the Software and any related documentation in confidence
and will not use the same except for the purpose of the department's own 
research. The Software will be used only by such of our officers or
employees to whom it must reasonably be communicated to enable us to
undertake our research and who agree to be bound by the same confidence.
The department shall procure and enforce such agreement from its staff for
the benefit of the authors.

2. The publication of research using the Software must reference "McDonald
IK, Naylor DN, Jones DT & Thornton J M (1993), 'HBPLUS', Computer Program,
Department of Biochemistry and Molecular Biology, University College
London." or successor references as defined by the authors.

3. Research shall take place solely at the department's premises at

    __________________________________________________________________

4. All forms of the Software will be kept in a reasonably secure place to
prevent unauthorised access.

5. Each copy of the Software or, if not practicable then, any package
associated therewith shall be suitably marked (and such marking maintained)
with the following copyright notice: "Copyright 1991-3 Ian McDonald, Dorica
Naylor, David Jones and Janet M Thornton All Rights Reserved".

6. The Software may be modified but any changes made shall be made
available to the authors.

7. The Software shall be used exclusively for academic teaching and
research. The Software will not be used for any commercial research or
research associated with an industrial company.

8. The confidentiality obligation in paragraph one shall not apply:

   (i)  to information and data known to the department at the time of
	receipt hereunder (as evidenced by its written records);

  (ii)	to information and data which was at the time of receipt in the 
	public domain or thereafter becomes so through no wrongful act of
	the department;

 (iii)	to information and data which the department receives from a third
	party not in breach of any obligation of confidentiality owed to
	the authors.



Please sign this Undertaking and return a copy of it to indicate that you 
have read, understood and accepted the above terms.



		      For and on behalf of _____________________________

		      _________________________________________________
		     
		      ..................................................

		      Dated ............................................


                      Address __________________________________________

		      _________________________________________________
		     
		      _________________________________________________
		     
		      Country _________________ Postcode ______________
		     
                      Telephone ________________________________________

                      Electronic Mail Address to which HBPLUS shall be

                      sent _____________________________________________

		      _________________________________________________
		     

APPENDIX 1 - DETAILS OF THE HBPLUS PROGRAM PROVIDED
---------------------------------------------------

Files to be included
--------------------

         1. hbplus.h              }
         2. hbp_gen.c             } Source program files,
         3. hbp_inpdb.c           } ,formerly hbplus.c
         4. hbp_findh.c           }
         5. hbp_hhb.c             }
         6. hbp_main.c            }
         7. hbplus.man         } Documentation

*/

