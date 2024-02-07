
/* TESTING on 1a23 */


/* HBPLUS - Hydrogen Bond Calculator - General v 3.2 */
/* Copyright I K Mcdonald, D T Jones, D Naylor and J M Thornton 1993
             All Rights Reserved
   Copies are free to Academic Users - contact mcdonald@uk.ac.ucl.bsm.bioc.bsm

 The publication of research using the Software must reference "McDonald
 IK, Naylor DN, Jones DT & Thornton J M (1993), 'HBPLUS', Computer Program,
 Department of Biochemistry and Molecular Biology, University College
 London." or successor references as defined by the authors.

 Unless informed otherwise, you should be an academic user and have sent a
 signed confidentiality agreement to the authors (address at the end of the
 source code).  If you do not have a copy of HBPLUS and would like to
 receive one (free to academic users), please detach the confidentiality
 agreement from the end of this document, sign it, and send to the address
 given.  Please allow other people in your department to use your copy of
 HBPLUS, but do not allow them to make their own copy. */

/* Contents */
/* 1. Copyright **/
/* 2. Version Log **/
/* 3. Tokens as abbreviations **/

/* 4. Tokens as flags and customisables */
/* 5. Globals as flags and customisables */
/* 6. Globals as name arrays */

/* 8. Vector/Matrix Routines **/

/* 8.5 Globals as globals **/

/* 8.7Minor useful functions **/

/* 9. Hydrogen Position Calculation Subroutines */
/*10. Hydrogen Position Calculation Routine */
/*11. PDBOUT */
/*12. Connectivity Calculation Subroutines */
/*13. Connectivity Calculation Routine */
/*13.2 Domain calculation Routines */
/*13.5 Adding new residues to the selection */
/*14. find_hb */
/*15. inpdb_file */
/*16. main */

/****************************************************************************/
/* 2. Version Log */

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
  HBPLUS.C                and limit hydrogen bonds to those bonds where the angle
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
                   found in find_hb

/* Version 1.0l    IM 31st August 1992
                   throw out all bonds where the hydrogen cannot be positioned
                   but the donor-acceptor distance is more than the allowed
                   hydrogen-acceptor distance plus one Donor-H bond length.
                   (set at one Angstrom)

/* Version 1.0m    IM 19th September 1992
                   allow the -c option which refers to CYS SG atoms as either
                   CSS SG or CYH SG depending on whether they are Cystines or
                   Cysteines

/* Version 1.0n    IM  9th October 1992
                   Check /data/pdb/prerelease/pdb????.ent as well as
                   /data/pdb/p????.pdb

/* Version 1.0p    IM  7th November 1992
                   Tighten up the positioning of NHs on atoms with insertion
                   codes and the listing of pdbout (ie including said Hydrogen
                   positions) files.

/* Version 1.0q    IM 23rd December 1992 (<- Hard Worker, eh ?)
                   Redo the lines to read hydrogens from files
                   change oracle to idata
                   change angle at OH of Tyr from 120 to 110

/* Version 1.0r    IM 27th December 1992
                   Remove un-needed debugging line */

/* Version 1.0s    IM 28th December 1992 redo OH Tyr angle */

/* Version 1.0t    IM  7th April 1993 fix dha,haaa,daaa separately */

/* Version 1.0u    IM 7th May 1993
                   The smaller angle at the acceptor is used, not the larger
                   Command line argument (-x) allow for H-Bonds with "wrong" 
                   atoms of Asn, Gln and His.

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

   Version 2.23    IM 4th January 1994
                   Separating into header files
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

/****************************************************************************/
/* 8. Vector/Matrix Routines */


/* these are the vector arithmetic routines I wrote on 31/1/91 */

/* firstly a vector structure is important */
/* struct vect - see hbplus.h */
/* and vectors will sometimes be put in arrays */

struct vect vector_array (float * array)
{
    struct vect vector;
    vector.x = *array;
    vector.y = *(array+1);
    vector.z = *(array+2);
    return (vector);
}

struct vect vector_of(float x, float y, float z) /*1/4/92 */
{
    struct vect vector;
    vector.x = x;
    vector.y = y;
    vector.z = z;
    return (vector);
}


/* subtracting */
struct vect from_a_to_b (struct vect a, struct vect b)
{
    struct vect difference;
    difference.x = b.x - a.x;
    difference.y = b.y - a.y;
    difference.z = b.z - a.z;
    return (difference);
}

/* dot product */
double dot_product (struct vect a, struct vect b)
{
    return (a.x*b.x+a.y*b.y+a.z*b.z);
}

/* vector length */

float vector_length(struct vect a)
{
    double square, root;
    square= a.x*a.x + a.y*a.y + a.z*a.z;
    root = sqrt(square);
    return ((float)root);
}

double double_vector_length(struct vect a)
{
    double square;
    square= a.x*a.x + a.y*a.y + a.z*a.z;
    return (sqrt(square));
}  

float length_squared(struct vect a)
{
    return (a.x*a.x + a.y*a.y + a.z*a.z);
}

double double_length_squared(struct vect a)
{
    return (a.x*a.x + a.y*a.y + a.z*a.z);
}


/* add two vectors */
struct vect vector_plus(struct vect a, struct vect b)
{
    struct vect sum;
    sum.x=a.x+b.x;
    sum.y=a.y+b.y;
    sum.z=a.z+b.z;
    return (sum);
}

/* multiply a vector by a scalar */
struct vect float_times_vect(float f, struct vect v)
{
    struct vect product;
    product.x= f*v.x;
    product.y= f*v.y;
    product.z= f*v.z;
    return (product);
}

/* return the unit vector of . . . */
struct vect unit_vector(struct vect v)
{
    return (float_times_vect( 1/vector_length(v),v));
}

/* angle - give the angle between two lines */
/* watchit - this returns an angle in radians */
float angle(struct vect a, struct vect b, struct vect c)
{
    double ab2,bc2, dot;
    float result ;
    
    ab2=double_length_squared(to(a,b));
    bc2=double_length_squared(to(b,c));
    dot=dot_product(to(b,a), to(b,c));
    dot=dot/sqrt(ab2*bc2);
    if(dot< -1)
        dot = -1;
    if(dot>1)
        dot = 1;
    result = acos(dot);
    if (result > 3.1415927 || result < -3.1415927)
        printf("BUG: Arithmetic error in function angle. result = %f\n",
               result);
    
    return (result);
}

/* cross product - gives a vector perpendicular to both the arguments*/
struct vect cross_product(struct vect a, struct vect b)
{
    /*struct vect cross_product,origin={0.0,0.0,0.0};
    
    cross_product =  (vector_of (a.y*b.z-a.z*b.y, -a.x*b.z+a.z*b.x, a.x*b.y-a.y*b.x));
    printf("Cross_product gives angles %5.3f and %5.3f.\n",angle(a,origin,cross_product),angle(b,origin,cross_product));
    return( cross_product );*/
    return  (vector_of (a.y*b.z-a.z*b.y, -a.x*b.z+a.z*b.x, a.x*b.y-a.y*b.x));
    
}

/* find the point of intersection of a line and a plane */
struct vect intersect(struct vect normal, /* line perp to the plane */
                      struct vect point, /* point in the plane */
                      struct vect first, /* of the line */
                      struct vect last /* of the line */)
{
    float k /* as in normal dot point equals v */,micro, denom;
    struct vect line;
    line= from_a_to_b(first,last);
    denom= dot_product(normal, line);
    if (denom==0.0) denom=0.0000001 /*this is cheating*/;
    k=dot_product(normal,point);
    micro=( (k-dot_product(normal,first))/denom );
    return (vector_plus(first, float_times_vect(micro,line) ) );
}

/* translate a point onto the surface of a sphere */
struct vect onto_sphere(struct vect point, struct vect centre, float radius)
{
    struct vect centre_to_point;
    centre_to_point=from_a_to_b(centre,point);
    centre_to_point=float_times_vect( radius, unit(centre_to_point) );
    return (vector_plus (centre,centre_to_point) );
}




struct vect mlt_mtx_vct(struct mat3x3 pre, struct vect post)
{
    return ( vector_of( dot_product(pre.r1,post), dot_product(pre.r2,post), 
              dot_product(pre.r3,post) ) );
}

struct mat3x3 axs_rot(struct vect x_axs, struct vect y_axs, struct vect z_axs)
{
    struct mat3x3 result;
    result.r1=vector_of(x_axs.x, y_axs.x, z_axs.x );
    result.r2=vector_of(x_axs.y, y_axs.y, z_axs.y );
    result.r3=vector_of(x_axs.z, y_axs.z, z_axs.z );
    return (result);
}

struct vect perpendicular( struct vect a, struct vect b, struct vect c)
{
/*    struct vect var, o_a, o_c, origin={0.0,0.0,0.0};*/
    
    return ( ( cross_product( to(b,a) , to(b,c) ) ) );
/*    printf("Perpendicular gives angles %5.3f and %5.3f.\n",angle(a,b,vector_plus(b,var)),angle(c,b,vector_plus(b,var)));
    return(var);*/
    
}

/**************************************************************************/
/* 8.5 Globals as global variables */

struct pdbatm atom[MAXNATM];
struct hb_prtnr hb_prtnr[MAXNHB_ATOMS];

/* h_ptr is the element number in the hydrogen array of the first hydrogen
   connected to that hydrogen - nh is the number of Hydrogens.*/

struct loci h_atm[MAXNATM]; /*array of loci of all hydrogen */

/* struct pdbres . . . residue - hbplus.h */
struct pdbres residue[MAXNRES];

struct pdbsstruc sstrucseg[MAXNRES];

/* main-chain amino acid co-ordinates */
struct main_chain mc[MAXNRES];

/* the atom array ordinates of the main_chain heavy atoms */


short   brake[MAXBRKS][2];
/* allof these are used by in order to find chain-breaks */

/* used to help calculate connectivities */

/*char    *pdbfn[160], csdfn[160], logfn[160], keyword[40], buf[160];
  not logical as global variables, csdfn and logfn are not even mentined*/
char    brcode[5] /*, inpdbfn[160]*/;

FILE      *dbgfp=NULL; /* new in 0.9m */

short exchangeflg=0; /* =0 HIS, ASN, GLN hydrogen bond normally */
                   /* =1 HIS, ASN, GLN hydrogen bond as if swappable */
                   /* =2 Only the bonds with swappable atoms are counted */
short bsmoptflg=0; /*v 3.11 =1 MSW options*/
short OMLINSERT = 0; /*v 3.13 =1 replace '-' with 'h' for HETATM */

int       natoms, debug=0, nbrakes, nconrecs, reskount ;
int       numsstruc;
int       num_chains=0; /*v2.29*/
int       aareskount=0 ; /* new v2.2, count of aminoacid residues */
int       icon[MAXNATM][MAXCON]; /* 0.9m */ /*the REAL connectivity*/

/****************************************************************************/

/****************************************************************************/
/* 8.7 Minor useful functions */

/***************************************************************************/

char * atomid(int a,char * string)
{
    /*new treatment of inscode at 3.13*/
    char inscode;
    inscode = (OMLINSERT && atom[a].hetflg)?atom[a].inscode:'h';
    
    sprintf(string,"%.4s/%c%.04d%c%.3s%.4s",brcode,atom[a].chnid,
                  atom[a].aanum, inscode ,atom[a].resnam,
            atom[a].atmnam);
    return(string);
}



short iswater(char * resid)
{
    char tmpstr[4];
    strcpy(tmpstr,resid);

    return( NULL != instr("H2O:HHO:OHH:HOH:OH2:SOL:WAT:DOD:DOH:HOD:D2O:DDO:ODD:OD2:HO1:HO2:HO3:HO4",tmpstr) );
}

short isnucleotide(char * resid)
{
    char tmpstr[4];
    tmpstr[3]='\0';
    
    strncpy(tmpstr,resid,3);
    return( NULL != instr("  C:  A:  U:  G:  T:",tmpstr) );
}

short isplanar(int i)
{
    int j;
    float mean_angle=0;
    char atomidbuf[20];
    
    for(j=0;j<atom[i].ncon;j++)
    {
        if (icon[i][j] < 0 || icon[i][j] > natoms)
        {
            printf("BUG: Connectivities not set up for %s.\n",atomid(i,atomidbuf));
            return(0);
        }
    }
    switch( atom[i].ncon )
    {
    case 1:
        printf("BUG: Only one bond to %s in isplanar.\n",atomid(i,atomidbuf));
        break;
    case 2:
        mean_angle = angle( atom[ icon[i][0] ].p, atom[ i ].p, atom[ icon[i][1] ].p);
        break;
    case 3:
        mean_angle = angle( atom[ icon[i][0] ].p, atom[ i ].p, atom[ icon[i][1] ].p);
        mean_angle += angle( atom[ icon[i][1] ].p, atom[ i ].p, atom[ icon[i][2] ].p);
        mean_angle += angle( atom[ icon[i][2] ].p, atom[ i ].p, atom[ icon[i][0] ].p);
        mean_angle /= 3;
#ifdef BSM
        if (debug)
            printf("NOTE: mean_angle %5.3f %s\n",mean_angle,atomid(i,atomidbuf));
#endif        
        break;
    default:
        printf("BUG: %d bonds to %s in isplanar.\n",atom[i].ncon, atomid(i,atomidbuf));
        break;
    }
    return (mean_angle > 2.065); /* cos 120= -1/2, cos 109.5= -1/3
                                       cos 2.00rads = - 5/12 */
}


void fail(char * errstr)
{
    printf("\n");
    printf("\nFATAL ERROR: %s\n\n", errstr);
#ifndef BSM    
/* Amendment. RAL 9 Mar 2009 --> 
    remove("hbdebug.dat");
 <-- RAL 9 Mar 2009 */
#endif    
    exit(-1);
}

/**********************************************************/

FILE *
    fmfopen(char* fn, char * mode) /* fopen with autodestruct */
{
    FILE * fp = NULL;
    fp=fopen(fn,mode);
    if (fp==NULL)
    {
        fprintf(stderr,"\n*** Cannot open %s file in mode %s.\n",fn,mode);
#ifndef BSM    
/* Amendment. RAL 9 Mar 2009 -->
        remove("hbdebug.dat");
 <-- RAL 9 Mar 2009 */
#endif    

        exit(-1);
    }
    return fp;
}


char *instr(char *ct, char *cs)
{
    int debug=0;
    
    int   l = strlen(cs);

    /* Input to INSTR consists of two strings CT and CS. It returns the address
       of the position in CT of the occurrence of the first character of the
       substring CS (if it is present) or zero. */
    
    if (debug) printf("'%s' '%s' %d\n", ct, cs, l); 

    for (; *ct; ct++)
        if (!strncmp(ct, cs, l))
            return (ct);
    return (NULL);
}

/****************************************************************************/
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
int find_atom(int caindex, char chnid, int atmtyp, int atomnum)
{
    int i,debug=0;
    atomnum=(atomnum<0)?0:atomnum;
    /*if (atomnum == 11)
    {
        printf("XXXX/%c%04dX%03d\n",chnid,caindex,atmtyp);
        debug=0;
    }
    else debug = 0;*/
    
    
    if (debug) printf("Chain ID %c, comparing with %c.\n",atom[atomnum].chnid,chnid);
    
    if (atom[atomnum].caindex == caindex && atom[atomnum].atmtyp == atmtyp && atom[atomnum].chnid == chnid )
        return (atomnum);
    for(i=(atomnum<10)?0:atomnum-10; i < atomnum+10 && i < natoms; i++)
        if (atom[i].caindex == caindex && atom[i].atmtyp == atmtyp && atom[i].chnid == chnid)
            return (i);
    for(i=0; i < natoms ; i++)
        if (atom[i].caindex == caindex && atom[i].atmtyp == atmtyp && atom[i].chnid == chnid)
            {
                if (debug) printf("0000/%c%04d %03d\n",atom[i].chnid,atom[i].caindex,atom[i].atmtyp);
                return (i);
            }
    

    if (debug) printf("Atom %d not found\n",atomnum);
    
    return (-99); /* if the atom has not been found */
}

int find_atom2(int caindex, char chnid,char * atmnam, int atomnum)
{
    int i,debug=0;
    atomnum=(atomnum<0)?0:atomnum;
    chnid=(chnid==' ')?'-':chnid;
    /* this assumes that the blank chainids are stored as '-' */
    if (debug) printf("Inside find_atom2\n");
    if (debug) printf("Chain ID %c, comparing with %c.\n",atom[atomnum].chnid,chnid);
    if (atom[atomnum].caindex == caindex && !strncmp(atom[atomnum].atmnam+2,atmnam+2,2) && atom[atomnum].chnid == chnid)
        return (atomnum);
/* Amendment. RAL 11 Dec 2013 --> */
    /* Loop over the nearby atoms */
    for(i=(atomnum<10)?0:atomnum-10; i < atomnum+10 && i < natoms; i++)
      {
        /* Identify heavy atom bound to hydrogen using PDB v.3.0 format */
        if (atom[i].caindex == caindex && atom[i].chnid == chnid)
          {
            /* Check atom name if 'H' is in first position */
            if (atmnam[0] == 'H')
              {
                if (!strncmp(atom[i].atmnam+2,atmnam+1,2))
                  {
                    return (i);
                  }
              }
            else if (atmnam[1] == 'H')
              {
                if (atom[i].atmnam[0] == atmnam[0] &&
                    (atom[i].atmnam[2] == atmnam[2] ||
                     (atom[i].atmnam[3] == atmnam[3] && atmnam[3] == ' ')))
                  {
                    return (i);
                  }
              }
          }
      }
/* <-- RAL 11 Dec 2013 */
    for(i=(atomnum<10)?0:atomnum-10; i < atomnum+10 && i < natoms; i++)
        if (atom[i].caindex == caindex && !strncmp(atom[i].atmnam+2,atmnam+2,2) && atom[i].chnid == chnid)
            return (i);
    for(i=0; i < natoms ; i++)
        if (atom[i].caindex == caindex && !strncmp(atom[i].atmnam+2,atmnam+2,2) && atom[i].chnid == chnid)
            return (i);
    if (debug) printf("Atom %d, %s in residue %d%c not found\n",atomnum,atmnam,caindex, chnid);
    
    return (-99); /* if the atom has not been found */
}    

int find_atom3(int aanum, char chnid,char * atmnam, int atomnum)
{
    int i,debug=0;
    atomnum=(atomnum<0)?0:atomnum;
    chnid=(chnid==' ')?'-':chnid;
    /* this assumes that the blank chainids are stored as '-' */
    if (debug) printf("Inside find_atom3\n");
    
    if (debug) printf("Chain ID %c, comparing with %c.\n",atom[atomnum].chnid,chnid);
    if (debug) printf("Other: %s on res %d comparing with %s on res %d\n",atom[atomnum].atmnam,atom[atomnum].aanum,atmnam,aanum);
    
    if (atom[atomnum].aanum == aanum && !strncmp(atom[atomnum].atmnam,atmnam,4) && atom[atomnum].chnid == chnid)
        return (atomnum);
    for(i=(atomnum<10)?0:atomnum-10; i < atomnum+10 && i < natoms; i++)
        if (atom[i].aanum == aanum && !strncmp(atom[i].atmnam,atmnam,4) && atom[i].chnid == chnid)
            return (i);
    for(i=0; i < natoms ; i++)
        if (atom[i].aanum == aanum && !strncmp(atom[i].atmnam,atmnam,4) && atom[i].chnid == chnid)
            return (i);
    if (debug) printf("Atom %d, %s in residue %d%c not found\n",atomnum,atmnam,aanum, chnid);
    
    return (-99); /* if the atom has not been found */
}    

int find_atom4(char chnid, int aanum, char inscode,char * resnam, char * atmnam, int startindex)
{
    int i;
    int debug=0;
/*    if (startindex<12) debug=0;
    else debug=0;*/
    if (debug) printf("Looking for /%c%04d%c%.3s%.4s, starting at %d\n",chnid,aanum,inscode,resnam,atmnam,startindex);
    
    for(i=startindex;i>=0;i--)
    {
        if (debug) printf("Comparing with /%c%04d%c%.3s%.4s, number %d\n",atom[i].chnid,atom[i].aanum,atom[i].inscode,atom[i].resnam,atom[i].atmnam,i);
        if (atom[i].chnid == chnid && atom[i].aanum == aanum && atom[i].inscode == inscode && !strncmp(atom[i].resnam,resnam,3) && !strncmp(atom[i].atmnam,atmnam,4) )
            return(i);
    }
    for(i=startindex;i<natoms;i++)
    {
        if (atom[i].chnid == chnid && atom[i].aanum == aanum && atom[i].inscode == inscode && !strncmp(atom[i].resnam,resnam,3) && !strncmp(atom[i].atmnam,atmnam,4) )
            return(i);
    }
    return (-99);
}
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


char * atm2molaa(int i)
{
    static char answer[7] ;
    answer[0]='\0';
    
    if(atom[i].chnid!='-')
        sprintf(answer+strlen(answer),"%c",atom[i].chnid);
    sprintf(answer+strlen(answer),"%d",atom[i].aanum);
    if(atom[i].inscode !='-')
        sprintf(answer+strlen(answer),"%c",atom[i].inscode);
    
    return (answer);
}

char * strnremspc(char * string, int n)
{
    static char answer[256]="\0";
    char * i, *j;
    j=answer;
        
    for(i=string;*i&&i<string+n;i++)
        if(*i!=' ')
            *j ++ = *i;
    *j=0;
    return(answer);
}

int hbpabs(int x) /* defined only bc I don't trust the header file */
{
/*    printf("Absolute of %d is %d\n",x,x>0?x:-x);*/
    return ((x > 0) ? x : -x) ;
}


int find_gap(int ri, int rj)
{
    int gap, brk_i;

    if ((ri == -99) || !residue[ri].ca || (rj == -99) || !residue[rj].ca  )
    {/*
        printf("Setting gap to -2\n");
        */
        gap = -2;
    }
    
    else
    {
/*        printf("Proper gap\n");
        */
        if (nbrakes<=0)
        {
          /*  printf("Absolute hbpabs\n");
            */
            gap = hbpabs(ri - rj);
        }
        else
        {
/*            printf("Need to check for brakes\n");
        */    
            for (brk_i = 0; brk_i < nbrakes; brk_i++)
            {
                if (brake[brk_i][0] >= min(ri,rj) &&
                    brake[brk_i][1] <= max(ri,rj))
                {
                /*    printf("Setting gap to -1\n");
                  */  
                    gap = -1;
                    break;
                }
                else
                {
/*                    printf("Called hbpabs\n");
        */            
                    gap = hbpabs(ri - rj);
                }
                
            } /* end for */
            
        }
        
    }
    if (gap < -2 || gap > 2000)
        printf("BUG: Gap between %d & %d = %d\n",ri,rj,gap);
    
    return(gap);
}
/*11. PDBOUT */
/***************************************************************************/

short pdboutflg=0;


void pdbout(char * outfn,char * pdbname,short int filter_flg,short int * filter)
{
     int i,j,atmnum=1,hetflg;
     int debug = 0;
     FILE * outfp = NULL;

     char      calendar_s[128],cappdbname[128],*p;

     time_t    calendar;
     struct tm * calendar_tm_p;

     printf("Outputting PDB format file into %s . . .\n",outfn);

     if(debug) printf("Entering pdbout\n");
     outfp=fmfopen(outfn,"w");
     if (debug) printf("Called the file-opening routine\n");

     time (&calendar); /* whats the time . . .*/
     calendar_tm_p = localtime( &calendar ); /* in a human-readable format */
     strftime(calendar_s,128,"%d-%b-%y\0",calendar_tm_p); /* in PDB format */

     if (debug==1) printf("About to set date to upper case\n");
     for (p=calendar_s; *p; p++)
        *p = toupper(*p);
     if (debug==1) printf ("Just set date to upper case\n");

     fprintf(outfp,"HEADER    %-40s%9s   %4s \n","PDB FORMAT FILE", calendar_s,brcode);
     fprintf(outfp,"REMARK   1 REFERENCE 1\n");
     fprintf(outfp,"REMARK   1  AUTH   I.K.MCDONALD,D.NAYLOR,D.T.JONES,J.M.THORNTON\n");
     fprintf(outfp,"REMARK   1  TITL   HBPLUS - HYDROGEN BOND CALCULATOR\n");
     fprintf(outfp,"REMARK   1  REF    COMPUTER PROGRAM - BIOCHEM. DEPT., UCL,\n");
     fprintf(outfp,"REMARK   1  REF  2 GOWER STREET, WC1E 6BT 1993\n");
     fprintf(outfp,"REMARK   2 ONLY ATOM AND HETATM RECORDS ARE INCLUDED.  HYDROGENS ARE\n");
     fprintf(outfp,"REMARK   2 INCLUDED WHERE POSSIBLE, AS GENERATED BY PROGRAMHBPLUS\n");
     strcpy(cappdbname,pdbname);

     for(p=cappdbname;*p;p++)
        *p = toupper(*p);

     fprintf(outfp,"REMARK   3 ORIGINAL PDB FILE %s.\n",pdbname);
     fprintf(outfp,"REMARK   4 HBPLUS IS FREE TO ALL ACADEMIC USERS.\n");
     fprintf(outfp,"REMARK   4 CONTACT I.K.MCDONALD AT ABOVE ADDRESS OR ON ELECTRONIC\n");
     fprintf(outfp,"REMARK   4 MAIL AT \\MCDONALD@UK.AC.UCL.BSM$ FOR INFORMATION\n");

     if (debug)
     {
        printf("Outputted the header bit\n");
        fclose(outfp);
        outfp=fmfopen(outfn,"a");
     }
     for (hetflg=0; hetflg<2; hetflg++)
        for (i=0;i<natoms;i++)
        {
            if (atom[i].hetflg==hetflg &&
                (!filter_flg || (atom[i].mc_p && filter[ (int)(atom[i].mc_p-mc) ] ) 
#ifdef dduc
                 && (!atom[i].mc_p && atom[i].atmtyp < 4)
#endif
                 )
                )
            {
                {
                    fprintf(outfp, "%-6s", (hetflg)?"HETATM":"ATOM");
                    fprintf(outfp, "%5d", atmnum);
                    fprintf(outfp, " %4s", atom[i].atmnam);
                    fprintf(outfp, "%c", (atom[i].altcode=='-')?' ':atom[i].altcode);
                    fprintf(outfp, "%3s ", atom[i].resnam);
                    fprintf(outfp, "%c", (atom[i].chnid=='-')?' ':atom[i].chnid);
                    fprintf(outfp, "%4d", atom[i].aanum);
                    fprintf(outfp, "%c   ", (atom[i].inscode=='-')?' ':atom[i].inscode);
                    fprintf(outfp, TF, VXYZ(atom[i].p) );
                    if ( atom[i].occ > -99 )
                        fprintf(outfp, "%6.2f",atom[i].occ);
                    else
                        fprintf(outfp, "      " );
                    if (atom[i].b > -99 )
                        fprintf(outfp, "%6.2f", atom[i].b);
                    else
                        fprintf(outfp, "      " );
 /*                 fprintf(outfp, " %4d %2d", atom[i].h_ptr, atom[i].nh);*/
                    if (atom[i].acc > -99.9)
                        fprintf(outfp, "%7.3f", atom[i].acc);
                    else
                        fprintf(outfp, "       ");
                    
                    fprintf(outfp, "\n");
                }
                    if (debug==2)
                    {
                        printf("%-6s", (hetflg)?"HETATM":"ATOM");
                        printf("%5d", atmnum);
                        printf(" %4s", atom[i].atmnam);
                        printf("%c", (atom[i].altcode=='-')?' ':atom[i].altcode);
                        printf("%3s ", atom[i].resnam);
                    printf("%c", (atom[i].chnid=='-')?' ':atom[i].chnid);
                    printf("%4d", atom[i].aanum);
                    printf("%c   ", (atom[i].inscode=='-')?' ':atom[i].inscode);                    printf(TF, VXYZ(atom[i].p) );
                    printf(" %d->%d",atom[i].nh,atom[i].h_ptr);
                    printf("\n");
                }

                if (debug==1) printf("out atom %d ",i);

                atmnum++;
                if (atom[i].nh)
                {
                    for(j=0;j<atom[i].nh;j++)
                        if (h_atm[atom[i].h_ptr+j].typ){
                            if (debug==2)
                            {
                                printf("%-6s", (hetflg)?"HETATM":"ATOM");
                                printf("%5d", atmnum);
                                printf(" %cH%c%c",(atom[i].nh==1)?' ':j+'1',atom[i].atmnam[2], atom[i].atmnam[3]);
                                printf("%c", (atom[i].altcode=='-')?' ':atom[i].altcode);
                                printf("%3s ", atom[i].resnam);
                                printf("%c", (atom[i].chnid=='-')?' ':atom[i].chnid);
                                printf("%4d", atom[i].aanum);
                                printf("%c   ", (atom[i].inscode=='-')?' ':atom[i].inscode);                    printf(TF, VXYZ(h_atm[atom[i].h_ptr+j].p) );
                                printf(" \n");
                            }
                        {
                            fprintf(outfp, "%-6s", (hetflg)?"HETATM":"ATOM");
                            fprintf(outfp, "%5d", atmnum);
                            fprintf(outfp, " %cH%c%c",(atom[i].nh==1)?' ':j+'1',atom[i].atmnam[2], atom[i].atmnam[3]);
                            fprintf(outfp, "%c", (atom[i].altcode=='-')?' ':atom[i].altcode);
                            fprintf(outfp, "%3s ", atom[i].resnam);
                            fprintf(outfp, "%c", (atom[i].chnid=='-')?' ':atom[i].chnid);
                            fprintf(outfp, "%4d", atom[i].aanum);
                            fprintf(outfp, "%c   ", (atom[i].inscode=='-')?' ':atom[i].inscode);
                            fprintf(outfp, TF, VXYZ(h_atm[atom[i].h_ptr+j].p) );
                            if ( atom[i].occ > -99 )
                                fprintf(outfp, "%6.2f",atom[i].occ);
                            else
                                fprintf(outfp, "      " );
                            if (atom[i].b > -99 )
                                fprintf(outfp, "%6.2f", atom[i].b);
                            else
                                fprintf(outfp, "      " );
                            if (atom[i].acc > -99.9)
                                fprintf(outfp, "%7.3f", atom[i].acc);
                            else
                                fprintf(outfp, "       ");
                            
                            
                            fprintf(outfp, " \n");
                        }
                            if (debug ==2) printf("out h %d \n",j);
                            
                            /* output a line of brookhaven */
                            /* the H atom branch numbering is not unique - Hs
                               attached to two zeta atoms would not be unique */
                            atmnum++;
                        } /* end of j (Hydrogen atom) loop */
                } /* end of hydrogen block */
                } /* end of "if it is right one from atom or hetatm" */
            
        } /* end of i loop that goes through atom array once */
     /*end of hetflg loop that executes i loop twice and has no brackets*/
     fclose(outfp);
     return;
 }

/*****************************************************************************/
/*===========================================================================*/

void newparsefn(char * fname,char * rootnam,int rootlen)

/* Takes the filename in FNAME and strips off any leading path together with
   any trailing extension. The result is returned in ROOTNAM. For example,
   the string "/usr/fred/filename.ext" is returned as "filename".  
   ROOTLEN should specify the actual size of the string array ROOTNAM. */

{
    char  *p, *dotpos, *slashpos;
    int    i, len, rlen, dbgtemp;

    dbgtemp = debug;
    debug = 0;
    if (debug == 1)
    {
        fprintf(dbgfp, "FNAME   passed to newparsefn is '%s' with length %d\n", 
                        fname, strlen(fname)); 
        fprintf(dbgfp, "ROOTNAM passed to newparsefn is '%s' with length %d\n", 
                        rootnam, strlen(rootnam)); 
    }
    dotpos   = 0;
    slashpos = 0;
    len      = strlen(fname);

    /* Begin by blanking out current/previous contents of ROOTNAM */

    for (i=0; i<rootlen; i++)
    {
        p = rootnam + i;
        *p = '\0';
    }
    if (debug == 1)
        fprintf(dbgfp, "ROOTNAM after blanking is '%s' with length %d\n",
                        rootnam, strlen(rootnam));

    for (i=1; i<=len; i++)
    {
        p = fname + len - i;
        if (*p == '/')
        {
            if (slashpos == 0)
                slashpos = p;
            break;               /* once a "/" is found stop looking for "."s */
        }
        if (*p == '.')
        {
            if (dotpos == 0)
                dotpos = p;
        }
    }

    if (slashpos == 0)             /* no leading path */
        slashpos = fname - 1;

    if (dotpos == 0)               /* no trailing extension */
        dotpos = fname + len;

    rlen = dotpos - slashpos - 1;
    strncpy(rootnam, slashpos+1, rlen);
  
    if (debug == 1)
        fprintf(dbgfp, "ROOTNAM being returned is  '%s'\n", rootnam); 
    debug = dbgtemp;

}
/*============================================================================*/

void parsefn(char * fname,char * rootnam)

{ 

    /* takes the filename in fname, strips it of any leading path together with
       any trailing extension.  The result is returned in rootnam. eg
       /home/bsm/mcdonald/prog/test.pdb becomes test. */
       /*NOTE: including the full stop */
    /* I think this is a DN routine that I edited to cope with names without
        suffixes */

/* Amendment. RAL 13 Aug 2009 --> */
   char slash;
/* <-- RAL 13 Aug 2009 */

    char *p, *dotpos, *slashpos;
    int len;
    
    int debug=0;

    dotpos = 0;
    slashpos = 0;
    len = strlen(fname);
    
    if (debug)
        printf("DEBUG: fname is %s\n",fname);
    
    p = fname + strlen(fname);
/*    while (*--p != '.' && p != '/' && p > fname);*/
    
/*        while (isalnum(*--p))
        if (p == fname)
        {
            p--;
            break;
        }*/
/* Amendment. RAL 13 Aug 2009 --> */
/*        while (*--p != '/' && p >= fname); */
        while (*p != '\\' && *--p != '/' && p >= fname);
/* <-- RAL 13 Aug 2009 */
    strcpy(rootnam, p + 1);

    
    p = rootnam + strlen(rootnam);
 /*   while (*--p != '.');*/
    
    while (*--p != '.' && p >= rootnam);
    if(*(p) == '.' )
        *(p +1) = '\0';
    else
    {
        strcpy(rootnam + strlen(rootnam), ".");
        /* I haven't use *address = char, because I also want to append a '\0' */
    }
    
    if (debug)
        printf("DEBUG: root is %s\n",rootnam);   
    return;
    
}

/* Amendment. RAL 3 Jul 2012 --> */
/***********************************************************************

string_truncate  -  Truncate the given string at the last non-blank
                    character, or at its maximum length

***********************************************************************/

int string_truncate(char *string,int max_length)
{
  char ch;

  int iend, ipos, len;

  /* Initialise variables */
  iend = -1;
  ch = string[0];
  for (ipos = 0; ipos < max_length && ch != '\0'; ipos++)
    {
      /* Get the current character and check for end of string */
      ch = string[ipos];
      if (ch == '\n' || ch == 13)
        ch = '\0';
      if (ch != '\0')
        {
          /* If not a space, then mark end of string */
          if (ch != ' ')
            iend = ipos;
        }
    }

  /* End the string after the last non-blank character */
  if (iend + 1 < max_length)
    {
      string[iend + 1] = '\0';
      len = iend + 1;
    }
  else
    {
      string[max_length - 1] = '\0';
      len = max_length - 1;
    }

  /* Return the string length */
  return(len);
}
/* <-- RAL 3 Jul 2012 */
