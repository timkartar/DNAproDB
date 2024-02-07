/* HBPLUS - Hydrogen Bond Calculator v 3.2 */
/* Copyright I K Mcdonald and J M Thornton 1993
             All Rights Reserved
   Copies are free to Academic Users - contact mcdonald@uk.ac.ucl.bsm.bioc.bsm

 The publication of research using the Software must reference "McDonald
 IK, Naylor DN, Jones DT & Thornton J M (1993), 'HBPLUS', Computer Program,
 Department of Biochemistry and Molecular Biology, University College
 London." or successor references as defined by the authors.

 Unless informed otherwise, you should be an academic user and have sent a
 signed confidentiality agreement to the authors (address at the end of the
 source code hbplus.h).  If you do not have a copy of HBPLUS and would like
 to receive one (free to academic users), please detach the confidentiality
 agreement from the end of this document, sign it, and send to the address
 given.  Please allow other people in your department to use your copy of
 HBPLUS, but do not allow them to make their own copy. */

/* Contents */
/* 1. Copyright - all**/
/* 2. Version Log - all**/
/* 3. Tokens as abbreviations**/

/* 4. Tokens as flags and customisables */
/* 5. Globals as flags and customisables */
/* 6. Globals as name arrays */
/* 8. Vector/Matrix Routines */
/* 8.5 Globals as arrays */
/* 8.7Minor useful functions */

/* 9. Hydrogen Position Calculation Subroutines **/

/*10. Hydrogen Position Calculation Routine **/

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

   Version 2.25    IM 21st January 1994
                   Separating into files - creation of hbp_findh.c

   Version 3.0     IM 15th June 1994

*/

#include <string.h>
#include <ctype.h>
#include <math.h>
#include "hbplus.h"


/****************************************************************************/
/**  9. Hydrogen Position Calculation Subroutines                          **/
/****************************************************************************/

/* Struct loci is defined in the mathematical routines */
/* HLOCI - see hbplus */


/* this is a table of all the atoms needed for predicting the positions of
   hydrogens.  This time round, I'll give them all as letters.  Numbers,
   perhaps later.  I will need donor, dd, ddd and possibly dd1 or dd2.  The
   kinds of "frames" I will model are:
	 planar 1H (eg peptide NH), d, dd1, dd2
	 planar 2H (eg amide)       d, dd, ddd1, ddd2
	  also applies to Tyr
	 hydroxyl H                 d, dd,
	 tetrahedral 3H (eg Lys)    d, dd, ddd
 */

enum HTYPE 
{
    planH=1, planH2, planOH, tetrOH, tetrH3
    }
;

#define TOTDONATM 32

const struct posH
{
    int aac,d,typ;
    int dd,ddd,d1,d2; /* d1 and d2 hold any pair of "frame" atoms that
	       	   branch out at either side */
    float dh, ang, sinang, cosang;
} posH[TOTDONATM] /* (how to find the) position (of the) Hydrogen */ =
{
    { -99,  0, planOH,-99,  0, -1,  2,  1.00 ,   4, 0.069756, 0.997564},
/* Backbone NH ele 0 - angles and lengths from Pauling et al 1951 PNAS V37
   P205 */
    { -99,  0, tetrH3,  1,  4,  0,  2,  1.00 , 110, 0.93969, -0.34202}, 
/* Backbone Terminus ele 1*/
    { Cys,  5, tetrOH,  4,  1,  0,  0,  1.33 ,  96, 0.99452, -0.10453},
    { Cyh,  5, tetrOH,  4,  1,  0,  0,  1.33 ,  96, 0.99452, -0.10453},
    { His,  6, planH,   0,  0,  5,  8,  1.00 },
    { His,  9, planH,   0,  0,  7,  8,  1.00/*?*/},
/*      " N   CA  C   O   CB  CG  ND1 CD2 CE1 NE2 OXT", */
    { His,  7, planH,   0,  0,  5,  9,  1.00}, /* if his swappable */
    { His,  8, planH,   0,  0,  9,  6,  1.00}, /* if his swappable */
    { Lys,  8, tetrH3,  7,  6,  0,  0,  1.014, 112, 0.92718, -0.37461},
    { Asn,  7, planH2,  5,  0,  4,  6,  1.00,  120, 0.86603, -0.50000},
    { Gln,  8, planH2,  6,  0,  5,  7,  1.00,  120, 0.86603, -0.50000 },
    { Asn,  6, planH2,  5,  0,  4,  7,  1.00,  120, 0.86603, -0.50000}, /*sw*/
    { Gln,  7, planH2,  6,  0,  5,  8,  1.00,  120, 0.86603, -0.50000 }, /*sw*/
    { Arg,  7, planH,   0,  0,  6,  8,  1.00},
    { Arg,  9, planH2,  8,  0,  7, 10,  1.00,   120, 0.86603, -0.50000},
    { Arg, 10, planH2,  8,  0,  7,  9,  1.00,   120, 0.86603, -0.50000},
    { Ser,  5, tetrOH,  4,  1,  0,  0,  1.00,  110, 0.93969, -0.34202},
    { Thr,  5, tetrOH,  4,  1,  0,  0,  1.00,  110, 0.93969, -0.34202},
    { Trp,  8, planH,   0,  0,  6,  9,  1.00},
    { Tyr, 11, planOH, 10,  0,  9,  8,  1.00,  110, 0.93969, -0.34202},
    { Mpr,  5, tetrOH,  4,  1,  0,  0,  1.33,   96, 0.99452, -0.10453},
    { Lym,  8, tetrH3,  7,  6,  0,  0,  1.014, 112, 0.92718, -0.37461},
    { Lov, 13, tetrOH, 12, 11,  0,  0,  1.00,  110, 0.93969, -0.34202},
    { Sta,  9, tetrOH,  8,  1,  0,  0,  1.00,  110, 0.93969, -0.34202},
    { Cal, 12, tetrOH, 11,  1,  0,  0,  1.00,  110, 0.93969, -0.34202},
    { Ahs, 12, tetrOH, 11,  1,  0,  0,  1.00,  110, 0.93969, -0.34202},
    { Chs, 12, tetrOH, 11,  1,  0,  0,  1.00,  110, 0.93969, -0.34202},
    { Pca,  7, tetrOH,  6,  5,  0,  0,  1.00,  110, 0.93969, -0.34202},
    { Asx,  6, planH2,  5,  0,  4,  6,  1.00,   120, 0.86603, -0.50000},
    { Asx,  7, planH2,  5,  0,  4,  7,  1.00,   120, 0.86603, -0.50000},
    { Glx,  8, planH2,  6,  0,  5,  7,  1.00,  120, 0.86603, -0.50000 },
    { Glx,  7, planH2,  6,  0,  5,  8,  1.00,  120, 0.86603, -0.50000 } /*sw*/
};

/* return the position of a single proton attached to a sp2 atom *
 * in other words, on the main chain, Trp, His ND1, Arg NE Nitrogens*
 * where the donor atom is attached to two other atoms */


struct vect position_of_single_dh(struct vect dd1, struct vect donor, struct vect dd2, struct posH donrow)
{
    struct vect dd1_donor,dd2_donor,mean, h;
    dd1_donor=unit_vector( (from_a_to_b(dd1,donor)));
    dd2_donor=unit_vector( (from_a_to_b(dd2,donor)));
    mean= unit_vector(float_times_vect( 0.5, vector_plus(dd1_donor,dd2_donor)));
    h= vector_plus (donor, float_times_vect(donrow.dh, mean));
    return (h);
}

void position_of_planar_dh2_2
    (struct vect h[2],struct vect donor,struct vect dd, struct vect ddd1, struct vect ddd2, struct posH donrow )
{
    struct vect crossbar;
    crossbar=unit_vector( from_a_to_b(ddd1,ddd2) );
    crossbar=unit_vector(crossbar);
    h[0]=vector_plus(donor, float_times_vect(donrow.dh * cos30 ,crossbar));
    h[0]=vector_plus(h[0], float_times_vect(donrow.dh*0.5, unit_vector( from_a_to_b(dd,donor)) ) );
    h[1]=vector_plus(donor, float_times_vect(-donrow.dh* cos30 ,crossbar));
    h[1]=vector_plus(h[1], float_times_vect(donrow.dh*0.5, unit_vector( from_a_to_b(dd,donor)) ) );
    return;
	       }

void position_of_tetra_dh3(struct vect position[3],struct vect donor,struct vect dd, struct vect ddd, struct posH donrow)
		   /* used for Lys and Thr, Ser, Cys, but not Tyr */
{
    struct vect newx,newy,newz;
    struct vect h1, h2, h3; 
    struct mat3x3 new_axs;
    h1 = vector_of(-donrow.cosang, donrow.sinang, 0.0);
    h2 = vector_of(h1.x, -donrow.sinang*.5, donrow.sinang*0.86603);
    h3 = vector_of( h2.x, h2.y,-h2.z );
    newx=unit_vector( to(dd,donor) );
    newz=unit_vector(cross_product( newx, to(ddd,dd) ) );
    newy=cross_product( newz, newx );
    new_axs=axs_rot(newx, newy, newz);
    h1= mlt_mtx_vct(new_axs,h1);
    h2= mlt_mtx_vct(new_axs,h2);
    h3= mlt_mtx_vct(new_axs,h3);
    position[0]=vector_plus(donor,float_times_vect(donrow.dh,h1));
    position[1]=vector_plus(donor,float_times_vect(donrow.dh,h2));
    position[2]=vector_plus(donor,float_times_vect(donrow.dh,h3));
    return;
}    

struct vect position_of_tetra_oh(struct vect donor, struct vect dd, struct vect ddd, struct posH donrow)
{
    struct vect newx,newy,newz;
    struct vect h1 ;
    struct mat3x3 new_axs;
    h1 = vector_of( -donrow.cosang, donrow.sinang, 0.000 );
    newx=unit_vector( to(dd,donor) );
    newz=unit_vector(cross_product( newx, to(ddd,dd) ) );
    newy=cross_product( newz, newx );
    new_axs=axs_rot(newx, newy, newz);
    h1= float_times_vect(donrow.dh,mlt_mtx_vct(new_axs,h1));
    return (vector_plus(donor,h1) );    
}

struct loci loci_of_tetra_oh(struct vect donor, struct vect dd, struct vect ddd, struct posH donrow)
{
    struct loci result;
    result.typ = circle;
    result.p = position_of_tetra_oh(donor, dd, ddd, donrow);
    result.a = float_times_vect(donrow.dh*-donrow.cosang, unit_vector( to(dd,donor) ) );
    result.a = vector_plus( result.a, donor);
    return (result);
}

void position_of_planar_dh2(struct vect h[2], struct vect donor, struct vect dd, struct vect ddd1, struct vect ddd2 /*not used*/,struct posH donrow)
{
    struct loci loci;
    loci=loci_of_tetra_oh(donor,dd, ddd1, donrow);
    h[0]=loci.p;
    h[1]=vector_plus(loci.p, float_times_vect(2,to(loci.p,loci.a) ) );
    return;
}

struct vect position_of_planar_oh(struct vect donor, struct vect dd, struct vect ddd1, struct vect ddd2, struct posH donrow)
{
    struct vect oh[2];
    position_of_planar_dh2(oh,donor,dd,ddd1,ddd2,donrow);
    return (oh[0]);
}



/****************************************************************************/
/*10. Hydrogen Position Calculation Routine*/
/* The routine to position Hydrogens - Ian McDonald's bit ! */

/* How to generate a donor type off hand                                      
                                                                              
                       N                   O                                  
   .ncon                                                                      
        0     1     2     3          0    1    2                              
       /      |     |      \        /     |     \                            
      /       |     |       X    none     |      X                            
 none         |     planH                 |                                   
              |                       [dd].ncon                               
           [dd].ncon                [must be C]                   
          [must be C]                                                         
                                 2     3      1     4                         
 1    2    3    4                              \    |                         
 |               \                              \   |                         
WARN  check       WARN           check the       \  |                         
      angle                         angle         \ |                         
      |    \                                       \|                         
      |     \                        / \            |                         
      |      \                      /   \           |                         
      |       |                    /     \          |                         
      |       |                   |       \         /                         
      |       |                   |        \       /                          
      |       |                   |         \     /                           
   planH2  tetrH3               planOH    tetrOH_/                            

*/    

void find_h()
{
    struct vect pos_buf[3];
    int i,j,c,ca; /*n has been deleted because atom[i] is used instead */
    int d, dd, ddd, d1, d2, dd1, dd2, ddd1, ddd2;
    int debug=0;
    short tmpatmtyp; /* used in place of calling the variable direct,
			  so that swapped atoms can be dealt with properly */
    char atomidbuf[18];
    struct posH donrow;
    
    
    if (debug==1) printf("\nPositioning Hs on all %d atoms\n",natoms);
    printf("\nAdding Polar Hydrogens.\n");
    
    
    for(i=0; i<natoms; i++)
    {
	if (debug==2) printf(" %d ",i);
	if (debug && atom[i].aanum==2)
	    printf("%d %c%d %s %s nh %d typ %d ",i,SHORTID(i),atom[i].nh,h_atm[atom[i].h_ptr].typ);
	
	if (atom[i].nh && !(h_atm[ atom[i].h_ptr ].typ) )
	    /* does it have any hydrogens on it anyway? */
	{
	    /*if (i<5) debug=1; else debug=0; */
	    
	    if (debug) printf("\n %d %s ",i,atomid(i,atomidbuf) );
	    
	    
	    if ((atom[i].atmtyp==0) /* main chain NH*/ && atom[i].nh == 1)
	    {
		
		if (debug) 
		    printf("Backbone H on res %d",atom[i].aanum);
		/* set for normal HBPLUS NH positions.  Unset for K&S positions */
		if( !kshflg )
		{
		    if (atom[i].ncon==2) /* the new technique */
		    {
			c=icon[i][0];
			ca=icon[i][1];
			if ( strncmp(" C  ",atom[c].atmnam,4) )
			{
			    c=icon[i][1];
			    ca=icon[i][0];
			}
			
		    }
		    else
		    {
			c=find_atom(atom[i].caindex-1,atom[i].chnid,2, i+10);
			ca=find_atom(atom[i].caindex,atom[i].chnid, 1, i+1);
			if (ca != -99 && c != -99)
			{
			    if (debug && atom[i].aanum==2) printf ("Got C=%c%d %s %s and CA=%c%d %s %s",atom[c].chnid,atom[c].aanum,atom[c].resnam,atom[c].atmnam,atom[ca].chnid,atom[ca].aanum,atom[ca].resnam,atom[ca].atmnam);
			}
		    }
		    if (ca == -99 || c == -99)
			printf("WARNING: Cannot locate both donor antecedents of backbone %s.\n", atomid(i,atomidbuf));
		    else
		    {
			pos_buf[0]=position_of_single_dh(atom[ca].p,atom[i].p,atom[c].p, /*donrow*/ posH[0]);
			/* this puts the 'bisecting' position in pos_buf[0] */
			
			h_atm[atom[i].h_ptr].p= position_of_planar_oh( 
								      atom[i].p, /*donor*/
								      /*artificial dd*/      pos_buf[0]
								      , atom[c].p, atom[ca].p, posH[0]);
			/* this is unusual, because the H is actually 4 degrees lobsided (Pauling et
			   al 1951 (yes, 1951.  It isn't a misprint), so I'm calling an unexpected 
			   routine to put it on a stalk that is created from the midpoint of the ca
			   and c atoms to the N atom, similarly to SER OG or THR OG1 */
			h_atm[atom[i].h_ptr].typ= fixed;
		    }
		}/* pairs if !kshflg */
		else
		{ /* Kabsch-Sander Hydrogen Positions, parralel to CO */
		    h_atm[atom[i].h_ptr].p = 
			vector_plus ( atom[i].p , unit_vector( from_a_to_b 
							      ( *(residue[ atom[i].caindex-1 ].o), *( residue[ atom[i].caindex-1 ].c) ) ) );
		    h_atm[atom[i].h_ptr].typ = fixed;
		}
		
	    }
	    else /* ie will be in the posH table */
	    {
		donrow.d= -99;
		donrow.typ=0;
		
		for(j=2;j<TOTDONATM;j++)
		{
		    if (atom[i].aacode == posH[j].aac && atom[i].atmtyp == posH[j].d )
		    {
			donrow=posH[j];
		    }
		}
		if (atom[i].nh == 3 && atom[i].atmtyp ==  0)
		    donrow=posH[1];/*Backbone Terminus*/
		
		if(debug) printf(" donrow = %d ",donrow.aac);
		
		d=i;
		dd=dd1=dd2=ddd=ddd1=ddd2= -99;
		
		switch (atom[i].ncon) /*Fix Antecedents*/
		{
		case 0:
		    continue;
		case 1:
		    dd=icon[i][0];
		    if (dd < 0 || dd > natoms)
			printf("BUG: %d in icon array for %s.\n",dd,atomid(i,atomidbuf));
		    switch ( atom[dd].ncon )
		    {
		    case 1:
			printf("WARNING: %s and",atomid(dd,atomidbuf));
			printf(" %s form 1 molecule.\n",atomid(i,atomidbuf));
			continue;
		    case 2:
			if (icon[dd][0]==d)
			    ddd=icon[dd][1];
			else
			    ddd=icon[dd][0];
			break;
		    case 3:
			ddd1= icon[dd][0];
			ddd2= icon[dd][1];
			if (ddd1 == d)
			    ddd1=icon[dd][2];
			if (ddd2 == d)
			    ddd2=icon[dd][2];
			break;
		    case 4:
			printf("WARNING: 4 bonds to donor antecedent %s.\n",atomid(i,atomidbuf));
			continue;
		    default:
			printf("WARNING: Unexpected number of bonds to dd %s.\n",atomid(dd,atomidbuf));
			continue;
		    }
		    break;
		case 2:
		    if ( atom[i].atmnam[1]=='O' || atom[i].atmnam[1]=='S' )  /*ie otherwise bonded*/
			continue;
		    dd1=icon[i][0];
		    dd2=icon[i][1];
		    break;
		default:
		    printf("WARNING: %s forms %d covalent bonds.\n",atomid(d,atomidbuf),atom[d].ncon);
		    continue;
		}/*Fix Antecedents*/
		
		/*fix ddd, if it has not already been fixed*/
		if (ddd== -99)
		{
		    if (ddd1>-99 && ddd2 > -99)
			if (atom[ddd1].ncon != atom[ddd2].ncon)
			    ddd = (atom[ddd1].ncon > atom[ddd2].ncon)?ddd1:ddd2;
			else if (atom[ddd1].atmtyp != atom[ddd2].atmtyp)
			    ddd = (atom[ddd1].atmtyp > atom[ddd2].atmtyp)?ddd1:ddd2;
			else ddd = ddd1;
		    else
			if (ddd1== -99)
			    ddd=ddd2;
			else
			    ddd=ddd1;
		}
		/* this could leave ddd unset if both ddd1 and ddd2 were also unset*/
		
		
			    
		tmpatmtyp = atom[i].atmnam[1] ;
		
		if (exchangeflg &&
		    ( (atom[i].aacode == His && atom[i].atmnam[1] == 'C') ||
		      ( (atom[i].aacode == Gln || atom[i].aacode == Asn) &&
		         atom[i].atmnam[1] == 'O' &&
		         atom[i].atmnam[3] == '1'                       ) ) ) 
		    tmpatmtyp = 'N';

/*Version 2.28*/
		if (isdigit(tmpatmtyp))
		    tmpatmtyp = atom[i].atmnam[0];
		
		switch( tmpatmtyp )/*Find the donor role*/
		{
		case 'N':
		    donrow.dh=1.00;
		    switch (atom[i].ncon)
		    {
		    case 1:
			if (atom[dd].ncon != 2 && atom[dd].ncon !=3)
			{
			    printf("WARNING: Unexpected ncon for %s.\n",atomid(dd,atomidbuf));
			    break;
			}
			if (atom[i].nh==2)
			{
			    donrow.ang=120;
			    donrow.sinang=0.86603;
			    donrow.cosang= -0.50000;
			    donrow.typ=planH2;
			}
			if (atom[i].nh==3)
			{
			    donrow.ang=112;
			    donrow.sinang=0.92718;
			    donrow.cosang= -0.37461;
			    donrow.dh=1.014;
			    donrow.typ=tetrH3;
			}
			if (atom[i].nh==1)
			    printf("WARNING: Only one H on %s.\n",atomid(d,atomidbuf));
			
			break;
		    case 2:
			donrow.typ=planH;
			break;
		    case 3:
			continue;
		    default:
			printf("WARNING: %s forms %d covalent bonds.\n",atomid(i,atomidbuf),atom[i].ncon);
			break;
		    }
		    break;
		    
		case 'S':
		case 'O':
		    if (atom[i].ncon!=1)
			continue;/*can only donate if ncon < 2*/
		    if (atom[dd].ncon ==1 || atom[dd].ncon ==4)
		    {
			printf("WARNING: %s forms %d covalent bonds.\n",atomid(dd,atomidbuf),atom[dd].ncon);
			/*continue;*/
		    }
		    if (atom[d].nh !=1)
		    {
			printf("WARNING: %s has %d hydrogens.\n",atomid(d,atomidbuf),atom[d].nh);
			continue;
		    }
		    
		    if ( atom[i].atmnam[1] == 'O' )
		    {
			donrow.dh=1.00;
			donrow.ang=110;
			donrow.sinang=0.93969;
			donrow.cosang= -0.34202;
		    }
		    else
		    {
			donrow.dh=1.33;
			donrow.ang=96;
			donrow.sinang=0.99452;
			donrow.cosang= -0.10453;
		    }
		    if (isplanar(dd) && atom[dd].ncon>=3)
			donrow.typ=planOH;
		    else
			donrow.typ=tetrOH;
		    break;

		default:
		    printf("WARNING: %s not recognised as N, O or S\n",atomid(i,atomidbuf));
		    printf("WARNING: This atom is listed as a donor but it's hydrogen positions cannot be\n         calculated.  This may cause BUG messages.\n");
		    break;
		}
		
		    
		if (donrow.typ == -99)
		{
		    if (atom[i].aacode < TOTNAA)
		      printf("WARNING: Unforseen donor type for %s in find_h.\n", atomid(i,atomidbuf));
			/*printf("Please mail mcdonald@uk.ac.ucl.biochemistry.bsm\n");*/
		    

		} /* if donrow.d == -99 */

		else
		    switch ( donrow.typ )
		    {
			
		    case planH: 
		    {
			
			if (debug) printf(" planH ");
			d = i;
			/*d1= find_atom(atom[i].caindex, atom[i].chnid, donrow.d1, i);
			d2= find_atom(atom[i].caindex, atom[i].chnid, donrow.d2, i);*/
			d1 = dd1;
			d2 = dd2;
			
			if (debug){
			    if (d1!= -99) printf ("d1 %c%d %s %s ",SHORTID(d1));
			    if (d2!= -99) printf ("d2 %c%d %s %s ",SHORTID(d2));
			}
			
			if ( d== -99 || d1== -99 || d2== -99 )
			{
			    printf("\nWARNING: Cannot locate all donor antecedents of sp2 1H 2DD %s\n", atomid(i,atomidbuf));
			}
			
			else
			{    
			    h_atm[atom[i].h_ptr].p=position_of_single_dh
				(atom[d1].p,atom[d].p,atom[d2].p, donrow);
			    h_atm[atom[i].h_ptr].typ=fixed;
			}
		    }
			break;
		    case planH2: 
		    {
			if(debug) printf("\n planH2 ");
			
			/*dd= find_atom(atom[i].caindex, atom[i].chnid, donrow.dd, i);
			d1= find_atom(atom[i].caindex, atom[i].chnid, donrow.d1, i);
			d2= find_atom(atom[i].caindex, atom[i].chnid, donrow.d2, i);*/
			d1=ddd1;
			d2=ddd2;
			
			if (debug){
			    if (dd != -99) printf(" dd %c%d %s %s ",SHORTID(dd));
			    if (d1 != -99) printf(" d1 %c%d %s %s ",SHORTID(d1));
			    if (d2 != -99) printf(" d2 %c%d %s %s ",SHORTID(d2));
			}
			
			if ( i== -99 || d1== -99 || d2== -99 || dd== -99)
			    printf("\nWARNING: Cannot locate all donor antecedents of sp2 2H 1DD %s\n", atomid(i,atomidbuf));
			else
			{
			    position_of_planar_dh2(pos_buf,
						   atom[i].p,atom[dd].p,atom[d1].p,atom[d2].p,donrow);
			    for(j=0;j<2;j++)
			    {
				h_atm[atom[i].h_ptr+j].p=pos_buf[j];
				h_atm[atom[i].h_ptr+j].typ=fixed;
			    }
			}
		    } 
			break;
		    case planOH:
			d = i;
			/*dd= find_atom(atom[i].caindex, atom[i].chnid, donrow.dd, i);
			d1= find_atom(atom[i].caindex, atom[i].chnid, donrow.d1, i);
			d2= find_atom(atom[i].caindex, atom[i].chnid, donrow.d2, i);*/
			d1=ddd1;
			d2=ddd2;
			
			if (debug)
			{
			    if(dd != -99) printf("dd %c%d %s %s ",SHORTID(dd));
			    if(d1 != -99) printf("d1 %c%d %s %s ",SHORTID(d1));
			    if(d2 != -99) printf("d2 %c%d %s %s ",SHORTID(d2));
			}
			
			if ( d== -99 || d1== -99 || d2== -99 || dd== -99)
			    printf("\nWARNING: Cannot locate all donor antecedents of sp2 1H 1DD %s\n", atomid(i,atomidbuf));
			else
			{
			    position_of_planar_dh2(pos_buf, atom[d].p,
						   atom[dd].p,atom[d1].p,atom[d2].p, donrow);
			    h_atm[atom[i].h_ptr].typ=alternatives;
			    h_atm[atom[i].h_ptr].p=pos_buf[0];
			    h_atm[atom[i].h_ptr].a=pos_buf[1];
			}
			break;
			
		    case tetrOH:
			/*dd= find_atom(atom[i].caindex, atom[i].chnid, donrow.dd, i);
			ddd= find_atom(atom[i].caindex, atom[i].chnid, donrow.ddd, i);*/
			if (debug)
			{
			    printf(" tetrOH ");
			    if (dd!= -99) printf("dd %c%d %s %s ",SHORTID(dd));
			    if (ddd!= -99)printf("ddd %c%d %s %s ",SHORTID(ddd));
			}
			
			if ( d== -99 || dd== -99 || ddd== -99)
			    printf("\nWARNING: Cannot locate all donor antecedents of sp3 1H 1DD %s\n", atomid(i,atomidbuf));
			else
			{
			    if (debug)
				printf("tetrOH on atom %d\n",d);
			    h_atm[atom[i].h_ptr]=loci_of_tetra_oh
				(atom[i].p,atom[dd].p,atom[ddd].p, donrow);
			}
			break;
			
		    case tetrH3:
			d = i;
			/*dd= find_atom(atom[i].caindex, atom[i].chnid, donrow.dd, i);
			ddd= find_atom(atom[i].caindex, atom[i].chnid, donrow.ddd, i);*/
			if (ddd== -99 && donrow.d2 > 0)
			    ddd = find_atom(atom[i].caindex, atom[i].chnid, donrow.d2, i);
			if (debug)
			{
			    printf(" tetrH3 ");
			    if (dd!= -99) printf("dd %c%d %s %s ",SHORTID(dd));
			    if (ddd!= -99)printf("ddd %c%d %s %s ",SHORTID(ddd));
			}
			if ( dd== -99 || ddd== -99)
			{
			    if (dd== -99)
				printf("\nWARNING: Cannot locate %.4s donor antecedant of terminyl %s\n", necatm[atom[i].aacode]+4*donrow.dd,atomid(i,atomidbuf));
			    if (ddd== -99)
				printf("\nWARNING: Cannot locate %.4s donor antecedant of terminyl %s\n", necatm[atom[i].aacode]+4*donrow.dd,atomid(i,atomidbuf));
			    
			}
			else
			{
			    if (debug)
				printf("tetrH3 on atom %d\n",d);
			    position_of_tetra_dh3(pos_buf,
						  atom[d].p,atom[dd].p,atom[ddd].p,donrow);
			    for(j=0;j<3;j++)
			    {
				h_atm[atom[i].h_ptr+j].p=pos_buf[j];
				h_atm[atom[i].h_ptr+j].typ=fixed;
			    }
			}
			break;
		    default:
			printf("\nBUG: Unforseen posh[].typ for %s in find_h.\n", atomid(i,atomidbuf));
			printf("Please mail mcdonald@uk.ac.ucl.biochemistry.bsm\n");
		    }/* end switch*/

	    }/*end the "if not amide H" */
	    
	}/* "if it is a donor with unset H" block ends */
	if(debug && i>natoms-5)printf("%d of %d, ",i, natoms);
	
    }/* i loop */
    if (debug) printf("Found Hs from the last of all %d atoms", natoms);
    return;
}

