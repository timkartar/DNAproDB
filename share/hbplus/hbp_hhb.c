/* HBPLUS - Hydrogen Bond Calculator v 3.2 */
/* Hydrogen Bond Calculation Routines */
/* Copyright I K Mcdonald, D T Jones, D Naylor and J M Thornton 1993
             All Rights Reserved
   Copies are free to Academic Users - contact mcdonald@uk.ac.ucl.bsm.bioc.bsm

 The publication of research using the Software must reference "McDonald
 IK, Naylor DN, Jones DT & Thornton J M (1993), 'HBPLUS', Computer Program,
 Department of Biochemistry and Molecular Biology, University College
 London." or successor references as defined by the authors.

 Unless informed otherwise, you should be an academic user and have sent a
 signed confidentiality agreement to the authors (address at the end of the
 hbplus.h source code).  If you do not have a copy of HBPLUS and would like
 to receive one (free to academic users), please detach the confidentiality
 agreement from the end of this document, sign it, and send to the address
 given.  Please allow other people in your department to use your copy of
 HBPLUS, but do not allow them to make their own copy. */

/* Contents */
/* 1. Copyright - all*/
/* 2. Version Log - all*/
/* 3. Tokens as abbreviations */

/* 4. Tokens as flags and customisables **/

/* 5. Globals as flags and customisables **/

/* 6. Globals as name arrays **/

/* 8. Vector/Matrix Routines **/

/* 8.5 Globals as arrays **/

/* 8.7Minor useful functions */
/* 9. Hydrogen Position Calculation Subroutines */
/*10. Hydrogen Position Calculation Routine */

/*11. PDBOUT **/

/*12. Connectivity Calculation Subroutines */
/*13. Connectivity Calculation Routine */
/*13.2 Domain calculation Routines */
/*13.5 Adding new residues to the selection */

/*14. find_hb **/

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
		   by compilers that did not have BSMSOL1 predefined (ie
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
                   Tell the bsmsol1 version to look in ~/data/pdb for 
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
   
   Version 2.23    IM 4th January 1993
                   Separating code into >1 file.

   Version 2.25    IM 21st January 1993
                   Further Separation
*/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "hbplus.h"

/* 3. Tokens as abbreviations */

/* see hbplus.h */

/* 4. Tokens as flags and customisables */
/*********************************************************************/

float   sinMIN_DHA = 1;
float    cosMIN_DHA = 0; /* this means the dot product of the two vectors
		       would come above zero if it below this minimum*/
double    mindha = 90;    /*angle */
    
float    sinMIN_DAAA = 1;
float    cosMIN_DAAA = 0; /* this means the dot product of the two vectors
                        would come above zero if it below this minimum*/
double    mindaaa = 90;   /*angle */
    
float    sinMIN_HAAA = 1;
float    cosMIN_HAAA = 0; /* this means the dot product of the two vectors
                        would come above zero if it below this minimum*/
double    minhaaa = 90;   /*angle */
    
    
float    sinMAX_DAAX = 0.34202014;
float    cosMAX_DAAX = 0.93969262;
double    maxdaax = 20;
    
float    sinMAX_HAAX = 0.34202014;
float    cosMAX_HAAX = 0.93969262;
double    maxhaax = 20;
    
float    MAX_HA = 2.5;
float    HBDIST = 3.9;
    
float    CAWARN = 12.0;



/* 5. Globals as flags and customisables */

#ifdef BSM
short longoutflg=1; /* =0 *.hb2 format =1 *.hhb format */
#else
short longoutflg=0;
#endif /* =0 for hb2/nb2 & 1 for hhb, nnb */

short numcovbonds = 2; /* number of covalent bonds that count as "nearly bonded" */
char    hbdn[160]; /* only used in main and procfile. HB directory */
FILE      *ofp=NULL; /*HB2 filename defined*/

short nnbflg = 0;


/* 6. Globals as name arrays */

extern short donors[MAXNAA][TOTNATM]; /* set in in_pdb */
extern short accepts[MAXNAA][TOTNATM];


/* 8. Vector/Matrix Routines */
/* 8.5 Globals as arrays */

 void enable_exchange(void)
 {
     accepts[His][7]=accepts[His][8]=donors[His][7]=donors[His][8]=1;
     accepts[Asn][7]=donors[Asn][6]=2;
     accepts[Gln][8]=donors[Gln][7]=2;
     return;
 }

 void enable_aromatic(void)
 {
     accepts[Phe][5]=accepts[Phe][6]=accepts[Phe][7]=accepts[Phe][8]=accepts[Phe][9]=accepts[Phe][10]= -1;
     accepts[Tyr][5]=accepts[Tyr][6]=accepts[Tyr][7]=accepts[Tyr][8]=accepts[Tyr][9]=accepts[Tyr][10]= -1;
     accepts[Trp][5]=accepts[Trp][6]=accepts[Trp][7]=accepts[Trp][9]=accepts[Trp][10]=accepts[Trp][11]=accepts[Trp][12]=accepts[Trp][13]= -1;
     /*printf("enable_aromatic has been called");*/

     return;
 }

 void disable_aromatic(void)
 {
     accepts[Phe][5]=accepts[Phe][6]=accepts[Phe][7]=accepts[Phe][8]=accepts[Phe][9]=accepts[Phe][10]=0;
     accepts[Tyr][5]=accepts[Tyr][6]=accepts[Tyr][7]=accepts[Tyr][8]=accepts[Tyr][9]=accepts[Tyr][10]=0;
     accepts[Trp][5]=accepts[Trp][6]=accepts[Trp][7]=accepts[Trp][9]=accepts[Trp][10]=accepts[Trp][11]=accepts[Trp][12]=accepts[Trp][13]=0;
     return;
 }

/* 8.7Minor useful functions */
/* 9. Hydrogen Position Calculation Subroutines */
/*10. Hydrogen Position Calculation Routine */


/*12. Connectivity Calculation Subroutines */
/*13. Connectivity Calculation Routine */
/*13.2 Domain calculation Routines */
/*13.5 Adding new residues to the selection */
/*14. find_hb */

/*****************************************************************************/
/*14. find_hb */

/* I'm making an exception to my own rules and declaring all my routines in
   advance.  This is because I didn't think moving all my routines around
   the file would be unhelpful. */

int is_donacc(int i, int j);
void find_hb (char * fname);
/*
int calc_hb(char * o_str,int i, int j);
void sprintf_hb(char * o_str, int i, int j, struct vect nearest, float d, float ha, float daaa_ang, float haaa_ang);
void fprintf_hb(char * o_str, int i, int j); True Globals so hbp_qnh.c can reach them */

int nhb = 0; /* Number of hydrogen bonds.  File Global so that fprint_hb can increment it.*/


void find_hb (char * fname)   
{ 
    int       i, ii, j, jj, k ;
    
    int       debug=0;
    char buf1[20],buf2[20]; /* buffers for atomid routines */
    
    char      calendar_s[128];
    
    /* void find_ss()*/
    /* {*/
    /*     int i,j;*/
    /*     float d;*/
    time_t    calendar;
    struct tm * calendar_tm_p;
    
    char output_string[256]="\0";
    
    /*     /* Check for disulphide bridges */
	
    if (!longoutflg)
    {
 	if (debug) printf("If !longoutflg\n");
 	
 	time ( &calendar );
 	if (debug) printf("  time\n");
 	
 	calendar_tm_p = localtime ( &calendar );
 	if (debug) printf("  calendar\n");
 	
 	strftime(calendar_s,128,"%b %d %X %Z %Y\0",calendar_tm_p);
 	if (debug) printf("  strftime\n");
 	
 	fprintf(ofp,"HBPLUS Hydrogen Bond Calculator v %s            %s\n",VERSION,calendar_s );
 	fprintf(ofp,"(c) I McDonald, D Naylor, D Jones and J Thornton 1993 All Rights Reserved.\n");
 	fprintf(ofp,"Citing HBPLUS in publications that use these results is condition of use.\n");
 	fprintf(ofp,"%4.4s <- Brookhaven Code \"%s\" <- PDB file\n",brcode,fname);
 	
 	fprintf(ofp,"<---DONOR---> <-ACCEPTOR-->    atom                        ^               \n");
 	fprintf(ofp,"c    i                          cat <-CA-CA->   ^        H-A-AA   ^      H- \n");
 	fprintf(ofp,"h    n   atom  resd res      DA  || num        DHA   H-A  angle D-A-AA Bond\n");
 	fprintf(ofp,"n    s   type  num  typ     dist DA aas  dist angle  dist       angle   num\n");
 	/*
	   <---DONOR---> <-ACCEPTOR-->    atom                        ^               
	   c    i                          cat <-CA-CA->   ^        H-A-AA   ^      H- 
	   h    n   atom  resd res      DA  ||   num      DHA   H-A  angle D-A-AA Bond
	   n    s   type  num  typ     dist DA aas  dist angle  dist       angle   num
	   E0010-LEU N   E0006-TYR O   3.07 MM   4  6.17 157.2  2.13 147.8 154.9    17
	   */
    }
    
    /* First Pass : Check for disulphide bridges */
    
    
    puts("Checking for hydrogen bonds . . . ");
    
    for (ii = 0; ii < (natoms-1); ii++)
    {
	/* ignore atom incapable of H-bonding immediately */
	if (atom[ii].aacode <TOTNAA &&
	    donors[atom[ii].aacode][atom[ii].atmtyp] == 0 &&
	    accepts[atom[ii].aacode][atom[ii].atmtyp] == 0)
	    if (!nnbflg)
 		continue;
	for (jj = ii + 1; jj < natoms; jj++)
	{
 	    debug=0;
 	    
 	    if (atom[jj].aacode <TOTNAA &&
		donors[atom[jj].aacode][atom[jj].atmtyp] == 0 &&
		accepts[atom[jj].aacode][atom[jj].atmtyp] == 0)
		if (!nnbflg)
 		    continue;
	    
 	    /* if the exchangeflg has been set to 2 by the -X option, only
 	       do the hydrogen bonds for the His, Asn, and Gln side-chains */
 	    if (exchangeflg==2)
 		if ( !( (atom[ii].atmtyp > 4 && (atom[ii].aacode == His 
						 || atom[ii].aacode == Asn || atom[ii].aacode == Gln) ) 
		       || ( atom[jj].atmtyp > 4 && (atom[jj].aacode == Asn 
 		    || atom[jj].aacode == Gln || atom[jj].aacode == His) ) 
		       ) )
 		    continue;
 	    
	    /* ignore very distant atom pairs immediately */
	    if (fabs(atom[ii].p.x - atom[jj].p.x) > HBDIST || 
		fabs(atom[ii].p.y - atom[jj].p.y) > HBDIST || 
		fabs(atom[ii].p.z - atom[jj].p.z) > HBDIST )
	    {
		if (debug == 2)
		    printf("                      very distant, ignoring this pair!\n");
		continue;
	    }
	    
	    
	    /* Check each IJ atom pair twice. I.e. check for I(don) / J(acc) 
	       and J(don) / I(acc) */
	    for (k = (nnbflg)?1:0; k < 2; k++) /* slightly opaque this one */
 		/* if nnbflg is set, then this loop only executes once,
 		   because it doesn't matter which is the donor and which is
 		   the acceptor */
	    {
		if (debug == 2)
		    printf("   k = %1d\n",k);
		/* Load/swap for correct donor/acceptor order */
		if (k)
		{
		    i = jj;
		    j = ii;
		}
		else
		{
		    i = ii;
		    j = jj;
		}
		
/*		if (!strncmp(atomid(i,buf1),"6XIA/-0002-HOH O  ",10) ||
		    !strncmp(atomid(i,buf2),"6XIA/-0038-SER OG ",10) )
				debug=2;
				else
				debug=0;*//*debug setting line*/
 		if (debug == 2)
 		    printf("Checking between %s and %s.\n",atomid(i,buf1),atomid(j,buf2));
 		
 		if( is_donacc(i,j) )
		{
		    if (debug==2)
			printf(". . . are Donor/Acceptor pair.\n");
		    
		    if( calc_hb(output_string,i,j))
		    {
			fprintf_hb(output_string,i,j);
			if (debug==2)
			    printf(". . . are H-Bonded.\n");
			
		    }
		    else
			if(debug==2)
			    printf(". . . but fail the calc_hb test\n");
		    
		}
		else
		    if (debug==2)
			printf(". . . NOT Donor/Acceptor pair.\n");
		
	    } /* end k loop */
 	}     /* end jj loop */
    }         /* end ii loop */
    if (nnbflg)
	printf("%d neighbour interactions found.\n", nhb);
    else
 	printf("%d hydrogen bonds found.\n", nhb);
}

/* }*/
/* */
/* hhb format
   3CHY 3CHY/-0007-LYS CB ? 3CHY/-0007-LYS O  ? 2.88SM   0-1.00 -1.0-1.00 -1.0
   123456789 123456789 123456789 123456789 123456789 123456789 123456789 12345
   */


int is_donacc(int i, int j)
{
    int debug=0;
    
    if (!nnbflg)
    {
	if (atom[i].aacode<TOTNAA && donors[atom[i].aacode][atom[i].atmtyp] == 0)
	    return(0); /*atom i is not recognized donor from a recognized residue*/
	if (debug==2)
	    printf("First atom is a donor\n");
	
	if (atom[j].aacode<TOTNAA && accepts[atom[j].aacode][atom[j].atmtyp] == 0)
	    return(0); /*atom j is not recognized acceptor from a recognized residue*/
	if (debug==2)
	    printf("Second atom is an acceptor\n");

	if ( atom[i].ssflg )
	    return(0); /* one or the other atom is in an S-S bridge */
	/* adjusted by IM to 'the donor is in a disulphide bridge' */
	if (debug==2)
	    printf("Not a ss bridge\n");
    }
    
    if (alreadybonded(i,j))
	return(0);
    	/* atoms i and j are bonded */
    if (debug==2)
	printf("Not already bonded\n");
    
    if (bondedwithin(i,j,numcovbonds))/* was  (nearlybonded(i,j)) */
	return(0); /* atoms i and j are 1-3 or 1-4 */
    /* adjusted by IM to be only 1-3 thrown out */
    if (debug==2)
	printf("pair not already bonded or nearly bonded\n");
    
    if (!nnbflg) 
    {
	if (atom[i].aacode == TOTNAA) /* atom i is from ligand/solvent */
	{
	    if ( atom[i].atmnam[1] == 'N'  ||
		(atom[i].atmnam[1] == 'O'  && iswater(atom[i].resnam)))
	    { 
		/* this ligand atom is potential donor */
	    }
	    else
		return(0); /* atom i not a donor, so next k */
	} 
	if (atom[j].aacode == TOTNAA && atom[j].atmnam[1] != 'O')
	    return(0); /* ignore ligand atoms as acceptors if they
			 are not oxygen */
    }
    return(1);
    
}

int calc_hb(char o_str[256],int i, int j)
{
    int       l, /*ri, rj, gap,*/ h, near_h;
    
    /*char      bndtyp[3], acctyp[4], dontyp[4], space=' ';*/
    float     /*ca_d,*/ d, ha2,/*IM*/tmp_d, hd, ha, haaa_ang, daaa_ang;
    struct vect nearest, tmp_point, aromatic_axis, aa[MAXCON] /* acceptor antecedants*/;
    
    
    char buf1[20],buf2[20]; /* buffers for atomid routines */
    
    /* void find_ss()*/
    /* {*/
    /*     int i,j;*/
    /*     float d;*/
    short aromatic_flg;
    

    *o_str=0;
    
/*    if(i==203 && j==210)
    {
	printf("Testing from %s",atomid(i,buf1));
	printf(" to %s at start of calc_hb.\n",atomid(j,buf1));
    	debug=0;
    }
    else
	debug=0;
  */  

    d = length_squared( to(atom[i].p, atom[j].p) );
    if (debug==2)
	printf("don-acc distance ^2 = %f\n",d);
    
    if (d <= SQR(HBDIST))
	/* these are nearby ! */; /*IM begins*/
    else return(0);
    
    if (debug) printf("Acceptable distance\n");
    
    /* time to fill the acceptor array */
    daaa_ang= -99.9; /* reset it . . . . */
    
    aromatic_flg =  ( atom[j].aacode < TOTNAA )?(accepts[atom[j].aacode][atom[j].atmtyp] < 0):0;
    
    if (!aromatic_flg)    
    {
	for(l=0;l<MAXCON;l++)
	{
	    if (icon[j][l]>-1) 
	    {
		float daaa_tmp,aaa;
		char buf1[20],buf2[20];
		
		aa[l]=atom[icon[j][l] ].p;
		
		daaa_tmp= dot_product(
 						  to( aa[l], atom[j].p ), to( atom[i].p, atom[j].p )       ); /* = DA * A-AA * cos D-A-AA */
 			    aaa=vector_length ( to(aa[l], atom[j].p) );
 			    
 			    daaa_tmp/= aaa * sqrt( d); /* = cos D-A-AA */
 			    /* d= DA^2 */
 			    if (daaa_tmp>daaa_ang) daaa_ang=daaa_tmp; /* select lowest */
 			    if (aaa<0.8 || aaa > 6.0) 
				printf("WARNING: %5.2fA between %s and %s ((D..)A-AA bond).\n", aaa,atomid(icon[j][l],buf1),atomid(j,buf2));
 			   
 			}
 		    }
 		    if (daaa_ang > cosMIN_DAAA && daaa_ang>-99.9) 
 			if (!nnbflg)
 			    return(0);
 		    
 		    /*NB -99 was not originally used as a default whilst checking
 		      because we were wrongly trying to find a lowest possible 
 		      value */
 		    /*	printf("lowest daaa_ang %6.1f\n",daaa_ang);*/
 		}
 		else /*aromatic_hbond*/
 		{
 		    float daaa_tmp;
 		    /*if( !strncmp("4DFR/A0030-TRP CE3",atomid(j,buf1),18))
		      debug=1;*/
 		    
 		    if (debug){
 			printf("NOTE: %s is aromatic.\n",atomid(j,buf1));
 			printf("      connects to ");
 			for(l=0;l< atom[j].ncon;l++)
 			    printf("%s ",atomid(icon[j][l],buf1));
 			printf("\n");
 			
 		    }
 		    
 		    if ( accepts[atom[j].aacode][atom[j].atmtyp] != -1 )
 			printf("BUG: %s not registered as aromatic acceptor in \"aromatic_hbond\".\n",atomid(j,buf1));
 		    if ( atom[j].ncon!=2 && atom[j].ncon!=3 )
 			printf("BUG: %s has connectivity %d but is registered as aromatic.\n",atomid(j,buf1),atom[j].ncon);
 		    if ( atom[j].ncon==3 )
 		    {
 			aromatic_axis = unit_vector(
						    vector_plus( perpendicular( atom[ icon[j][0] ].p, atom[j].p, atom[ icon[j][1] ].p),
								vector_plus( perpendicular( atom[ icon[j][0] ].p, atom[j].p, atom[ icon[j][2] ].p),
									    perpendicular( atom[ icon[j][1] ].p, atom[j].p, atom[ icon[j][2] ].p) ) ) );
 			
 		    }
 		    else
 			aromatic_axis = unit_vector(perpendicular( atom[ icon[j][0] ].p, atom[j].p, atom[ icon[j][1] ].p));
 		    
 		    aa[0]= vector_plus( aromatic_axis, atom[j].p );
 		    aa[1]= vector_plus( float_times_vect( -1.0, aromatic_axis), atom[j].p );
 		    if (debug) printf("Perpendicular gives angles %5.3f and %5.3f with 1,2 / %d.\n",angle( atom[ icon[j][0] ].p, atom[j].p, aa[0] ), angle( atom[ icon[j][1] ].p, atom[j].p, aa[0] ) , atom[j].ncon);
 		    
 		    for(l=0;l<2;l++)
 		    {
 			daaa_tmp= dot_product( to( aa[l], atom[j].p), to( atom[i].p, atom[j].p ) ); /*=DA*A-AA*cosD-A-AA */
 			daaa_tmp /= sqrt(d); /* = cos D-A-AX */
 			if ( daaa_tmp > daaa_ang ) daaa_ang=daaa_tmp; /* select high values of cos D-A-AX ie low, better angles*/
 			if (debug) printf("Cosine of D-A-AX = %5.3f\n",daaa_tmp);
 		        if (debug) {
 			    printf("AX co-ods");
 			    printf(TF,VXYZ( aa[l] ) );
 			    printf("\n");
 			}
 			
 		    }
 		    if (daaa_ang < cosMAX_DAAX)
 		    {
 			if (debug)
 			    printf("Rejected because DAAX too high.\n");
 			if (!nnbflg)
 			    return(0);
 		    }
 		}
 		
 		
#define nearer(vector) { tmp_d=length_squared( to( atom[j].p, vector ) ); if (ha2>tmp_d) { ha2=tmp_d; nearest=vector; near_h = h;} }
 		
 		/* find the nearest possible H position on donor i */
 		ha2=100;
		/*	if (i==362) debug=1;
			else debug=0;*/
 		
 		if (debug && atom[i].hetflg) printf("H on %s %d  ",atom[i].resnam, atom[i].aanum);
 		if (debug) printf("%d Hs positioned\n",atom[i].nh);
 		
 		for(h=atom[i].h_ptr;h<atom[i].nh+atom[i].h_ptr;h++)
 		{
 		    if (debug)
 			printf("Hydrogen %d - type %d\n",atom[i].h_ptr,h_atm[h].typ);
 		    if (debug) printf("atom[i].nh %d, ].h_ptr %d, h %d\n",atom[i].nh,atom[i].h_ptr,h);
 		    
 		    switch (h_atm[h].typ)
 		    {
 		    case 0:
 			if (debug) printf("Atom %d (%s%s %d) has H %d missing\n",i,atom[i].atmnam,atom[i].resnam, atom[i].aanum, h);
 			/* no H present, must be HETATM */
 			;/* used to be ha2=1, cannot decide why */
 			break;
 		    case alternatives:
 			nearer( h_atm[h].a );
 			if (debug) printf("Checked alternative position: ha2 = %f\n",ha2);			    			    
 		    case fixed:
 			nearer( h_atm[h].p );
 			if (debug) printf("Checked default position: ha2 = %f\n",ha2);
 			break;
 		    case circle:
 			tmp_point= intersect(to(atom[i].p,h_atm[h].a),
 					     h_atm[h].p, atom[i].p, atom[j].p );
 			tmp_point= onto_sphere(tmp_point, h_atm[h].a,
 					       vector_length( to(h_atm[h].a, h_atm[h].p) ));
 			nearer(tmp_point);
 			tmp_point= vector_plus( tmp_point,
 					       float_times_vect(2,to(tmp_point, h_atm[h].a)));
 			nearer(tmp_point);
 			break;
 		    default:
 			printf("\nBUG: Unforseen htyp for atom %s in find_hb.\n", atomid(i,buf1));
 			printf("Please mail mcdonald@uk.ac.ucl.biochemistry.bsm\n");
			
 		    }/* end of switch*/
 		}/* end of h loop */
 		
 		if (ha2<100) /* ie if any Hydrogens have been found and
 				their positions taken.  Or, to put it 
 				another way, if it anything other than a hetatm 
 				which is flagged in .nh as having hydrogens but
 				has no positons in h_atm[].typ. */
 		{
 		    
 		    hd=vector_length( to(nearest, atom[i].p));
 		    ha=vector_length( to(nearest, atom[j].p));
 		    if (hd<0.9 || hd > 2.0) 
 		    {
 			printf("WARNING: %5.2fA between %s and %d (H).\n", hd, atomid(i,buf1),near_h-atom[i].h_ptr+1);		
 			if (debug) printf("atom[i].nh %d, ].h_ptr %d, h %d\n",atom[i].nh,atom[i].h_ptr,h);
 		    }
 		    
		    if (debug) printf("hd=%f,  ha=%f\n",hd,ha);
 		    
 		    if (!nnbflg)
 		    {
 			if (ha>MAX_HA)
 			    return(0);
 			if (debug)
 			    printf("HA length OK\n");
 			
 			if (dot_product( to(nearest, atom[i].p), 
 					to(nearest, atom[j].p) ) >= hd*ha*cosMIN_DHA)
 			    return(0);
 		    }
 		    
 		    haaa_ang=100.0; /* reset it . . . . */
 		    
 		    if (!aromatic_flg)
 		    {	
 			for(l=0;l<MAXCON;l++)
 			{
 			    if (icon[j][l]>-1) 
 			    {
 				float haaa_tmp,aaa;
 				aa[l]=atom[icon[j][l] ].p;
 				haaa_tmp= dot_product(
 						      to(atom[j].p, aa[l] ), to( atom[j].p, nearest )       );
 				aaa= vector_length( to(aa[l],atom[j].p) );
 				
 				haaa_tmp/= aaa *ha;
 				if (haaa_tmp<haaa_ang) haaa_ang=haaa_tmp; /* select lowest */
 				if (aaa<0.8 || aaa > 6.0) 
 				    printf("WARNING: %5.2fA between %s and %s ((H..)A-AA).\n", aaa,atomid(j,buf1),atomid(icon[j][l],buf2));		    
 			    }
 			}
 			if ((haaa_ang > cosMIN_HAAA && haaa_ang<100.0) && !nnbflg) 
 			    return(0);
 			else
 			    if (haaa_ang==100.0)
 				haaa_ang= -99.9;
 			/*NB -99 could not have been used as a default whilst checking
 			  because we are trying to find a lowest possible value */
 			
 			if (debug)
 			    printf("Angle OK\n");
 		    }
 		    else
 		    {
 			float haaa_tmp;
 			haaa_ang= -99; /*dealing with hAax, prefer small, ie big cosines*/
 			for(l=0;l<2;l++)
 			{
 			    
 			    haaa_tmp=dot_product( to(atom[j].p,aa[l]),to(atom[j ].p,nearest));
 			    haaa_tmp/= ha;
 			    if (debug) printf("Cosine of H-A-AX = %5.3f\n",haaa_tmp);
 			    if (haaa_tmp > haaa_ang) {
 				if (debug) printf("%5.3f > %5.3f.\n",haaa_tmp,haaa_ang);
 				
 				haaa_ang=haaa_tmp; /* select the low ie more meaningful */
 				if (debug) printf("Lowest yet.\n");
 			    }
 			    
 			}
 			if (haaa_ang < cosMAX_HAAX /* && haaa_ang<100.0 */ )
 			{
 			    if (debug) {
 				printf("Rejected because of wide HAAX angle.\n");
 				printf("Cos haaa_ang = %5.3f, cosMAX_HAAX = %5.3f.\n", haaa_ang, cosMAX_HAAX);
 			    }
 			    if (!nnbflg)
 				return(0);
 			}
 		    }
 		    
 		    
 		}/* end of 'if any Hs actually found' clause */
 		else
 		    if (d>SQR(MAX_HA+1.0))
 			if (!nnbflg)
 			    return(0); /* if the hydrogens are unfixed, knock
					 anything more than 1.0 Ang - the 
					 normal DH distance - back. */

    sprintf_hb(o_str,i,j,nearest,d,ha,daaa_ang,haaa_ang);
    return(1);
}

void sprintf_hb(char * o_str, int i, int j, struct vect nearest, float d, float ha, float daaa_ang, float haaa_ang)
{
    char      bndtyp[3], acctyp[4], dontyp[4], space=' ';
    int       ca_d, ri, rj, gap;
    char buf1[20],buf2[20]; /* buffers for atomid routines */
    float     tmp;
        
    /*we have a hydrogen bond!*/
    
    bndtyp[0] = bndtyp[1] = '?';
    acctyp[3] = dontyp[3] = bndtyp[2] = '\0';
    
    if (atom[i].aacode == TOTNAA || atom[i].hetflg )
	bndtyp[0] = 'H';
    else
    {
	if (atom[i].aacode <STDAA)
	{
	    bndtyp[0] = (atom[i].atmtyp <=3) ? 'M' : 'S';
	}
	else
	{
	    bndtyp[0] = (atom[i].atmtyp <=3) ? 'm' : 's';
	}
    }
    if (atom[j].aacode == TOTNAA || atom[j].hetflg )
	bndtyp[1] = 'H';
    else
    {
	if (atom[j].aacode <STDAA)
	{
	    bndtyp[1] = (atom[j].atmtyp <=3) ? 'M' : 'S';
	}
	else
	{
	    bndtyp[1] = (atom[j].atmtyp <=3) ? 'm' : 's';
	}
    }
 		
    ri = atom[i].caindex;
    rj = atom[j].caindex;
    
    gap=find_gap( ri, rj);
		
#if 0 		
    if (atom[i].aacode==TOTNAA || (ri == -99) || !residue[ri].ca ||                       atom[j].aacode==TOTNAA || (rj == -99) || !residue[rj].ca  )
	gap = -1;
    else
    {
	gap = cagap(ri,rj);
    }
#endif 		
    if (atom[i].aacode==TOTNAA || (ri == -99) || !residue[ri].ca ||  
	atom[j].aacode==TOTNAA || (rj == -99) || !residue[rj].ca  )
	ca_d = -1.0;
    else
    {
	ca_d = length_squared( to(*(residue[ri].ca),*(residue[rj].ca)) );
	
	if (ca_d > SQR(CAWARN) && debug)
	{
	    printf("WARNING: unusual CA-CA distance of %5.2f for",sqrt(ca_d));
	    printf("%s number ??\n", (nnbflg)?" interaction":" H-bond"/*,nhb*/);
	    printf("         between %s and %s\n", atomid(i,buf1),atomid(j,buf2));
	    
	}
    }
#if 0
    printf("%4s",   brcode);
    printf("%c",    (atom[i].chnid == '-') ? space : atom[i].chnid);
    printf("%4s/",  brcode);
    printf("%c",    atom[i].chnid);
    printf("%04d",  atom[i].aanum);
    if(OMLINSERT && atom[i].hetflg)
	printf("h");
    else
	printf("%c",    atom[i].inscode);
    /* printf("%3s",   dontyp);               */
    printf("%3s",   atom[i].resnam);
    printf("%4s",   atom[i].atmnam);
    printf("%c",    atom[i].strucsum);
    printf("%c",    (atom[j].chnid == '-') ? space : atom[j].chnid);
		printf("%4s/",  brcode);
		printf("%c",    atom[j].chnid);
		printf("%04d",  atom[j].aanum);
    if(OMLINSERT && atom[j].hetflg)
	printf("h");
    else
	printf("%c",    atom[j].inscode);
    /* printf("%3s",   acctyp);               */
		printf("%3s",   atom[j].resnam);
		printf("%4s",   atom[j].atmnam);
		printf("%c",    atom[j].strucsum);
		printf("%4.1f", sqrt(d));
		printf(" %2s",   bndtyp);
    printf("%4d",   gap);
    printf("%4.1f", (ca_d <= 0.0) ? -1.0 : sqrt(ca_d));
    printf("%6.1d", 180/3.1415927*angle(atom[i].p,nearest,atom[j].p) );
    printf("%4.1d", d);
    printf("%6d\n", nhb);
		
#endif
    *o_str = 0 ; 

/*This isn't neccessary if the longoutput format is used, but if the short
  format is used, then the first output must (i) include a 'strlen' in case
  it is part of a 'long format' output, and (ii) have some way of knowing
  where to begin in case of a 'short format' output. */
    
    if (longoutflg)  
 		{
 		    /* start to output the line, *.hhb format */
 		    sprintf(o_str              , "%4s",   brcode);
 		    sprintf(o_str+strlen(o_str), "%c",    (atom[i].chnid == '-') ? space : 
 			    atom[i].chnid);
 		    sprintf(o_str+strlen(o_str), "%4s/",  brcode);
 		}
 		sprintf(o_str+strlen(o_str), "%c",    atom[i].chnid);
 		sprintf(o_str+strlen(o_str), "%04d",  atom[i].aanum);
if(OMLINSERT && atom[i].hetflg) /*v3.13*/
    sprintf(o_str+strlen(o_str), "h");
else
    sprintf(o_str+strlen(o_str), "%c",    atom[i].inscode);
		/* sprintf(o_str+strlen(o_str), "%3s",   dontyp);               */
/* Amendment. RAL 14 Jun 2012 --> */
// 		if (atom[i].aacode==Cys && atom[i].atmtyp == 5 & cssflg)
 		if (atom[i].aacode==Cys && atom[i].atmtyp == 5 && cssflg)
/* <-- Amendment. RAL 14 Jun 2012 */
 		{
 		    if (atom[i].ssflg) 
 			sprintf(o_str+strlen(o_str),"CSS");
 		    else 
 			sprintf(o_str+strlen(o_str),"CYH");    }
 		else
 		    sprintf(o_str+strlen(o_str), "%3s",   atom[i].resnam);
 		sprintf(o_str+strlen(o_str), "%4s",   atom[i].atmnam);
 		if (longoutflg)  
 		{
 		    
 		    sprintf(o_str+strlen(o_str), "%c",    atom[i].strucsum);
 		    sprintf(o_str+strlen(o_str), "%c",    (atom[j].chnid == '-') ? space : 
 			    atom[j].chnid);
 		    sprintf(o_str+strlen(o_str), "%4s/",  brcode);
 		}
 		else
 		    sprintf(o_str+strlen(o_str), " ");
 		
 		sprintf(o_str+strlen(o_str), "%c",    atom[j].chnid);
                sprintf(o_str+strlen(o_str), "%04d",  atom[j].aanum);
    if(OMLINSERT && atom[j].hetflg) /*v3.13 start*/
	sprintf(o_str+strlen(o_str),"h");
    else /*v3.13 end*/
 		sprintf(o_str+strlen(o_str), "%c",    atom[j].inscode);
 		/* sprintf(o_str+strlen(o_str), "%3s",   acctyp);               */
 		if (atom[j].aacode==Cys && atom[j].atmtyp == 5 && cssflg)
 		{
 		    if (atom[j].ssflg) 
 			sprintf(o_str+strlen(o_str),"CSS");
 		    else 
 			sprintf(o_str+strlen(o_str),"CYH");    }
 		else
 		    sprintf(o_str+strlen(o_str), "%3s",   atom[j].resnam);
 		sprintf(o_str+strlen(o_str), "%4s",   atom[j].atmnam);
 		if (longoutflg) 
 		    sprintf(o_str+strlen(o_str), "%c",    atom[j].strucsum);
 		sprintf(o_str+strlen(o_str), "%5.2f", sqrt(d));
 		if (!longoutflg) 
 		    sprintf(o_str+strlen(o_str), " ");
 		sprintf(o_str+strlen(o_str), "%2s",  bndtyp);
 		sprintf(o_str+strlen(o_str), "%4d",   gap);
 		if (!longoutflg) 
 		    sprintf(o_str+strlen(o_str), " ");
 		sprintf(o_str+strlen(o_str), "%5.2f", (ca_d <= 0.0) ? -1.0 : sqrt(ca_d));
 		if (!longoutflg) 
 		    sprintf(o_str+strlen(o_str), " ");
 		if (atom[i].nh && h_atm[atom[i].h_ptr].typ)
 		{
 		    /* if it has hydrogens and if they have been placed */
		    tmp = angle(atom[i].p,nearest,atom[j].p);
		    
 		    sprintf(o_str+strlen(o_str), "%5.1f", (180.0/3.1415927)* tmp );
		    if (!longoutflg) 
			sprintf(o_str+strlen(o_str), " ");
 		    sprintf(o_str+strlen(o_str), "%5.2f", ha);
		    if (!longoutflg)
			sprintf(o_str+strlen(o_str), " ");
 		    if (haaa_ang>-50.0)
 			sprintf(o_str+strlen(o_str),"%5.1f", 180/3.1415927* (float) acos(haaa_ang) );
 		    else
 			sprintf(o_str+strlen(o_str),"%5.1f", -1.0);
		    
 		}
 		else
 		    if (longoutflg) 
 			sprintf(o_str+strlen(o_str), "%5.1f%5.2f%5.1f", -1.0, -1.0, -1.0);
 		    else
 			sprintf(o_str+strlen(o_str), "%5.1f %5.2f %5.1f", -1.0, -1.0, -1.0);
 		if (daaa_ang>-50.0)
 		    sprintf(o_str+strlen(o_str),"%6.1f", 180/3.1415927*acos(daaa_ang) );
 		else
 		    sprintf(o_str+strlen(o_str),"%6.1f", -1.0);
		
 		if (debug)
 		{
 		    sprintf(o_str+strlen(o_str), "D " );
 		    sprintf(o_str+strlen(o_str), TF,VXYZ(atom[i].p) );
 		    sprintf(o_str+strlen(o_str), " H " );
 		    sprintf(o_str+strlen(o_str), TF,VXYZ(nearest) );
 		    sprintf(o_str+strlen(o_str), " A " );
 		    sprintf(o_str+strlen(o_str), TF, VXYZ(atom[j].p) );
 		}
 		
/* 		sprintf(o_str+strlen(o_str), "%6d\n", nhb);
		THIS IS DONE BY THE CALLING ROUTINE */

 		/* end of outputting the line, *.hhb format */

    return;

}

void fprintf_hb(char * o_str, int i, int j) /*old code moved v3.08.  Name explanatory*/
{
    static int nda = 0;
    fprintf(ofp,"%s%6d\n",o_str,nhb+1);
    /* *output_string=0 I can return this if I need to, but I don't see why it's needed v3.08 */;
    if(!nhb) /*If this is the first HB of the file ,- 3.10 */
	nda = 0; /*Start at zero again <- 3.10*/

      /*Bugfix.  Without this, any lengthy multiple pdbfile run sent the
        value of nda through the ceiling and ended up over-writing the next
        file.  Which I hadn't foolproofed against*/
 
    nhb++; /*File global*/
    if( !atom[i].don_p )
	atom[i].don_p = &hb_prtnr[nda++];
    if( !atom[j].acc_p )
	atom[j].acc_p = &hb_prtnr[nda++];
    if(nda >= MAXNHB_ATOMS )
	fail("FATAL ERROR: Too Many Hydrogen Bonding Atoms");
        
    atom[i].don_p -> n[ atom[i].ndonhb ++] = j;
    if(atom[i].ndonhb > MAXNHB_PER_ATM)
	fail("FATAL ERROR: Too Many Hydrogen Bonds per Single Atom");
    atom[j].acc_p -> n[ atom[j].nacchb ++] = i;
    if(atom[j].nacchb > MAXNHB_PER_ATM)
	fail("FATAL ERROR: Too Many Hydrogen Bonds per Single Atom");
    
}

/*15. procfile */
/*16. main */

