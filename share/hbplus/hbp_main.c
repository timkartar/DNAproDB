/* HBPLUS - Hydrogen Bond Calculator control file v 3.2 */
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

   Version 2.23    RAL 1 May 1997
                   Bug-fix to PDB-reading routine which was ignoring all
                   HETATM records beyond atom number 9999 because of use
                   of scanf to read in the key-word.
                   Amendment for NMR structure so that program stops
                   reading ATOM and HETATM records after it has encountered
                   an ENDMDL record. Similarly to ignore CONECT records
                   relating to any but the first model.

   Version 2.24    RAL 14 May 1997
                   Amendment to routines reading in hbplus.rc which fail on
                   atom names containing a hash, "#" (eg in 2hwb, 2hwc, 2hwd
                   and 2hwe).
                   RAL 25 May 1997
                   Amendment to routines reading in hbplus.rc which fail on
                   atom names containing a double quote (eg 1gac). Program
                   hbadd altered to replace double-quotes by @'s which are
                   here converted back.

   Version 2.25    RAL 24 Jun 1997
                   Increased size of MAXBRKS from 2000 to 5000 as program
                   considers there to be "chain-breaks" between waters and
                   some structures have over 2000 waters (eg 1wap and
                   1vps).

   Version 2.26    RAL 17 Aug 1999
                   Bug-fix to prevent crashing under linux when closing
                   SST file with a NULL pointer.

   Version 2.26a   RAL 17 Nov 2007
                   Amendment for new format for DNA ATOM records.

   Version 2.27    RAL 24 Apr 2008
                   Bug-fix to check for NULL environment variable HOME.
                   Failure to check this was causing program to crash when
                   run from a CGI-script by w3nobody

   Version 2.28    RAL 9 Mar 2009
                   Amendment for call by LigPlus s.t. output directory is
                   supplied on the command line and used for naming the
                   hbdebug.dat file
                   RAL 13 Aug 2009
                   Change to parsefn to cope with backslashes in filenames

   Version 2.29    RAL 3 Mar 2011
                   Amendment for new format for RNA ATOM records.

   Version 2.30    RAL 3 Dec 2013
                   Amendment for new format for hydrogen atom names.

/* Contents */
/* 1. Copyright - all*/
/* 2. Version Log - all*/
/* 3. Tokens as abbreviations */

#include <time.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "hbplus.h"

/* Amendment. RAL 3 Jul 2012 --> */
/* Prototypes */
short inpdb_file(char * fname, char * inpdbfn);
int string_truncate(char *string,int max_length);
/* <-- RAL 3 Jul 2012 */


/* 4. Tokens as flags and customisables */
/* 5. Globals as flags and customisables */
/* 6. Globals as name arrays */
/* 8. Vector/Matrix Routines */
/* 8.5 Globals as arrays */
/* 8.7Minor useful functions */

 /* token dividing routine, similar to strtok, but allows for ", ' & ` quotes*/
 char * strqtok(char * buffer,char * brkchr, char * quote)
 {
     static char* pointer; /* when does the next one start ? */

     if (!buffer)
/* Amendment. RAL 14 Jun 2012 --> */
       {
/* <-- Amendment. RAL 14 Jun 2012 */
         if (pointer)
           buffer=pointer;
         else
           printf("BUG: strqtok called with no buffer and pointer unset\n");
/* Amendment. RAL 14 Jun 2012 --> */
       }
/* <-- Amendment. RAL 14 Jun 2012 */

     buffer+=strspn(buffer,brkchr);/*buffer=argstart*/

     if (!*buffer)
        return(NULL);

     if (strchr("\"\'",(int)*buffer))
     {
        pointer=strchr(buffer+1,(int)*buffer); /*pointer=matching quote*/
        buffer++; /*obviously, the start of the argument*/

        *pointer='\0';
        if (*(pointer+1))/*if closing quote isn't last character . . .*/
            pointer++;
     }
     else
     {
        pointer=buffer+strcspn(buffer,brkchr);/*pointer=next space start*/
        if (*pointer){
            *pointer='\0';
            pointer++;
        }
     }
     return(buffer); /*buffer holds the address of the next argument*/
 }

 /* parsing routine, dividing "buffer" into the equivalent of argv[]. */
 int parse_line(char buffer[256], char * values[128])
 {
     int i=0;
     if ( (values[i])= (strqtok(buffer,"\t\n ","")) )
     {
        i++;
        while ( (values[i])= (strqtok(NULL,"\t\n ","")) )
            i++;
     }
     return (i);
 }


/* 9. Hydrogen Position Calculation Subroutines */
/*10. Hydrogen Position Calculation Routine */
/*11. PDBOUT */
/*12. Connectivity Calculation Subroutines */
/*13. Connectivity Calculation Routine */
/*13.2 Domain calculation Routines */
/*13.5 Adding new residues to the selection */
/*14. find_hb */
/*15. inpdb_file */
void do_file(char * fname, char * inpdbfn)
{
    char      outfn[128]="\0", pdbhfn[128]="\0", basename[128]="\0";
    
    printf("\n");
    
/*    printf("DEBUG: About to parse fname\n");*/
    
    parsefn(fname,basename);
/*    printf("DEBUG: Parsed fname\n");*/
    
    if ( inpdb_file(fname,inpdbfn) )
    {
        find_h();

        if( asaflg)
        {
            /*parsefn(fname, pdbhfn); -> 3.01*/
            strcpy(pdbhfn, basename); /* <- 3.01 */
            
            strcat(pdbhfn, "asa");
            inasa_file(pdbhfn);
        }
        
        /*parsefn(fname, pdbhfn); -> 3.01*/
        strcpy(pdbhfn, basename); /* <- 3.01 */
        
        strcat(pdbhfn, "h\0");

        if (pdboutflg) {
            /*        printf("Calling pdbout with %s and %s\n",pdbhfn,fname);*/
            
            pdbout(pdbhfn,fname,0,(short int *)NULL);
        }

        strcpy(outfn, hbdn);
        strcat(outfn, basename); /* <- 3.01 */
        
        /*newparsefn(fname, outfn+strlen(hbdn), 128); /* -> 3.01 */
        /*I think this puts the 4 letter code into outfn*/
        if (longoutflg)
        if (nnbflg)
            strcat(outfn,"nnb");
        else
            strcat(outfn,"hhb");
        else
            if (nnbflg)
                strcat(outfn,"nb2");
            else
                strcat(outfn,"hb2");
        if (debug) printf("About to open output file\n");

        ofp = fopen(outfn, "w");
        /*ofp output file pointer - external*/
        if (!ofp)
        {
            printf("^GFailed to open output file %s.\n",outfn);
            return;
        }
        else
        {
            printf("Opened output file \"%s\".\n",outfn);
            find_hb(fname);
        }
#ifdef DOM        
        if (domainflg)
            find_domains(fname);
#endif        
#ifdef BSM
        if (bsmoptflg==1)
        {
            printf("DEBUG: About to do_msw\n");
            
            do_msw(basename);
        }
        
#endif
        if (exchangeflg==2)
            chkqnh();
        
        if (ofp) /*this used to be (!hbdn[0]) for some strange reason*/
            fclose(ofp);
        /*printf_numhb();*/ /* debugging line */

    }
    
}

void procinfilename(char * infilename,char * inpdbfn)
{
    char infilename2[128];
    FILE * bfp=NULL;
     
    if (infilename[0] == '@')
    {
         bfp = fopen(infilename + 1, "r");
         if (!bfp)
         {
             printf("Failed to open file of file names as input\n");
             return;
         }
         while (!feof(bfp))
         {
             if (!fgets(infilename, 128, bfp))
             {
                 fclose(bfp);
                 return;
             }
             sscanf(infilename," %s",infilename2);
             debug=0;
             if (debug) printf("from %s to %s\n",infilename,infilename2);
             debug=0;
             
             /*infilename[strlen(infilename) - 1] = '\0';*/
             do_file(infilename2, inpdbfn);
         }
         fclose(bfp);
     }
     else
     {
         do_file(infilename, inpdbfn);
     }
     return;
     
 }

/*16. main */

/****************************************************************************/
/*16. main */
/* Amendment. RAL 3 Jul 2012 --> */
//    char            tmpstr[128]="\0",inpdbfn[128]="\0";
    char            tmpstr[FILENAME_LEN], inpdbfn[FILENAME_LEN];
/* <-- RAL 3 Jul 2012 */

/* tmpstr is the name of the 'new' file.  inpdbfile is the name of the 'pdb' */
/* file.  They are global variables to allow them to be manipulated by more  */
/* than one routine - parse_options and main(), principally. */

/****************************************************************************/
/* but first, the routine to parse the options                              */
/* Syntax hbplus.exe [-a min] [-A dha haaa daaa] [-hH ha] [-dD da] [-Ss ss]
                     [-oO] (pdb file with Hs) [-cC] (cys->css or cyh)
                     [-xX] (exchange) -[rR] (aromatic) 
                     [-bB aromatic_angle ]
                     [-E 'donname' num] [-e 'accname' num]
                     [-nN] (neighbours) -[vV] (covalent bonds)
                     [-uUTM . . . defining your own atoms]
                     [-f optionfile ]
                     [-P] (print out atom, etc list)
                     [-Z] (calculate domains)
                     [-Y] (pepbndval vdvval)
                     [-K] (Kabsch/Sander NHs, default off)
                     [-q] Input asa file
                     [-Q] Input asa file & output HB data
                     [-W ?] BSM option - produce various outputs for BSM projects.
*/
/*This must be declared, as we have two functions that can call each other*/
void parse_options_file(char * filename);

/* for my reference, those options that have and have not been used . . */
/* ====== == = === === == =   =====  == = ==== =========
   abcdefghijklmnopqrstuvwxyz ABCDEFGHIJKLMNOPQRSTUVWXYZ */
/* Amendment. RAL 9 Mar 2009 --> */
/* void parse_options(int optc, char ** optv, int start) */
void parse_options(int optc, char ** optv, int start, char *wkdir)
{
/* <-- RAL 9 Mar 2009 */
  int i,j;
  float minang, cosMIN_ANG, sinMIN_ANG;

  /* Amendment. RAL 9 Mar 2009 --> */
  // Initialise the working directory
  wkdir[0] = '\0';
  /* <-- RAL 9 Mar 2009 */

  // Loop over the command line arguments
  for (i=start;i<optc;i++)
    {
      int jflg=1;
        
      if (optv[i][0]=='-')
        {
/* Amendment. RAL 9 Mar 2009 --> */
          // Check the -wkdir option
          if (start == 1 && !strcmp(optv[i],"-wkdir"))
            {
              // Retrieve the working directory from the next argument
              if (i+1>=optc)
                fail("SYNTAX ERROR: Not Enough Arguments for -wkdir option\n");
              strcpy(wkdir,optv[i+1]);
              i++;

              // Unset flag
              jflg = 0;
            }
/* <-- RAL 9 Mar 2009 */

          // Loop over the characters in this token
          for(j=1;optv[i][j] && jflg;j++)
            {
              switch(optv[i][j])
                { 
                case 'a':
                  if (i+1>=optc) fail("SYNTAX ERROR: Not Enough Arguments for -a option\n");
                    
                  sscanf(optv[i+1],"%f",&minang);
                  sinMIN_ANG = (float)sin((minang*3.1415927)/180);
                  cosMIN_ANG = (float)cos((minang*3.1415927)/180);
                  /*j=strlen(optv[++i]);*/
                  /* have to do something to break out of both case & loop */
                  i++;jflg=0;
                    
                  mindha=minang;
                  minhaaa=minang;
                  mindaaa=minang;
                  sinMIN_DHA=sinMIN_ANG;
                  sinMIN_HAAA=sinMIN_ANG;
                  sinMIN_DAAA=sinMIN_ANG;
                  cosMIN_DHA=cosMIN_ANG;
                  cosMIN_HAAA=cosMIN_ANG;
                  cosMIN_DAAA=cosMIN_ANG;
                  break;
                case 'A':
                  if (i+3>=optc) fail("SYNTAX ERROR: Not Enough Arguments for -A option\n");
                  sscanf(optv[i+1],"%lf",&mindha);
                  sinMIN_DHA = (float)sin((mindha*3.1415927)/180);
                  cosMIN_DHA = (float)cos((mindha*3.1415927)/180);
                  sscanf(optv[i+2],"%lf",&minhaaa);
                  sinMIN_HAAA = (float)sin((minhaaa*3.1415927)/180);
                  cosMIN_HAAA = (float)cos((minhaaa*3.1415927)/180);
                  sscanf(optv[i+3],"%lf",&mindaaa);
                  sinMIN_DAAA = (float)sin((mindaaa*3.1415927)/180);
                  cosMIN_DAAA = (float)cos((mindaaa*3.1415927)/180);
                  i+=3;/*j=strlen(optv[i]);*/
                  jflg=0;
                    
                  break;
                case 'b':
                  if (i+1>=optc) fail("SYNTAX ERROR: Not Enough Arguments for -b option\n");
                  sscanf(optv[i+1],"%lf",&maxdaax);
                  sinMAX_DAAX = (float) sin((maxdaax*3.1415927)/180);
                  cosMAX_DAAX = (float) cos((maxdaax*3.1415927)/180);
                  i++;
                  jflg=0;
                  maxhaax=maxdaax;
                  sinMAX_HAAX = sinMAX_DAAX;
                  cosMAX_HAAX = cosMAX_DAAX;
                  break;

                case 'B':
                  if (i+2>=optc) fail("SYNTAX ERROR: Not Enough Arguments for -B option\n");
                  sscanf(optv[i+1],"%lf",&maxhaax);
                  sinMAX_HAAX = (float) sin((maxhaax*3.1415927)/180);
                  cosMAX_HAAX = (float) cos((maxhaax*3.1415927)/180);
                  sscanf(optv[i+2],"%lf",&maxdaax);
                  sinMAX_DAAX = (float) sin((maxdaax*3.1415927)/180);
                  cosMAX_DAAX = (float) cos((maxdaax*3.1415927)/180);
                  i+=2;
                  jflg=0;
                  break;
                    
                case 'h':
                case 'H':
                  if (i+1>=optc) fail("SYNTAX ERROR: Not Enough Arguments for -h/H option\n");
                  sscanf(optv[i+1],"%f",&MAX_HA);
                  CAWARN = MAX_HA + 8.5;
                    
                  /*j=strlen(optv[++i]);*/
                  jflg=0;
                  i++;
                  break;
                case 'd':
                case 'D':
                  if (i+1>=optc) fail("SYNTAX ERROR: Not Enough Arguments for -d/D option\n");
                  sscanf(optv[i+1],"%f",&HBDIST);
                  /*j=strlen(optv[++i]);*/
                  jflg=0;
                  i++;
                  break;
                case 'e': /* set an acceptor */
                case 'E': /* set a donor */
                  if (i+3>=optc) fail("SYNTAX ERROR: Not Enough Arguments for -e/E option\n");
                  /*int add_donacc ( char * residue, char * atoms, short number, short donflg )*/
                  jflg=0;
                  {
                    int hbnum,oldi;
                    oldi=i;
                    i+=3;
                                        
                    if (!sscanf(optv[oldi+3],"%d",&hbnum))
                      {
                        printf("ERROR: ARGUMENT # %d \"%s\" NOT A VALID INTEGER.\n",oldi+2,optv[oldi+2]);
                        break;
                      }
                    
                    add_donacc(optv[oldi+1], optv[oldi+2],hbnum,optv[oldi][j]=='E');
                    break;
                  }
                    
#if 0
                  {
                    char atmstr[8];
                    int hbnum,res,atm,oldi;
                    
                    jflg=0;
                    oldi=i;
                    i+=2;
                                        
                    if (strcpy(atmstr,optv[oldi+1])==NULL )
                      {
                        printf("ERROR: INVALID ARGUMENT # %d\n",oldi+1);
                        break;
                      }
                    if (!sscanf(optv[oldi+2],"%d",&hbnum))
                      {
                        printf("ERROR: ARGUMENT # %d \"%s\" NOT A VALID INTEGER.\n",oldi+2,optv[oldi+2]);
                        break;
                      }
                    
                    if (!scanresatmnam(atmstr,&res,&atm)) /* 0=success*/
                      {
                        printf("Residue %.3s atom %.4s is now extra H-Bond ",rnames[res],necatm[res]+4*atm);
                        
                        if (optv[oldi][j]=='E')
                          {
                            donors[res][atm]=hbnum;
                            printf("donor.\n");
                          }
                        else
                          {
                            accepts[res][atm]=hbnum;
                            printf("acceptor.\n");
                          }
                      }
                    else
                      printf("ERROR: ATMID \"%s\" UNRECOGNISED\n",atmstr);
                    
                    
                    break;
                  }
#endif
                case 'f':
                  if (i+1>=optc) fail("SYNTAX ERROR: Not enough Arguments for -f option\n");
                  jflg=0;
                  i++;
                                        
                  parse_options_file(optv[i]);
                  break;
                    
                case 'k':
                  kshflg=0;
                  break;
                case 'K':
                  kshflg=1;
                  break;
                    
                case 'M': /* add atoM */
                  if (i+2>=optc) fail("SYNTAX ERROR: Not Enough Arguments for -m/M option\n");
                  add_atoms(optv[i+1],optv[i+2]);
                  i+=2;
                  jflg=0;
                  break;
                case 'r': /*unset aromatic types*/
                  disable_aromatic();
                  break;
                case 'R': /*set aromatic types*/
                  enable_aromatic();
                  break;
                case 's':
                case 'S':
                  sscanf(optv[i+1],"%f",&SSDIST);
                  jflg=0;
                  i++;
                  break;
                case 'T': /*connecT atoms*/;
                  if (i+2>=optc) fail("SYNTAX ERROR: Not Enough Arguments for -T option\n");
                  add_bonds(optv[i+1],optv[i+2]);
                  printf("Adding bonds to residue %s.\n",optv[i+1]);
                    
                  i+=2;
                  jflg=0;
                  break;
                case 'U': /*add residUe*/;
                  if (i+2>=optc) fail("SYNTAX ERROR: Not Enough Arguments for -U option\n");
                  add_residue_type(optv[i+1],optv[i+2]);
                  i+=2;
                  jflg=0;
                  break;
                case 'u': /*add residUe without simile*/;
                  if (i+1>=optc) fail("SYNTAX ERROR: Not Enough Arguments for -U option\n");
                  add_residue_type(optv[i+1],"");
                  i++;
                  jflg=0;
                  break;
                case 'o':
                case 'O':
                  pdboutflg=1;
                  break;
                case 'c':
                case 'C':    
                  cssflg=1;
                  break;
                case 'x':
                  exchangeflg=1;
                  enable_exchange();
                  break;
                case 'X':
                  exchangeflg=2;
                  enable_exchange();
                  break;
                case 'I':
                  inputsstflg=1;
                  break;
                case 'i':
                  inputsstflg=0;
                  break;
                case 'L':
                  longoutflg=1;
                  break;
                case 'l':
                  longoutflg=0;
                  break;
                case 'n':
                  nnbflg=0;
                  break;
                case 'N':
                  nnbflg=1;
                  break;
                case 'P':
                  {
                    int i,j;
                    printf("Listing of residue types\n");
                    for(i=0;*rnames[i];i++)
                      {
                        printf("\n%s\n",rnames[i]);
                        printf("%s\n",necatm[i]);
                        for(j=0;*(necatm[i]+4*j);j++)
                          if (donors[i][j])
                            printf("%.4s can donate %d H-Bond%c\n", necatm[i]+4*j, donors[i][j], (donors[i][j]==1)?' ':'s');
                        for(j=0;strncmp(necatm[i]+4*j,"    ",4)&&j<TOTNAA;j++)
                          if (accepts[i][j])
                            printf("%.4s can accept %d H-Bond%c\n", necatm[i]+4*j, accepts[i][j], (accepts[i][j]==1)?' ':'s');
                      }
                  }
                  break;
                    
                case 'v':
                case 'V':
                  if (i+1>=optc) fail("SYNTAX ERROR: Not Enough Arguments for -v/V option\n");
                  sscanf(optv[i+1],"%hd",&numcovbonds);
                  jflg=0;
                  i++;
                  break;
                case 'q':
                  asaflg=1;
                  break;
                case 'Q':
                  asaflg=1;
                  exchangeflg=2;
                  break;
                    
#ifdef DOM
                case 'Y':
                  if (i+4>=optc) fail("SYNTAX ERROR: Not Enough Arguments for -Y option\n");
                  sscanf(optv[i+1],"%f",&PEPBNDVAL);
                  sscanf(optv[i+2],"%f",&HBVAL);
                  sscanf(optv[i+3],"%f",&VDVVAL);
                  sscanf(optv[i+4],"%f",&LTLDOMBON);
                    
                  jflg=0;
                  i+=4;
                  break;
                    
                case 'Z':
                  domainflg=1;
                  break;
#endif
#ifdef BSM
                  /*v3.01:begin*/
                case 'W':
                  if (i+1>=optc) fail("SYNTAX ERROR: Not Enough Arguments for -Y option\n");
                  sscanf(optv[i+1],"%hd",&bsmoptflg);
                  printf("BSMOPTFLG set to %d\n",bsmoptflg);
                    
                  jflg=0;
                  i++;
                  break;
                  /*v3.01:end*/                    

#endif
                  /*v3.13:begin . . . susan wants to run HBPLUS to produce 'h' 
                    replacing the insertion code for hetatms*/
                case 'm':
                    
                  OMLINSERT = 1;
                  break;
                  /*v3.13:end*/
                default:
                  printf("Error - invalid option %s.\n",optv[i]);
                  continue;
                }
            }
        }
      else
        {
          if(*tmpstr)
            strcpy(inpdbfn,optv[i]);
          else
            strcpy(tmpstr,optv[i]);
        }
    }
}

void parse_options_file(char * filename)
{
     FILE * optfp;
     char buf[256],*linev[128],*hash;
     int linec;
/* Amendment. RAL 9 Mar 2009 --> */
     char wkdir[FILENAME_LEN];
/* <-- RAL 9 Mar 2009 */

/* Amendment. RAL 14 May 1997 --> */
     int done, in_quote;
     int ipos;
/* <-- RAL 14 May 1997 */
/* Amendment. RAL 25 May 1997 --> */
     int line;
/* <-- RAL 25 May 1997 */

     printf("Opening options file %s.\n",filename);
     optfp=fmfopen(filename,"r");
     while(fgets(buf,256,optfp) && !feof(optfp))
     {
/* Amendment. RAL 14 May 1997 --> */
/*         if(hash=strchr(buf,'#'))
             *hash='\0'; */ /*allow # to be used for comments*/

       /* Initialise flags */
       done = FALSE;
       in_quote = FALSE;

       /* Loop through all the characters in the line to check for
          hashes representing a comment */
       for (ipos = 0; ipos < strlen(buf) && done == FALSE; ipos++)
         {
           /* If this is a double-quote, then check whether it is a
              start- or end-quote */
           if (buf[ipos] == '"')
             {
               if (in_quote == FALSE)
                 in_quote = TRUE;
               else
                 in_quote = FALSE;
             }

           /* If this is a hash, check whether it is inside or outside
              a string */
           else if (buf[ipos] == '#')
             {
               /* If outside a string, then terminate the line here */
               if (in_quote == FALSE)
                 {
                   buf[ipos] = '\0';
                   done = TRUE;
                 }
             }

           /* If end of string reached, then end here */
           else if (buf[ipos] == '\0' || buf[ipos] == '\n')
             done = TRUE;
         }
/* <-- RAL 14 May 1997 */
         linec=parse_line(buf,linev);

/* Amendment. RAL 25 May 1997 --> */
       /* Replace any @ symbols inserted by hbadd to signify double-quotes
          in atom-names */
       for (line = 0; line < linec; line++)
         {
           for (ipos = 0; ipos < (int) strlen(linev[line]); ipos++)
             if (linev[line][ipos] == '@')
               linev[line][ipos] = '"';
         }
/* <-- RAL 25 May 1997 */

/* Amendment. RAL 9 Mar 2009 --> */
         parse_options(linec,linev,0,wkdir);
/* <-- RAL 9 Mar 2009 */
     }
     fclose(optfp);
 }

/******************************************************************************/
/* Criteria - minimum angles at the Hydrogen and the Acceptor
            - maximum distances from the Hydrogen and from the Donor to the 
              Acceptor */

int main(int argc, char * argv[])
{
    int             i;
    char            calendar_s[128];
/* Amendment. RAL 24 Apr 2008 --> */
    char *err;
/* <-- RAL 24 Apr 2008 */
/* Amendment. RAL 9 Mar 2009 --> */
    char file_name[FILENAME_LEN], wkdir[FILENAME_LEN];

    int len;
/* <-- RAL 9 Mar 2009 */

    time_t    calendar;
    struct tm * calendar_tm_p;
    debug=0;

/* Amendment. RAL 3 Jul 2012 --> */
    tmpstr[0] = inpdbfn[0] = '\0';
/* <-- RAL 3 Jul 2012 */

    time ( &calendar );
    calendar_tm_p = localtime ( &calendar );
    strftime(calendar_s,128,"%b %d %X %Z %Y\0",calendar_tm_p);
    
    initialise_arrays();
    supplement_arrays();
#ifdef __unix
{
    char rcfn[256];
    FILE * rcfp;
    
/*    printf("Home is %s\n",getenv("HOME"));*/
/* Amendment. RAL 24 Apr 2008 --> */
    err = getenv("HOME");
    if (err != NULL)
      {

/* <-- RAL 24 Apr 2008 */
        strcpy(rcfn,getenv("HOME"));

/*    printf(".hbplusrc file at %s\n",rcfn);*/
        strcat(rcfn,"/");
/*    printf(".hbplusrc file at %s\n",rcfn);*/
        strcat(rcfn,".hbplusrc");
/*    printf(".hbplusrc file at %s\n",rcfn);*/
    
        rcfp = fopen(rcfn,"r"); /* check if the hbplusrc file is there.*/
        if (rcfp)
          {
            fclose(rcfp);
            
            parse_options_file(rcfn);
          }
/* Amendment. RAL 24 Apr 2008 --> */
      }
/* <-- RAL 24 Apr 2008 */
}

#endif
    
    printf("HBPLUS Hydrogen Bond Calculator v %s            %s\n",VERSION,calendar_s );
    printf("(c) I McDonald, D Naylor, D Jones and J Thornton 1993 All Rights Reserved.\n");
    
    printf("\nConfigured for %5d atoms and %5d residues.\n",MAXNATM , MAXNRES );
    
/* Amendment. RAL 9 Mar 2009 --> */
    parse_options(argc,argv,1,wkdir);
    hbdn[0] = '\0';
    if (wkdir[0] != '\0')
      strcpy(hbdn,wkdir);
/* <-- RAL 9 Mar 2009 */

/* Amendment. RAL 9 Mar 2009 --> */
    // Form name of hbdebug.dat file
    strcpy(file_name,"hbdebug.dat");
    if (wkdir[0] != '\0')
      {
        strcpy(file_name,wkdir);
        len = strlen(file_name);
        if (len > 0 && (file_name[len - 1] == '/' || file_name[len - 1] == '\\'))
          {
            // Add file separator
            strcat(file_name,"/");
          }
        strcat(file_name,"hbdebug.dat");
      }

    // Open the hbdebug.dat file
    dbgfp = fopen(file_name, "w");
/* <-- RAL 9 Mar 2009 */
    

    printf("\nCriteria\n\n");
    if(nnbflg)
        printf("Neighbour interactions, max atom-atom separation %3.1f\n", HBDIST);
    else
    {
        printf(  "Minimum Angles; DHA %5.2f, HAAA %5.2f, DAAA %5.2f\nMaximum Distances; D-A %3.1f, H-A %3.1f, S-S %3.1f\n",(mindha),(minhaaa),(mindaaa),HBDIST, MAX_HA, SSDIST);
        printf(  "Maximum angles at aromatic acceptors DAAX %5.2f, HAAAX %5.2f\n",maxdaax, maxhaax);
    }
    printf(  "Minimum covalent separation %d Covalent bonds\n",numcovbonds+1);
    if(cssflg)
        printf("CYS converted to CYH and CSS depending on presence of disulphides.\n");
    if(exchangeflg)
        printf("Alternative locations of His, Asn and Gln permissible\n");
    if(exchangeflg==2)
        printf("Only H-Bonds with His, Asn or Gln side-chains considered\n");
printf("\n");
        
if(debug) printf("About to process files %s and %s\n",tmpstr,inpdbfn);

    if(!*tmpstr)
    {
        printf("Enter directory name for output file (default current directory)\n");
/* Amendment. RAL 3 Jul 2012 --> */
//        gets(tmpstr);
        fgets(tmpstr,FILENAME_LEN,stdin);
        string_truncate(tmpstr,FILENAME_LEN);
/* <-- RAL 3 Jul 2012 */
        sscanf(tmpstr, "%s", hbdn);
        
        if (hbdn[strlen(hbdn) - 1] != '/' && *hbdn)
            strcat(hbdn,"/\0");
        
        for (;;)
        {
            printf("\nEnter the name of the next file to be processed:\n");
            printf("just a file name for a Brookhaven Format PDB file or\n");
            printf("@file name for a file holding names of Brookhaven");
            printf("format files - or a blank line to end.\n\n");
/* Amendment. RAL 3 Jul 2012 --> */
//            gets(tmpstr);
            fgets(tmpstr,FILENAME_LEN,stdin);
              string_truncate(tmpstr,FILENAME_LEN);
/* <-- RAL 3 Jul 2012 */
            if (tmpstr[0] == '\0' || tmpstr[0] == '\n')
                break;
            procinfilename(tmpstr,inpdbfn);
        }
    }
    else
        procinfilename(tmpstr,inpdbfn);
    

#ifndef BSM
/* Amendment. RAL 9 Mar 2009 --> 
    remove("hbdebug.dat");
 <-- RAL 9 Mar 2009 */
#endif
    return(0) ;
    
}


