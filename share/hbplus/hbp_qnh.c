/* HBPLUS Amide and Histidine Routines */
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
 
/*This includes input routines for accessibilties*/
#include <string.h>  
#include "hbplus.h"
  
short int asaflg=0;
void getcoord(float *x,float * y,float * z,float * occ, float *b, char * chain,int * n,int * aacode,char * resnam,char * atmnam,int * atmnum,char * buf);
enum CONFCLASS 
{ 
    u_u,u_i,u_s,i_i,i_s,s_s
};
enum SCCLASS
{ 
    h_sus, s_sus, indiff, s_opt, h_opt
};

int qn_tot[5], his_tot[5];

int inasa_file(char * asafname)
{
    FILE *afp;
    char buf[256],chain,resnam[4],atmnam[5],altcode,inscode;
    int n,aacode,atmnum,atomcount,thisnum,lastnum;
    
    float f_null;
    char c_null, *s_null;
    int i_null;
    
    fopen(asafname,"r");

    printf("Processing asa file\"%s\"\n", asafname);
    afp = fopen(asafname, "r");
    if (!afp)
    {
        printf("Failed to open specified input file %s!\n",asafname);
        return(0);
    }
    else
	printf("Opened asa file.\n");
    
    lastnum= -1;
    
    while(!feof(afp))
    {
	if (!fgets(buf, 160, afp))
            break;
	if (!strncmp("ATOM  ",buf,6)||!strncmp("HETATM",buf,6))
	{
	    getcoord(&f_null,&f_null,&f_null,&f_null,&f_null,
		     &chain,&n,&aacode,resnam,atmnam,&atmnum,buf);

	    altcode = buf[16];
            if (altcode == ' ')
                altcode = '-';

            inscode = buf[26];
            if (inscode == ' ')
                inscode = '-';

	    thisnum=find_atom4(chain,n,inscode,resnam,atmnam,lastnum+1);
	    if(thisnum>-99)
	    {
		lastnum=thisnum;
		sscanf(buf+65,"%f",&atom[thisnum].acc);
	    }
	    
	}
    }
    return(lastnum);
    
}/*end of the 'reading asa' */

void printf_numhb()
{
    int i;
    char buf[80];
    
    for(i=0;i<MAXNATM;i++)
    {
	if(atom[i].aacode<MAXNAA)
	{
	    if( donors[ atom[i].aacode ][ atom[i].atmtyp ] )
		printf("%s %2d%7.3f\n",atomid(i,buf),atom[i].ndonhb,atom[i].acc);
	    if( accepts[ atom[i].aacode ][ atom[i].atmtyp ] )
		printf("%s %2d%7.3f\n",atomid(i,buf),atom[i].nacchb,atom[i].acc);
	}
    }
}

struct qnhresidue 
{
    struct pdbatm * a[4];
    struct pdbres * r;
} qnhres[MAXNRES];

int initialise_qnhres(void)
{
    int aacode,i;
    struct pdbatm * a; /* A pointer to the first atom in the residue*/
    int qnhkount=0;
    int j,k,donnum,metnum ; /* Added v3.08, loop variables to look for interresidue covalent bonds to His side-chains and the donor/metal atom[] numbers */
    char o_str[128]; /* Storage for output line for metal ion pseudo-HBonds */
    
       
    for(i=0;i<reskount;i++)
    {
	a= &atom[residue[i].strtat];
	
	aacode=a->aacode;
	
	if(aacode==Asn||aacode==Gln||aacode==His)
	{
	    qnhres[qnhkount].r= &residue[i];
	    
	    if(aacode==Asn)
	    {
/*		printf("Asn %c%d%c found\n",a->chnid,a->aanum,a->inscode);
*/		
		qnhres[qnhkount].a[0]= &atom[ find_atom4(a->chnid, a->aanum, a->inscode, a->resnam, " OD1", residue[i].stopat) ];
		qnhres[qnhkount].a[1]= &atom[ find_atom4(a->chnid, a->aanum, a->inscode, a->resnam, " ND2", residue[i].stopat) ];
	    }		
		
	    if(aacode==Gln)
	    {
		qnhres[qnhkount].a[0]= &atom[ find_atom4(a->chnid, a->aanum, a->inscode, a->resnam, " OE1", residue[i].stopat) ];
		qnhres[qnhkount].a[1]= &atom[ find_atom4(a->chnid, a->aanum, a->inscode, a->resnam, " NE2", residue[i].stopat) ];
	    }		
		
	    if(aacode==His)
	    {
/*		printf(". . . identified as Histidine (%d==%d)\n",aacode,His);
	*/	
		qnhres[qnhkount].a[0]= &atom[ find_atom4(a->chnid, a->aanum, a->inscode, a->resnam, " ND1", residue[i].stopat) ];
		qnhres[qnhkount].a[1]= &atom[ find_atom4(a->chnid, a->aanum, a->inscode, a->resnam, " CD2", residue[i].stopat) ];
		qnhres[qnhkount].a[2]= &atom[ find_atom4(a->chnid, a->aanum, a->inscode, a->resnam, " CE1", residue[i].stopat) ];
		qnhres[qnhkount].a[3]= &atom[ find_atom4(a->chnid, a->aanum, a->inscode, a->resnam, " NE2", residue[i].stopat) ];
		for(j=0;j<4;j++)
		{
		    donnum = qnhres[qnhkount].a[j] - atom;
		    
		    for(k=0;k<atom[donnum].ncon ; k++)
		    {
			metnum= icon[donnum][k];
			
			if(a->caindex != atom[ metnum ].caindex)
			{
			    /* If there is a covalent bond between the N and a different residue then count as a donated H-bond. */ 
			    sprintf_hb(o_str, donnum, metnum, vector_of(0,0,0), length_squared( to(atom[donnum].p,atom[metnum].p)),-1.0,-1.0,-1.0);
			    /* add to the H-bond list . . . . NB: takes 'd^2', not 'd'. */
			    fprintf_hb(o_str, donnum, metnum);
			}
			
		    }
		}
			    
	    }		

	    qnhkount++;
	}
    }
    return(qnhkount);
}

int class_his(struct pdbatm *nd, struct pdbatm *ne)
{
    int i;
    if( (nd->ndonhb || nd->nacchb) && (ne->ndonhb || ne->nacchb))
    {
	if( !nd->ndonhb && nd->nacchb && !ne->ndonhb && ne->nacchb)
	{
	    if( nd->acc<=0 && ne->acc<= 0)
		return(u_s);
	    else
		return(i_s);
	}
	else
	    return(s_s);
    }
    else
    {
	int numhb=0;
	if( nd->ndonhb || nd->nacchb )
	    numhb+=100;
	else
	    if( nd->acc > 0.0)
		numhb+=80;

	if( ne->ndonhb || ne->nacchb )
	    numhb+=100;
	else
	    if( ne->acc > 0.0)
		numhb+=80;
	
	if(numhb==0)
	    return(u_u);
	if(numhb==80)
	    return(u_i);
	if(numhb==100)
	    return(u_s);
	if(numhb==160)
	    return(i_i);
	if(numhb==180)
	    return(i_s);
    }
    printf("ERROR: No classification of %c%d%cHis set\n",nd->chnid,nd->aanum,nd->inscode);
    return(-1);
    
}
    
/* Amendment. RAL 14 Jun 2012 --> */
//class_qn(struct pdbatm *o, struct pdbatm *n)
int class_qn(struct pdbatm *o, struct pdbatm *n)
/* <-- Amendment. RAL 14 Jun 2012 */
{
    int numhb=0;
/*    printf("Oxygen: don %d acc %d slv %5.1f\n",o->ndonhb,o->nacchb,o->acc);
    printf("Nitrgn: don %d acc %d slv %5.1f\n",n->ndonhb,n->nacchb,n->acc);
 */   
    
    if( o->nacchb > 0 )
    {
	/*printf("Oxygen accepts real HB - +100\n");*/
	numhb+=100;
    }
    else
	if( o->acc >0.0 )
	{
	    /*printf("Oxygen is exposed +80\n");*/
	    numhb+=80;
	}
    
    if( n->ndonhb > 0)					
    {							
	/*printf("Nitrogen donates real HB +100\n");*/	
	numhb+=100;					
    }
    else
	if( n->acc >0.0 )
	{
	    /*printf("Nitrogen is exposed +80\n");*/
	    numhb+=80;
	}
    
    if(numhb==0)
	return(u_u);
    if(numhb==80)
	return(u_i);
    if(numhb==100)
	return(u_s);
    if(numhb==160)
	return(i_i);
    if(numhb==180)
	return(i_s);
    if(numhb==200)
	return(s_s);
    printf("BUG: numhb %d unrecognised in class_qn\n",numhb);
    return(-99);
}


class_qn_tiebrake(struct pdbatm *o, struct pdbatm *n)
{
    int numhb=0;
/*    printf("Oxygen: don %d acc %d slv %5.1f\n",o->ndonhb,o->nacchb,o->acc);
    printf("Nitrgn: don %d acc %d slv %5.1f\n",n->ndonhb,n->nacchb,n->acc);
 */   
    
    if( o->nacchb > 1 )
    {
	/*printf("Oxygen accepts real HB - +100\n");*/
	numhb+=100;
    }
    else
	if( o->acc >0.0 )
	{
	    /*printf("Oxygen is exposed +80\n");*/
	    numhb+=80;
	}
    
    if( n->ndonhb > 1)					
    {							
	/*printf("Nitrogen donates real HB +100\n");*/	
	numhb+=100;					
    }
    else
	if( n->acc >0.0 )
	{
	    /*printf("Nitrogen is exposed +80\n");*/
	    numhb+=80;
	}
    
    if(numhb==0)
	return(u_u);
    if(numhb==80)
	return(u_i);
    if(numhb==100)
	return(u_s);
    if(numhb==160)
	return(i_i);
    if(numhb==180)
	return(i_s);
    if(numhb==200)
	return(s_s);
    printf("BUG: numhb %d unrecognised in class_qn\n",numhb);
    return(-99);
}


/* Amendment. RAL 14 Jun 2012 --> */
//compare_class(int pdb, int alt)
int compare_class(int pdb, int alt)
/* <-- Amendment. RAL 14 Jun 2012 */
{
    if(pdb==alt)
	return(indiff);
    if(pdb==u_u)
	return(h_sus);
    if((pdb==u_i || pdb==u_s) && (alt==i_i||alt==i_s||alt==s_s))
	return(h_sus);
    if(pdb<alt)
	return(s_sus);
    if(alt==u_u)
	return(h_opt);
    if((alt==u_i || alt==u_s) && (pdb==i_i||pdb==i_s||pdb==s_s))
	return(h_opt);
    if(alt<pdb)
	return(s_opt);
    printf("BUG: No classification in comp_class\n");
    return(-99);
}


void do_qnh(int qnh)
{
    char buf[256];
    int i,j,num_intatms,pdb,alt;
    int longoutflg_old;
    int result;
    
    short qnhgufo_flg=1;
    
    if( qnhres[qnh].a[0]->aacode == His )
	num_intatms=4;
    else
	num_intatms=2;
        
    if( qnhres[qnh].a[0]->aacode == His )
    {
	pdb=class_his(qnhres[qnh].a[0],qnhres[qnh].a[3]);
	alt=class_his(qnhres[qnh].a[1],qnhres[qnh].a[2]);
	result=compare_class(pdb,alt);
	
    }
    	else
    {
	pdb=class_qn(qnhres[qnh].a[0],qnhres[qnh].a[1]);
	alt=class_qn(qnhres[qnh].a[1],qnhres[qnh].a[0]);
	result=compare_class(pdb,alt);
	if(result==indiff)
	{
	    
	    result=compare_class(class_qn_tiebrake(qnhres[qnh].a[0],qnhres[qnh].a[1]),class_qn_tiebrake(qnhres[qnh].a[1],qnhres[qnh].a[0]));
	    if(result==h_sus)
		result=s_sus;
	    if(result==h_opt)
		result=s_opt;
	}
	
    }    

    if( !qnhgufo_flg )
    {
	printf("%s,",atomid(qnhres[qnh].r->strtat,buf));
	for(i=0;i<num_intatms;i++)
	    printf("%d,%d,",qnhres[qnh].a[i]->ndonhb,qnhres[qnh].a[i]->nacchb);
	for(i=0;i<num_intatms;i++)
	    if( qnhres[qnh].a[i]->acc > 0 )
		printf("SLV O  ,SLV O  ,");
	    else
		printf(",,");
	printf("%d,%d,%d",pdb,alt,result);
	
	printf("\n");
	
    }
    else
    {
	longoutflg_old=longoutflg;
	longoutflg=0;
	
	printf("\nRESIDUE %.14s : ",atomid((int)(qnhres[qnh].a[0] - atom),buf));
        switch( result )
	{
	case h_sus:
	    printf("HIGHLY SUSPECT");
	    break;
	case s_sus:
	    printf("SLIGHTLY SUSPECT");
	    break;
	case indiff:
	    printf("INDIFFERENT");
	    break;
	case s_opt:
	    printf("SLIGHTLY OPTIMAL");
	    break;
	case h_opt:
	    printf("HIGHLY OPTIMAL");
	    break;
	}
	printf("\n");
	printf("==============================================================================\n");
	printf("\nAccessibilities\n");
	printf("---------------\n");
	printf("Atom\t\t\tAccessibility\n");
	
	for(i=0;i<num_intatms;i++)
	{
	    printf("%s\t%5.2f\n",atomid((int)(qnhres[qnh].a[i] - atom),buf),qnhres[qnh].a[i]->acc );
	    *buf=0;
	}
	printf("\nHydrogen Bonds\n");
	printf("--------------\n");
	printf("<---DONOR---> <-ACCEPTOR--> D-A     gap   CC   DHA   H-A H-A-AA D-A-AA\n");
	
	for(i=0;i<num_intatms;i++)
	{
	    for(j=0;j<qnhres[qnh].a[i]->ndonhb;j++)
	    {
		/*printf("i=%d, j=%d,",i,j);
		  printf("DON to j=%s\n", atomid(qnhres[qnh].a[i]->don_p->n[j], buf));
		  *buf=0;*/
		
		
		if( calc_hb( buf, (int)(qnhres[qnh].a[i] - atom ), qnhres[qnh].a[i]->don_p->n[j] ) )
		    printf("%s\n",buf);
		else
		    if( alreadybonded( (int)(qnhres[qnh].a[i] - atom ), qnhres[qnh].a[i]->don_p->n[j]) )
			/* This is an inserted check v3.08 for His N-metal covalent bonds, in which case . . */
		    {
			sprintf_hb(buf,  (int)(qnhres[qnh].a[i] - atom ), qnhres[qnh].a[i]->don_p->n[j], vector_of(0,0,0), length_squared( to( qnhres[qnh].a[i]->p , atom[ qnhres[qnh].a[i]->don_p->n[j] ].p )),-1.0,-1.0,-1.0);
			printf("%s\n",buf);
		    }
		    else
			printf("BUG in do_qnh DON for%6d,%6d SC# %d - no H-bond # %d where expected\n",(int)(qnhres[qnh].a[i] - atom) , qnhres[qnh].a[i]->don_p->n[j],i,j);
		*buf=0;
		
	    }
	

	    for(j=0;j<qnhres[qnh].a[i]->nacchb;j++)
	    {
		/*printf("i=%d, j=%d",i,j);
		  printf(",ACC from j=%s\n", atomid(qnhres[qnh].a[i]->acc_p->n[j], buf));
		  *buf=0;*/
	       
		
		if( calc_hb( buf, qnhres[qnh].a[i]->acc_p->n[j],(int)(qnhres[qnh].a[i] - atom ) ) )
		    printf("%s\n",buf);
		else
		    printf("BUG in do_qnh ACC for%6d,%6d SC# %d - no H-bond # %d where expected\n", qnhres[qnh].a[i]->acc_p->n[j],(int)(qnhres[qnh].a[i] - atom ), i,j);
		*buf=0;
		debug=0;
		
	    }

	}
	
	longoutflg=longoutflg_old;
	
	
    }
    
    if( qnhres[qnh].a[0]->aacode == His )
	his_tot[result]++;
    else
	qn_tot[result]++;
    
    
}

void chkqnh(void)
{
    int qnhreskount,i;
/*    char buf[256]; */
    
    qnhreskount=initialise_qnhres();
    for(i=0;i<qnhreskount;i++)
    {
	
/*	printf("Checking i=%d - %s\n",i,atomid(qnhres[i].r->strtat,buf));
	printf("%c%d%c %s%s:%d: %d : His %d Asn %d Gln %d\n",qnhres[i].a[0]->chnid,qnhres[i].a[0]->aanum,qnhres[i].a[0]->inscode,qnhres[i].a[0]->resnam,qnhres[i].a[0]->atmnam, (int)(qnhres[i].a[0]-atom),qnhres[i].a[0]->aacode, His, Asn, Gln);*/
       
	do_qnh(i);
	
	
    }
#ifdef BSM
    printf("\nSUMMARY:%.4s,%d,%d,%d,%d,%d,-UNSAT-,%d,",brcode,qn_tot[4],qn_tot[3],qn_tot[2],qn_tot[1],qn_tot[0],qn_tot[0]+qn_tot[1]+qn_tot[2]+qn_tot[3]+qn_tot[4]);
    printf("%d,%d,%d,%d,%d,-UNSAT-,%d,\n",his_tot[4],his_tot[3],his_tot[2],his_tot[1],his_tot[0],his_tot[0]+his_tot[1]+his_tot[2]+his_tot[3]+his_tot[4]);
#endif
}

