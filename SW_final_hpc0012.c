/*-------------------------------------------------------------

Copyright (C) 2000 Peter Clote. 
All Rights Reserved.

Permission to use, copy, modify, and distribute this
software and its documentation for NON-COMMERCIAL purposes
and without fee is hereby granted provided that this
copyright notice appears in all copies.


THE AUTHOR AND PUBLISHER MAKE NO REPRESENTATIONS OR
WARRANTIES ABOUT THE SUITABILITY OF THE SOFTWARE, EITHER
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
PARTICULAR PURPOSE, OR NON-INFRINGEMENT. THE AUTHORS
AND PUBLISHER SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED
BY LICENSEE AS A RESULT OF USING, MODIFYING OR DISTRIBUTING
THIS SOFTWARE OR ITS DERIVATIVES.

-------------------------------------------------------------*/


/*************************************************
	Program: smithWaterman.c
	Peter Clote, 11 Oct 2000

Program for local sequence alignment, using the Smith-Waterman
algorithm and assuming a LINEAR gap penalty.
A traceback is used to determine the alignment, and
is determined by choosing that direction, for
which S(i-1,j-1)+sigma(a_i,b_j), S(i-1,j)+Delta and 
S(i,j-1)+Delta is maximum, i.e.  for which 

                    _
                   |
                   | H(i-1,j-1) + sigma(a_i,b_j)  (diagonal)
H(i,j) =  MAX of   | H(i-1,j)   + delta           (up)
                   | H(i,j-1)   + delta           (left)
                   | 0
                   | _


is a maximum.

*************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>   	// character handling
#include <stdlib.h>     // def of RAND_MAX 
#include "mpi.h"
   /* Just a note:                       */
   /* N must be the size of the arrays   */
   /* Here is assume that the two arrays */
   /* have the same size                 */


#define MAX_SEQ 50

#define CHECK_NULL(_check) {\
   if ((_check)==NULL) \
      fprintf(stderr, "Null Pointer allocating memory\n");\
   }

#define AA 20           // number of amino acids
#define MAX2(x,y)     ((x)<(y) ? (y) : (x))
#define MAX3(x,y,z)   (MAX2(x,y)<(z) ? (z) : MAX2(x,y))

// function prototypes
void error(char *);		/** error handling */

int char2AAmem[256];
int AA2charmem[AA];
void initChar2AATranslation(void);
int getRowCount(int rowsTotal, int mpiRank, int mpiSize);


int
main(int argc, char *argv[]) {

	// variable declarations
	FILE * in1, *in2, *pam;
	char ch;
	int temp;
	int i,j,jj,tempi,tempj,x,y,diag,down,right,DELTA;
	int topskip,bottomskip;
	char *aout,*bout;
	int Aend,Bend,Abegin,Bbegin;
	int max, Max, xMax, yMax;	
		// Max is first found maximum in similarity matrix H
		// max is auxilliary to determine largest of
		// diag,down,right, xMax,yMax are h-coord of Max
	short *a, *b;
	int *hptr;
	int **h;
	int sim[AA][AA];		// PAM similarity matrix
        int N;
        int nc;

        int B; // block size, as an argument from command line

	int rank, size;
	MPI_Status status;

	MPI_Init( &argc, &argv );
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int N_local;
	int N_sq;
	int N_sq2;


	/**** Error handling for input file ****/	
	if (( argc != 6) && (argc!=7)) {
	     	fprintf(stderr,"%s protein1 protein2 PAM gapPenalty [N]\n",argv[0]);
		exit(1);
	}
        else if (argc==6)  /* Maximum size of the proteins, they should have equal size */
        {
	   /***** Initialization of input file and pair array **/
	    in1   = fopen(argv[1],"r");
    	    in2   = fopen(argv[2],"r");
	    pam   = fopen(argv[3],"r");
	    DELTA = atoi(argv[4]);
	    B = atoi(argv[5]);
            N = MAX_SEQ;
        } else 
        {
	    in1   = fopen(argv[1],"r");
    	    in2   = fopen(argv[2],"r");
	    pam   = fopen(argv[3],"r");
	    DELTA = atoi(argv[4]);
	    B = atoi(argv[5]);
            N     = atoi(argv[6]);
        }

	N_local = getRowCount(N, rank, size);
	N_sq = (N+1)*(N+1);
	N_sq2 = (N+1)*(N_local+1);

        CHECK_NULL((aout = (char *) malloc(sizeof(char)*2*N)));
        CHECK_NULL((bout = (char *) malloc(sizeof(char)*2*N)));
        CHECK_NULL((a = (short *) malloc(sizeof(short)*(N+1))));
        CHECK_NULL((b = (short *) malloc(sizeof(short)*(N+1))));
        CHECK_NULL((hptr = (int *) malloc(sizeof(int)*(rank ? N_sq2 : N_sq) )));
        CHECK_NULL((h = (int **) malloc(sizeof(int*)*(rank ? N_local+1 : N+1) )));
        /* Mount h[N][N] */
        for(i=0;i<=N_local;i++) 
           h[i]=hptr+i*(N+1);


        initChar2AATranslation();

	/** read PAM250 similarity matrix **/	
        fscanf(pam,"%*s"); 
	for (i=0;i<AA;i++)
		for (j=0;j<=i;j++) {
		if (fscanf(pam, "%d ", &temp) == EOF) {
			fprintf(stderr, "PAM file empty\n");
			fclose(pam);
			exit(1);
			}
		sim[i][j]=temp;
		}
	fclose(pam);
	for (i=0;i<AA;i++)
		for (j=i+1;j<AA;j++) 
			sim[i][j]=sim[j][i]; 	// symmetrify

	if (!rank){
		/** read first file in array "a" **/	
	       
		i=0;
		do {
		   nc=fscanf(in1,"%c",&ch);
		   if (nc>0 && char2AAmem[ch]>=0) 
		   {
		      a[++i] = char2AAmem[ch]; 
		   }
		} while (nc>0 && (i<N));
		a[0]=i;  
		fclose(in1);
	}

	MPI_Bcast(a, N+1, MPI_SHORT, 0, MPI_COMM_WORLD);

	if (!rank){

		/** read second file in array "b" **/	
		i=0;
		do {
		   nc=fscanf(in2,"%c",&ch);
		   if (nc>0 && char2AAmem[ch]>=0) 
		   {
		      b[++i] = char2AAmem[ch]; 
		   }
		} while (nc>0 && (i<N));
		b[0]=i;  
		fclose(in2);

	}

	
	MPI_Bcast(b, N+1, MPI_SHORT, 0, MPI_COMM_WORLD);

	/** compute "h" local similarity array **/
	for (i=0;i<=N_local;i++) h[i][0]=0;
	if (!rank) for (j=0;j<=N;j++) h[0][j]=0;
	
	int bMax = 0;

	for (jj=1;jj<=N;jj+=B){
		if (rank){
			MPI_Recv(h[0]+jj, B, MPI_INT, rank - 1, 1, MPI_COMM_WORLD, &status);
		}
		for (i=1; i<=N_local; i++){
			for (j=jj;j<jj+B;j++) {
				diag    = h[i-1][j-1] + sim[a[i]][b[j]];
				down    = h[i-1][j] + DELTA;
				right   = h[i][j-1] + DELTA;
				max=MAX3(diag,down,right);
				if (max <= 0)  {
					h[i][j]=0;
					}
				else if (max == diag) {
					h[i][j]=diag;
					}
				else if (max == down) {
					h[i][j]=down;
					}
				else  {
					h[i][j]=right;
					}
				if (max > bMax){
					bMax=max;
					xMax=i;
					yMax=j;
					}
			}
		}
		if (rank < size - 1){
			MPI_Send(h[N_local]+jj, B, MPI_INT, rank + 1, 1, MPI_COMM_WORLD);	
		}
	}
	
	int finalMax;
	MPI_Reduce(&bMax, &finalMax, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

	int localmax[2];
	int globalmax[2];

	localmax[0] = bMax;
	localmax[1] = rank;

	MPI_Allreduce(localmax, globalmax, 1, MPI_2INT, MPI_MAXLOC, MPI_COMM_WORLD);

	int position[2];

	if (!rank && globalmax[1]==0){
		position[0]=xMax;
		position[1]=yMax;
	}
	
	if (rank==globalmax[1] && globalmax[1]!=0){
		position[0] = xMax;
		position[1] = yMax;
		MPI_Send(position, 1, MPI_2INT, 0, 2, MPI_COMM_WORLD);
	}
	if (!rank && globalmax[1]!=0){
		MPI_Recv(position, 1, MPI_2INT, globalmax[1], 2, MPI_COMM_WORLD, &status);
	}


	if (rank==0) {
		printf("Max score: %d \n", finalMax);
		printf("Rank max score: %d \n", globalmax[1]);
		printf("Position (i,j): (%d,%d) \n",position[0],position[1]);
	}

        MPI_Finalize( );
        return 0; 

	
}



void error(char * s) {
	fprintf(stderr,"%s\n",s);
	exit(1);
}

int getRowCount(int rowsTotal, int mpiRank, int mpiSize){
		return (rowsTotal / mpiSize) + (rowsTotal % mpiSize > mpiRank);
	}

void initChar2AATranslation(void)
{
    int i; 
    for(i=0; i<256; i++) char2AAmem[i]=-1;
    char2AAmem['c']=char2AAmem['C']=0;
    AA2charmem[0]='c';
    char2AAmem['g']=char2AAmem['G']=1;
    AA2charmem[1]='g';
    char2AAmem['p']=char2AAmem['P']=2;
    AA2charmem[2]='p';
    char2AAmem['s']=char2AAmem['S']=3;
    AA2charmem[3]='s';
    char2AAmem['a']=char2AAmem['A']=4;
    AA2charmem[4]='a';
    char2AAmem['t']=char2AAmem['T']=5;
    AA2charmem[5]='t';
    char2AAmem['d']=char2AAmem['D']=6;
    AA2charmem[6]='d';
    char2AAmem['e']=char2AAmem['E']=7;
    AA2charmem[7]='e';
    char2AAmem['n']=char2AAmem['N']=8;
    AA2charmem[8]='n';
    char2AAmem['q']=char2AAmem['Q']=9;
    AA2charmem[9]='q';
    char2AAmem['h']=char2AAmem['H']=10;
    AA2charmem[10]='h';
    char2AAmem['k']=char2AAmem['K']=11;
    AA2charmem[11]='k';
    char2AAmem['r']=char2AAmem['R']=12;
    AA2charmem[12]='r';
    char2AAmem['v']=char2AAmem['V']=13;
    AA2charmem[13]='v';
    char2AAmem['m']=char2AAmem['M']=14;
    AA2charmem[14]='m';
    char2AAmem['i']=char2AAmem['I']=15;
    AA2charmem[15]='i';
    char2AAmem['l']=char2AAmem['L']=16;
    AA2charmem[16]='l';
    char2AAmem['f']=char2AAmem['F']=17;
    AA2charmem[17]='L';
    char2AAmem['y']=char2AAmem['Y']=18;
    AA2charmem[18]='y';
    char2AAmem['w']=char2AAmem['W']=19;
    AA2charmem[19]='w';
}


