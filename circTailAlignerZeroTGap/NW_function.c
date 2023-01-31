
#define MATCH      1
#define MISMATCH  -3
#define INDEL     -100
// I have set INDEL penalty to 100 to prevent from any indel in the 
// Alignment, however I allow for indels in the beginnings of sequences
// by definding a separate START_INDEL penalty
#define START_INDEL -1



void NW(char* ReadsFile, char* refFile, int tr1)
{

// char      *x,*y;		/* the two sequences, x_1..M and y_1..N */
  int        M,N;               /* lengths of sequences x,y             */
  int        i,j;               /* residue coords in x and y            */
  int        S[500][500];       /* dynamic programming matrix           */
  int        BT[500][500];      /*Backtrack matrix                      */
  int        MM;
  
		
  char       ax[1000];          /* result: the two aligned sequences    */
  char       ay[1000];          /* result: the two aligned sequences    */
  int        K;			        /* length of pairwise alignment         */
  int        k;			        /* index in pairwise alignment, 0..K-1  */
  int        sc;                /* tmp variable used to find best path  */
  char       tmp;		        /* tmp variable used to swap chars      */
  /*****************************************************************
   * Some initial bookkeeping.
   *****************************************************************/
  /* Initialize the sequences x,y, using SEQX and SEQY.
   *
   * In C, arrays are indexed 0..N-1, but we want our sequences to be
   * indexed 1..N, and we need row/column 0 in the DP matrix for
   * initialization conditions. That's why we do some +1 fiddling
   * here, so our sequences are x_1..x_M and y_1..y_N.
   */
  char x[500];
  char y[500];
  char SEQY[500];  
  char Junk[500];
  
  char* filename = refFile;
  FILE *in_file = fopen(filename, "r");
  fscanf(in_file, "%[^\n] ", Junk);
  fscanf(in_file, "%[^\n] ", SEQY);
 // PosSpecificPenCalculator(SEQY);
  fclose(in_file);
  
  N = strlen(SEQY);
  strcpy(y+1, SEQY);		    /* y_1..y_N now defined; y_0 undefined */
  char ID[500];
  char SEQX[500];
  
  FILE *reads_file = fopen(ReadsFile, "r");
  char line_contents[500];
  int lineStat = 0;
  while (fscanf(reads_file, "%[^\n] ", line_contents) != EOF) {
		if (lineStat ==0){
			strcpy(ID, line_contents +1);
		}
		///////////////###################////////////
        else{ 
		    strcpy(SEQX, line_contents);
		
 
  M = strlen(SEQX);
  strcpy(x+1, SEQX);		    /* x_1..x_M now defined; x_0 undefined */


  /*****************************************************************
   * Recursion: the heart of the DP algorithm.
   *****************************************************************/
  /* Initialization.
   */
   
  /* BT symbols: 3 means the alignment is done
  just jump to the start position
  2 means diagonal movement and putting one nucleotide from each sequence
  1 means horizontal move
  0 means vertical move
  */
  S[0][0] = 0;
  BT[0][0] = 3;
  for (i = 1; i <= M; i++){
	  S[i][0] = i * START_INDEL;
	  BT[i][0] = 3;
  }
  
  for (j = 1; j <= N; j++){
	  
	  S[0][j] = j * START_INDEL;
	  BT[0][j] = 3;
	  }

  /* The dynamic programming, global alignment (Needleman/Wunsch) recursion.
   */
   
int end_onSeq1 = 0;
int end_onSeq2 = 0;
int maxByNow = -100;
int Tgap;
for (i = 1; i <= M; i++)
{
	
    for (j = 1; j <= N; j++)
    {
		Tgap = 0;
		
		MM = MISMATCH;
		
        /* case #1: i,j are aligned */
        if (x[i] == y[j]) S[i][j] = S[i-1][j-1] + MATCH;
        else S[i][j] = S[i-1][j-1] + MM;
		//printf("this nucleotide is being added from read %c\n this one is added from template %c\n this is the mismatch penalty %d\n", x[i], y[j], MM);
		BT[i][j] = 2;
		

        if(y[j] == 'T') sc = S[i][j-1] + Tgap;
		else sc=S[i][j-1] + INDEL; 
        /* case #3: j aligned to -  */
        if (sc > S[i][j]) {
			S[i][j] = sc;
			BT[i][j]= 1;
			}
		
		if(x[i] == 'T') sc = S[i-1][j] + Tgap;
		else sc = S[i-1][j] + INDEL; 
        /* case #2: i aligned to -  */
	    if (sc > S[i][j]){
			S[i][j] = sc;
			BT[i][j]= 0;
			}
			
		sc=(i+j)*START_INDEL;
		if (sc > S[i][j]){
			S[i][j] = sc;
			BT[i][j]= 3;
			}
			
		if (S[i][j] >= maxByNow){
			end_onSeq1 = i-1;
            end_onSeq2 = j-1;
            maxByNow = S[i][j];	
		}
      }
}	  

/* The result (optimal alignment score) is now in S[M][N]. */

    
int start_onSeq1 = 1;
int start_onSeq2 = 1;

  /* Allocate for our aligned strings (which will be 0..K-1)*/
  

  /* Start at M,N; work backwards */
  
  i = end_onSeq1 + 1;
  j = end_onSeq2 + 1;
  k = 0;
  while (2==2)
    {
      /* Case 1: best was x_i aligned to x_j?
       */
	   
	   if (BT[i][j]==2){
		   ax[k]= x[i];
           ay[k]= y[j];
           i--;
           j--;
		   k+=1;
			
	   }
	   else if (BT[i][j]==3){
		   start_onSeq1 = i;
		   start_onSeq2 = j;
		   break;
	   }
	   else if (BT[i][j]==1){
		   ax[k] = '-';
           ay[k]= y[j];
           j--; 
		   k+=1;
	   }
	   else if (BT[i][j]==0){
		   ax[k] = x[i];
           ay[k]= '-';
           i--; 
		   k+=1;
	   }   
	}
	ax[k] = '\0';
	ay[k] = '\0';
  /* Done with the traceback.
   * Now ax[0..k-1] and ay[0..k-1] are aligned strings, but backwards;
   * and k is our alignment length.
   * Reverse the strings and we're done. 
   */
  K = k;			/* K = remember the alignment length */
//  printf("this is k %d\n", K);
  
  for (k = 0; k < K/2; k++)
    {	// a sly, efficient in-place reversal by swapping pairs of chars: 
      tmp = ax[K-k-1]; ax[K-k-1] = ax[k]; ax[k] = tmp; 
      tmp = ay[K-k-1]; ay[K-k-1] = ay[k]; ay[k] = tmp; 
    }


  /*****************************************************************
   * Output
   *****************************************************************
   */
   
//  printf("Sequence X: %s\n", x+1);
//  printf("Sequence Y: %s\n", y+1);
/*
  printf("\n");
  
  printf("Dynamic programming matrix:\n");
  for (j = -1; j <= N; j++)
    printf("   %c ", j > 0 ? y[j] : ' ');
  printf("\n");
  for (i = 0; i <= M; i++)
    {
      printf("   %c ", i > 0 ? x[i] : ' ');
      for (j = 0; j <= N; j++)
	printf("%4d ", S[i][j]);
      printf("\n");
    }
  printf("\n\n\n\n");
*/

  printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n", ID, start_onSeq1+1, end_onSeq1+1, start_onSeq2+1, end_onSeq2+1, maxByNow, M, N, ax, ay);
  
//  printf("Alignment:\n");
//  printf("X: %s\n", ax);
//  printf("Y: %s\n", ay);

//  printf("Optimum alignment score: %d\n\n", maxByNow);

  /* Free the memory we allocated; then exit.
   */
    }
		lineStat = (lineStat +1)%2;
    }
  exit(0);
}