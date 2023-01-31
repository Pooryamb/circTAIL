
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>

#include "NW_function.c"


int main(int argc, char *argv[])
{

	
	char* ReadsFile = argv[1];
	char* refFile = argv[2];
	int tr1 = atoi(argv[3]);
//	int tr2 = atoi(argv[4]);
//	int tr3 = atoi(argv[5]);
	
	NW(ReadsFile, refFile, tr1);	
	
    
}

