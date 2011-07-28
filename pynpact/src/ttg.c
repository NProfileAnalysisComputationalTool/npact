# include <stdio.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>

// Direct strand start codons:
// atg= 1 * 0 + 5 * 3 + 25 * 2 = 65
//          a       t        g
// gtg= 1 * 2 + 5 * 3 + 25 * 2 = 67
//          g       t        g
// ttg= 1 * 3 + 5 * 3 + 25 * 2 = 68
//          t       t        g

// Complementary strand start codons:
// atg= 1 * 1 + 5 * 0 + 25 * 3 = 76
//          c       a        t
// gtg= 1 * 1 + 5 * 0 + 25 * 1 = 26
//          c       a        c
// ttg= 1 * 1 + 5 * 0 + 25 * 0 = 1
//          c       a        a

int numbers[]= { 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4 };

main(int argc, char *argv[])   // argv[1] is name of GenBank file
{
  int	n= 0, cod= 124, MYCOPLASMA;  // initializes codon to 'nnn'
char	c, longstr[200];
FILE	*input;

        if(argc == 1) { fprintf(stderr, "\nUsage: ttg input.gbk [ 0 / 1 (default: 0, not a Mycoplasma)]"); exit(1); }
  
        if(argc == 3) MYCOPLASMA= atoi(argv[2]);

	if((input= fopen(argv[1], "r")) == NULL) { fprintf(stderr, "\nInput file %s not found.\n", argv[1]); exit(1); }

	do fgets(longstr, 198, input); while(strncmp(longstr, "ORIGIN", 6) && !feof(input)); 

	if(feof(input)) { fprintf(stderr, "\nOrigin not found.\n"); exit(1); }

	while(!feof(input))
	{
	c= fgetc(input);
		if(c >= 'A' && c <= 'Z') c += 'a' - 'A';
		if(c >= 'a' && c <= 'z')
		{
		cod /= 5;
		cod += 25 * numbers[(int)(c - 'a')];
		++n;
			if(cod == 68) fprintf(stdout, "D %d\n", n - 2);
			else if(cod == 1) fprintf(stdout, "C %d\n", n);
		}
	}	

fclose(input);
}
