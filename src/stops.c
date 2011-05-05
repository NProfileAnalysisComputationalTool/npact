# include <stdio.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>

// Direct strand stop codons:
// taa= 1 * 3 + 5 * 0 + 25 * 0 = 3
//          a       t        g
// tag= 1 * 3 + 5 * 0 + 25 * 2 = 53
//          g       t        g
// tga= 1 * 3 + 5 * 2 + 25 * 0 = 13   if MYCOPLASMA == 0
//          t       t        g

// Complementary strand stop codons:
// tta= 1 * 3 + 5 * 3 + 25 * 0 = 18
//          c       a        t
// cta= 1 * 1 + 5 * 3 + 25 * 0 = 16
//          c       a        c
// tca= 1 * 3 + 5 * 1 + 25 * 0 = 8
//          c       a        a

int numbers[]= { 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4 };

main(int argc, char *argv[])   // argv[1] is name of GenBank file
{
int	n= 0, cod= 124, MYCOPLASMA= 0;  // initializes codon to 'nnn'
char	c, longstr[200];
FILE	*input;

	if(argc == 1) { fprintf(stderr, "\nUsage: stops input.gbk [ 0 / 1 (default: 0, not a Mycoplasma)]"); exit(1); }

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
			if(cod == 3 || cod == 53 || (!MYCOPLASMA && cod == 13)) fprintf(stdout, "D %d\n", n - 2);
			else if(cod == 18 || cod == 16 || (!MYCOPLASMA && cod == 8)) fprintf(stdout, "C %d\n", n);
		}
	}	

fclose(input);
}
