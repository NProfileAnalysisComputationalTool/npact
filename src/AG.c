# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <string.h>

# define ORDER 0

char usage[]= "\nUsage: AG input.gbk [start end [window_size step [period_of_frames]]]\nDefaults: start= 1; end= end of genome sequence; window_size= 201; step= 51; period_of_frames= 3.\n";

main(argc,argv)
int	argc;
char	*argv[];
{
int	i,j,n=0,m=0,start,end,window=201,step=51,period=3,tot=0;
char	longstr[200],filename[200],*box;
double	*S,k;
FILE	*input;

  if(!(argc==2 || argc==4 || argc==6 || argc==7)) { fprintf(stderr,usage); exit(1); }

  if((input= fopen(argv[1],"r"))==NULL) { fprintf(stderr,"\nERROR: Input file %s not found%s",argv[1],usage); exit(1); }

while(fgets(longstr,198,input) && strncmp(longstr,"ORIGIN",6));
  while(fgets(longstr,198,input) && !feof(input))
    for(i=0;i<strlen(longstr);++i) if(longstr[i]>='a' && longstr[i]<='z') ++tot;
fclose(input);

input= fopen(argv[1],"r");

  if(argc>=4) { start= atoi(argv[2]); end= atoi(argv[3]); }
  else { start= 1; end= tot; }

  if(argc>=6) { window= atoi(argv[4]); step= atoi(argv[5]); }
  if(argc==7) period= atoi(argv[6]);

S= (double *)calloc(period,sizeof(double));

  if(end-start+1<0 || end-start+1>tot) { fprintf(stderr,"\nERROR: Sequence must be longer than 0 and shorter than complete sequence (%d nt).\n%s",tot,usage); exit(1); }
  if(window%period) { fprintf(stderr,"\nERROR: Window_size must be divisable by period_of_frames\n%s",usage); exit(1); }

fprintf(stderr,"\nReading file %s, bases %d-%d with window %d nt, step %d nt and period %d nt.\n",argv[1],start,end,window,step,period);

while(fgets(longstr,198,input) && strncmp(longstr,"ORIGIN",6));

box= (char *)malloc(window*sizeof(char));

  for(i=0;i<window;++i) box[i]= ' ';

  while(fgets(longstr,198,input) && !feof(input) && (start+n)<=end)
  {
    for(i=0;i<strlen(longstr)-1 && (start+n)<=end;++i)
    {
      if(longstr[i]>='a' && longstr[i]<='z')
      {
      ++m;
        if(m>=start && m<=end)
        {
          if(longstr[i]=='a' || longstr[i]=='g') { if(box[n%window]!='S') { ++S[(m-1)%(period)]; box[n%window]='S'; } ++n; }
          else if(longstr[i]=='c' || longstr[i]=='t' || longstr[i]=='u') { if(box[n%window]=='S') { --S[(m-1)%(period)]; box[n%window]='W'; } ++n; }
          else { if(box[n%window]=='S') { --S[(m-1)%(period)]; box[n%window]='N'; } fprintf(stderr,"\nBase %c found at position %d\n",longstr[i],start+n); ++n; }
          if(n>=window && !((n-window)%step) && (start+n)<end)
          {
          fprintf(stdout,"%8d",start+n-1-window/2);
            for(j=0;j<period;++j) fprintf(stdout,"%8.1f",100.0/(float)(window/period)*(float)S[j]);
          fprintf(stdout,"\n");
          }
        }
      }
    }
  }

fclose(input);

fprintf(stderr,"\nbases %d-%d read (%d nt)\n",start,end,n);
}
