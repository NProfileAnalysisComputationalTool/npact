/*  -*- c-file-style:"linux" c-basic-offset:3 tab-width:3 -*-  */
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <math.h>
#include <string.h>

#define ORDER 0

const char usage[]= "\nUsage: CG input.gbk [start end [window_size step [period_of_frames]]]\nDefaults: start= 1; end= end of genome sequence; window_size= 201; step= 51; period_of_frames= 3.\n";


int mapfile(char* filename, char** addr, size_t* length) {
   struct stat sb;
   int fd;

   fd = open(filename, O_RDONLY);
   if (fd == -1) {
      fprintf(stderr, "ERROR: Error opening '%s': ", filename);
      perror(NULL);
      return 1;
   }

   if (fstat (fd, &sb) == -1) {
      perror("fstat");
      return 1;
   }

   if (!S_ISREG (sb.st_mode)) {
      fprintf (stderr, "%s is not a file\n", filename);
      return 1;
   }
   *length = sb.st_size;
   *addr = mmap (NULL, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);
   if (*addr == MAP_FAILED) {
      perror ("mmap");
      return 1;
   }

   if (close (fd) == -1) {
      perror ("close");
      return 1;
   }
}

int calculateProfile(const char* origin, const size_t length,
                     int start, int end, int window, int step, int period) {
   char *box, base;
   int offset, i, j, baseidx=0, n=0;
   double	*S,k;

   fprintf(stderr,
           "Reading bases %d-%d with window %d nt, step %d nt and period %d nt.\n",
           start, end, window, step, period);

   box = (char *)malloc(window*sizeof(char));
   S = (double *)calloc(period, sizeof(double));
   for(i=0; i<window; ++i) box[i]= ' ';

   //while(fgets(longstr,198,input) && !feof(input) && (start+n)<=end) {
   for(offset = 0; offset < length && start+n <= end; offset++) {
      base = origin[offset];
      if(base >= 'a' && base <= 'z') {
         ++baseidx;
         if(baseidx >= start && baseidx <= end) {
            if(base=='c' || base=='g') {
               if(box[n % window] != 'S') {
                  ++S[(baseidx-1)%(period)]; box[n%window]='S';
               }
               ++n;
            }
            else if(base=='a' || base=='t' || base=='u') {
               if(box[n%window]=='S') {
                  --S[(baseidx-1)%(period)]; box[n%window]='W';
               }
               ++n;
            }
            else {
               if(box[n%window]=='S') {
                  --S[(baseidx-1)%(period)];
                  box[n%window]='N';
               }
               fprintf(stderr,"\nBase %c found at position %d\n",base,start+n);
               ++n;
            }
            if(n>=window && !((n-window)%step) && (start+n)<end) {
               fprintf(stdout,"%8d",start+n-1-window/2);
               for(j=0;j<period;++j)
                  fprintf(stdout,"%8.1f", 100.0/(float)(window/period)*(float)S[j]);
               fprintf(stdout,"\n");
            }
         }
      }
   }
   fprintf(stderr,"Bases %d-%d read (%d nt)\n", start, end, n);
}

int countBases(const char* origin, const size_t length) {
   /* how many bases from there */
   int i, tot=0;
   for(i=0; i < length; i++)
      if(origin[i] >= 'a' && origin[i] <= 'z')
         ++tot;
   return tot;
}

main(argc,argv)
int	argc;
char	*argv[];
{
   int start,end,window=201,step=51,period=3,tot=0;
   char* mapped_file;
   char* origin;
   size_t length, olength;


   if(!(argc==2 || argc==4 || argc==6 || argc==7)) {
      fputs(usage, stderr);
      exit(1);
   }

   if ( mapfile(argv[1], &mapped_file, &length) )
      exit(1);


   fprintf(stderr, "Searching for coding sequence in %s. ", argv[1]);
   /* Skip until we see ORIGIN (the line that starts the coding seq)
    * If the file doesn't contain ORIGIN this may raise a segfault;
    * but at that point error is about appropriate anyways.
    */
   origin = strstr(mapped_file, "\nORIGIN");
   /* then saqve a pointer at the start of the next line */
   origin = strchr(origin + strlen("\nORIGIN"), '\n') + 1;
   olength = length - (origin - mapped_file);

   //while(fgets(longstr,198,input) && strncmp(longstr,"ORIGIN",6));

   tot = countBases(origin, olength);
   fprintf(stderr, "Found %d bases.\n", tot);

   if(argc >= 4) { start= atoi(argv[2]); end= atoi(argv[3]); }
   else { start= 1; end= tot; }

   if(argc >= 6) {
      window = atoi(argv[4]);
      step = atoi(argv[5]);
   }
   if(argc == 7)
      period = atoi(argv[6]);



   if(end-start+1 < 0 || end-start+1 > tot) {
      fprintf(stderr,"ERROR: Sequence must be longer than 0 and shorter than complete sequence (%d nt).\n%s",
              tot,usage);
      exit(1);
   }

   if(window % period) {
      fprintf(stderr,"ERROR: Window_size must be divisable by period_of_frames\n%s", usage);
      exit(1);
   }

   calculateProfile(origin, olength, start, end, window, step, period);

   munmap(mapped_file, length);


   exit(0);
}
