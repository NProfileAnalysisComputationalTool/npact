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

const char usage[]= "\nUsage: [-b bases] input.gbk [start end [window_size step [period_of_frames]]]\nDefaults: start=1; end= end of genome sequence; window_size= 201; step= 51; period_of_frames= 3.\n";


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
   return 0;
}

void calculateProfile(const char* origin,
                      const char* bases, const int start, const int end,
                      const int window, const int step, const int period) {

   char *box, base;
   int j, idx, n=0, baseidx=0;
   size_t *S;
   double normfactor = 100.0/(double)(window/period);

   fprintf(stderr,
           "Calculating [%s]-profile of bases %d-%d with window %d nt, step %d nt and period %d nt.\n",
           bases, start+1, end+1, window, step, period);

   box = (char *) calloc(window, sizeof(char));
   S = (size_t *) calloc(period, sizeof(size_t));

   for(idx=start; idx <= end; idx++) {
      base = origin[idx];
      baseidx = idx+1; // this is the 1 based index we use for printing.
      if(strchr(bases, base)) {
         if(! box[n % window]) {
            //if the flag wasn't set, then set it and increase the count.
            ++S[(idx) % period];
            box[n % window] = 1;
         }
      }
      else {
         if(box[n % window]) {
            //if the flag was set, clear it and decrease the count.
            --S[(idx) % period];
            box[n % window] = 0;
         }
         //unexpected base, print warning.
         if(! strchr("ACTGU", base))
            fprintf(stderr, "\nBase %c found at position %d\n", base, baseidx);
      }
      ++n;

      if(n >= window && ((n-window) % step) == 0 && idx <= end) {
         fprintf(stdout, "%8d", baseidx-window/2);
         for(j=0; j < period; ++j)
            fprintf(stdout, "%8.1f", normfactor*S[j]);
         fprintf(stdout,"\n");
      }
   }
   free(box);
   free(S);
   fprintf(stderr, "Bases %d-%d read (%d nt)\n", start, end, n);
}

int main(int argc, char *argv[]) {
   int start=0,end=0,window=201,step=51,period=3;
   int i;
   char* origin;
   size_t length;
   int argi = 1;
   char* filename;
   char* bases = "CG";

   if(argc < 2) {
      fprintf(stderr, "ERROR: not enough arguments");
      fprintf(stderr, usage);
      exit(1);
   }

   if(strcmp(argv[argi], "-b") == 0) {
      bases = (char*) argv[++argi];
      for (i = 0; bases[i]; i++) {
         if(! strchr("ACGT", bases[i])) {
            fprintf(stderr, "ERROR: must specify bases as subset of 'ACGT'\n");
         }
      }
      argi++;
   }
   filename = argv[argi++];
   if ( i = mapfile(filename, &origin, &length) )
      exit(i);


   if(argc > argi)
      start = atoi(argv[argi++]) - 1;

   end = (argc > argi) ? atoi(argv[argi++]) - 1 : length - 1;

   if(argc > argi)
      window = atoi(argv[argi++]);
   if(argc > argi)
      step = atoi(argv[argi++]);
   if(argc > argi)
      period = atoi(argv[argi++]);


   if(end >= length || start < 0 || end <= start) {
      fprintf(stderr, "ERROR: Sequence must be longer than 0 and shorter than complete sequence (%zd nt).\n%s",
              length, usage);
      exit(1);
   }
   if(window % period) {
      fprintf(stderr,"ERROR: Window_size must be divisable by period_of_frames\n%s", usage);
      exit(1);
   }

   calculateProfile(origin, bases, start, end, window, step, period);
   exit(munmap(origin, length));
}
