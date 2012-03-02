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

const char usage[]= "\nUsage: [-b bases] input.gbk [start end [window_size step [period_of_frames]]]\nDefaults: start= 1; end= end of genome sequence; window_size= 201; step= 51; period_of_frames= 3.\n";


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

void calculateProfile(const char* origin, const size_t length,
                      const char* bases, const int start, const int end,
                      const int window, const int step, const int period) {

   char *box, base;
   int offset =0 , j, baseidx=0, n=0;
   size_t *S;
   double normfactor = 100.0/(double)(window/period);

   fprintf(stderr,
           "Calculating [%s]-profile of bases %d-%d with window %d nt, step %d nt and period %d nt.\n",
           bases, start, end, window, step, period);

   box = (char *)calloc(window, sizeof(char));
   S = (size_t *)calloc(period, sizeof(size_t));


   /* skip up till start, but don't move offset passed start. */
   while (offset < length) {
      base = origin[offset];
      if(base >= 'a' && base <= 'z') {
         if(baseidx + 1 == start)
            break;
         baseidx++;
      }
      offset++;
   }

   for(; offset < length && baseidx < end; offset++) {
      base = origin[offset];
      if(base >= 'a' && base <= 'z') {
         ++baseidx;
         if(strchr(bases, base)) {
            if(! box[n % window]) {
               //if the flag wasn't set, then set it and increase the count.
               ++S[(baseidx-1) % period];
               box[n % window] = 1;
            }
         }
         else {
            if(box[n % window]) {
               //if the flag was set, clear it and decrease the count.
               --S[(baseidx-1) % period];
               box[n % window] = 0;
            }

            //unexpected base, print warning.
            if(! strchr("actgu", base))
               fprintf(stderr,"\nBase %c found at position %d\n",base,baseidx);
         }
         ++n;

         if(n >= window && ((n-window) % step) == 0 && baseidx < end) {
            fprintf(stdout, "%8d", baseidx-window/2);
            for(j=0; j < period; ++j)
               fprintf(stdout, "%8.1f", normfactor*S[j]);
            fprintf(stdout,"\n");
         }
      }
   }
   free(box);
   free(S);
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

int main(int argc, char *argv[]) {
   int start=1,end=0,window=201,step=51,period=3,tot=0;
   int i;
   char* mapped_file;
   char* origin;
   size_t length, olength;
   int argi = 1;
   char* filename;
   char* bases = "cg";


   if(!(argc==2 || argc==4 || argc==6 || argc==7)) {
      fputs(usage, stderr);
      exit(1);
   }
   if(strcmp(argv[argi], "-b") == 0) {
      bases = (char*) argv[++argi];
      for (i = 0; bases[i]; i++)
         if (bases[i] < 'a')
            bases[i] = bases[i] + 'a' - 'A';
      argi++;
   }
   filename = argv[argi++];


   if ( mapfile(filename, &mapped_file, &length) )
      exit(1);


   fprintf(stderr, "Searching for coding sequence in '%s'... ", filename);
   /* Skip until we see ORIGIN (the line that starts the coding seq)
    * If the file doesn't contain ORIGIN this may raise a segfault;
    * but at that point exiting with error is all we can do anyways.
    */
   origin = strstr(mapped_file, "\nORIGIN");
   /* then save a pointer at the start of the next line */
   origin = strchr(origin + strlen("\nORIGIN"), '\n') + 1;
   olength = length - (origin - mapped_file);

   tot = countBases(origin, olength);
   fprintf(stderr, "Found %d bases.\n", tot);

   if(argc > argi)
      start = atoi(argv[argi++]);

   end = (argc > argi) ? atoi(argv[argi++]) : tot;

   if(argc > argi)
      window = atoi(argv[argi++]);
   if(argc > argi)
      step = atoi(argv[argi++]);
   if(argc > argi)
      period = atoi(argv[argi++]);


   if(end > tot || start < 1 || end <= start) {
      fprintf(stderr,"ERROR: Sequence must be longer than 0 and shorter than complete sequence (%d nt).\n%s",
              tot,usage);
      exit(1);
   }

   if(window % period) {
      fprintf(stderr,"ERROR: Window_size must be divisable by period_of_frames\n%s", usage);
      exit(1);
   }

   calculateProfile(origin, olength, bases, start, end, window, step, period);

   munmap(mapped_file, length);
   exit(0);
}
