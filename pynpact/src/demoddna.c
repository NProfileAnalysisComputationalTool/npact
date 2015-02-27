#include <stdio.h>
#include <stdlib.h>


#include "util.h"

const char usage[]= "\nUsage: <input.ddna>\n";

int main(int argc, char *argv[]) {
  char* filename = "file.ddna";
  char* origin;
  size_t length;
  int i = 1;

  if(argc < 2) {
    fprintf(stderr, "ERROR: not enough arguments");
    fprintf(stderr, usage);
    exit(1);
  }
  filename = argv[i++];
  if( (i = mapfile(filename, &origin, &length)) ) {
    exit(i);
  }

  /**
   * In here references to origin[n] refer to the nth base (0 indexed)
   * of the genome.
   */

  exit(munmap(origin, length));
}
