extern int quiet;
#define logmsg(level,...) if(level >= quiet) { fprintf(stderr,__VA_ARGS__); }

/**
 * Read a line of arbitrary length from a file, returning a MALLOCed
 * null-terminated string with the \n already stripped.
 */
char * np_getl(FILE* f);




/**
 * This function will join a directory base and filename semi-intelligently.
 *
 *
 * Returns a malloced null-terminated string, should be freed by caller.
 */
char* join_paths(char* base, char* filename);


/**
 * mapfile opens a whole file as an mmap segment, RO, shared
 */
int mapfile(char* filename, char** addr, size_t* length);
