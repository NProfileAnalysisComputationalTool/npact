# include <stdlib.h>
# include <stdio.h>
# include <string.h>

/*
 * Read a line of arbitrary length from a file, returning a malloced
 * null-terminated string with the \n already stripped.
 */
char * np_getl(FILE * f) {
    size_t size = 0,last= 0, len = 0;
    char * buf  = NULL;
    if(feof(f))
        return NULL;
    do {
        /* BUFSIZ is defined as "the optimal read size for this platform" */
        size += BUFSIZ;
        /* realloc(NULL,n) is the same as malloc(n) */
        buf = (char*) realloc(buf,size);
        /* Actually do the read. Note that fgets puts a terminal '\0'
           on the end of the string, so we make sure we overwrite
           this */
        fgets(buf+len, BUFSIZ, f);
        len += strlen(buf+len);
        last = len - 1;
    } while (!feof(f) && buf[last] != '\n');

    if(buf[last] == '\n') {
        /* get rid of the newline */
        buf[last] = '\0';
        return (char*) realloc(buf,len);
    }
    else {
        return (char*) realloc(buf,len+1);
    }
}
