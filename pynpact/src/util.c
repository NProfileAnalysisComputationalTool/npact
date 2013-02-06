# include <stdlib.h>
# include <stdio.h>
# include <string.h>

#include "util.h"


/* For use with logmsg, controls how verbose the message stream is. */
int quiet = 0;

/*
 * Read a line of arbitrary length from a file; strips the \n off the end.
 *
 * Returns a malloced null-terminated string, should be freed by caller.
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


/**********************************/
/****  Function join_paths()   ****/
/**********************************/
/*
 * This function will join a directory base and filename semi-intelligently.
 *
 *
 * Returns a malloced null-terminated string, should be freed by caller.
 */
char* join_paths(char* base, char* filename) {
    char* returnfile = NULL;
    int len = strlen(filename) + 1;
    if(base) {
        // add an extra space in case we need a path separator.
        len += strlen(base) + 1;
    }

    returnfile = (char*) calloc(len, sizeof(char));
    //Start with the base
    if(base) {
        strcpy(returnfile, base);
        //Add a path separator if needed
        if (returnfile[strlen(returnfile)] != '/')
            strcat(returnfile, "/");
    }
    //Finally add the filename
    strcat(returnfile, filename);
    return returnfile;
}
