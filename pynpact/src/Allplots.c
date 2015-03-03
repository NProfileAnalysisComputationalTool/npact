/*  -*- c-file-style:"linux" c-basic-offset:4 tab-width:4 -*-  */
# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

#include "util.h"

# define DICT_SIZE 500
# define WIDTH 450.0
# define TITLES_POSITION 780
# define FIRST_TRANS -160
# define TRANSLATE -135
# define HIGHT 50.0
# define HIGHT_PUB 28.0  // Published annotation ORFs
# define HIGHT_MOD 28.0  // Modified ORFs
# define HIGHT_NEW 50.0  // New ORFs (bold-face) and Potential new ORFs (roman and light)
# define HIGHT_BIA 22.0
# define HIGHT_CP 8.0
# define HIGHT_SCP 8.0
# define HIGHT_HOM 52.5
# define HIGHT_MET 10.0
# define HAIRPIN_RADIUS 2.0
# define STEM_SEPARATION 1.0
# define PRINT_ZERO_LINE 1

# define VIEWER "showps"
# define TIC_W 3.0
# define NFS 6.0
# define NaFS 5.0
# define LFS 9.0
# define TFS 14.0
# define TFS2 14.0
# define HANGOVER 30
# define DSHORT_OFFSET 0.5
# define CSHORT_OFFSET -3.5
# define DLONG_OFFSET 6.5
# define CLONG_OFFSET -9.5

char	ANOPIAS= 0;
float	MAX= 100.00;
float   MAXLOG= 4.00;
float   GRID_Y;
float	TIC_X;
float	TIC_Y;
int     EVERY= 2;

char NAME[100];
char NUCLEOTIDES[10];

void printHelp() {
    fprintf(stderr,"\nUsage:\n\n  Allplots start interval lines [x-tics] [period_of_frames] [end_base]\n");
    fprintf(stderr,"start                 Genome interval first base.\n");
    fprintf(stderr,"interval              Number of bases.\n");
    fprintf(stderr,"lines                 Number of lines on page (one page).\n");
    fprintf(stderr,"x-tics                Number of subdivisions.\n");
    fprintf(stderr,"period_of_frame       Number of frames.\n");
    fprintf(stderr,"end_base              Genome interval of final base.\n");

    fprintf(stderr,"\nNeeds imput file \"Allplots.def\" of the form:\n\n");
    fprintf(stderr,"<Plot_title> <gbkfile length>\n");
    fprintf(stderr,"Nucleotide(s)_plotted (e.g.: C+G)\n");
    fprintf(stderr,"First-page title\n");
    fprintf(stderr,"Title of following pages\n");
    fprintf(stderr,"File_of_unbiased_CDSs\n");
    fprintf(stderr,"File_of_conserved_CDSs\n");
    fprintf(stderr,"File_of_new_CDSs\n");
    fprintf(stderr,"File_of_potential_new_CDSs\n");
    fprintf(stderr,"File_of_stretches_where_CG_is_asymmetric\n");
    fprintf(stderr,"File_of_published_accepted_CDSs\n");
    fprintf(stderr,"File_of_published_rejected_CDSs\n");
    fprintf(stderr,"File_of_blocks_from_new_ORFs_as_cds\n");
    fprintf(stderr,"File_of_blocks_from_annotated_genes_as_cds\n");
    fprintf(stderr,"File_of_GeneMark_regions\n");
    fprintf(stderr,"File_of_G+C_coding_potential_regions\n");
    fprintf(stderr,"File_of_met_positions (e.g.:D 432)\n");
    fprintf(stderr,"File_of_stop_positions (e.g.:D 432)\n");
    fprintf(stderr,"File_of_tatabox_positions (e.g.:105.73 D 432 TATAAAAG)\n");
    fprintf(stderr,"File_of_capbox_positions\n");
    fprintf(stderr,"File_of_ccaatbox_positions\n");
    fprintf(stderr,"File_of_gcbox_positions\n");
    fprintf(stderr,"File_of_kozak_positions\n");
    fprintf(stderr,"File_of_palindrom_positions_and_size\n");
    fprintf(stderr,"File_list_of_nucleotides_in_200bp windows.\n");
    fprintf(stderr,"File_list_of_reads.\n");
    fprintf(stderr,"Allplots.def should be in the current directory or provided on stdin with the --stdin option\n");
}


/**************************************/
/****** Function name_yoffset() *******/
/**************************************/

// pos and len are in the scale of line length, where
// 50,000 bp equal 300 text positions.

int name_yoffset(int *line, int *cline, int pos, int len, char strand)
{
int     i, j, n, s, k, y[4], cy[4], ctot[4]= { 0 }, tot[4]= { 0 };

pos -= 1;
len += 2;

    if(strand == 'D')
    {
        for(i= 0; i < len; ++i)
        {
                for(j= 0; j < 4; ++j) y[j]= 0;
        s= line[pos + i];
        s %= 16;
        n= (int)pow(2.0, (double)(3));
                for(j= 3; j > 0; --j)
                {
                y[j]= s / n;
                s %= n;
                n /= 2;
                }
        y[0]= s;
                for(j= 0; j < 4; ++j) if(!y[j]) ++tot[j];
        }

n= 0;
        while(tot[n] < len) ++n;

        for(i= pos; i < pos + len; ++i) line[i] +=  (int)pow(2.0, (double)n);
    }
    else
    {
        for(i= 0; i < len; ++i)
        {
                for(j= 0; j < 4; ++j) cy[j]= 0;
        s= cline[pos + i];
        s %= 16;
        n= (int)pow(2.0, (double)(3));
                for(j= 3; j > 0; --j)
                {
                cy[j]= s / n;
                s %= n;
                n /= 2;
                }
        cy[0]= s;
                for(j= 0; j < 4; ++j) if(!cy[j]) ++ctot[j];
        }

n= 0;
        while(ctot[n] < len) ++n;

        for(i= pos; i < pos + len; ++i) cline[i] +=  (int)pow(2.0, (double)n);
    }

return(n);
}

/*********************************************/
/****** End of function name_yoffset() *******/
/*********************************************/


int main(int argc, char *argv[]) {
    /* For cmd line argument parsing. */
    char* opt;
    int argi = 1;
    char use_stdin=0;
    int tstart, bp_per_page, tend, lines, period=3;

    /* other counters, values, variables for the program */

    int	i, j, k, d, n, N, nub, nc, nn, nnP, np, ne, nb, ns, ncg, nB, ncp, nScp,
        nm, nk, nt, ncap, ncca, ngcb, gs, ge, lp, len, swflag=0,
        start, end, pos, name_pos, name_len, line_range,
        expand=0,
        unbf=0, conf=0, newf=0, newPf=0, cgf=0, cpf= 0, Scpf= 0, npali= 0, wind,
        pub_line[300 + 2 * HANGOVER], new_line[300 + 2 * HANGOVER],
        exc_line[300 + 2 * HANGOVER], newP_line[300 + 2 * HANGOVER],
        pub_cline[300 + 2 * HANGOVER], new_cline[300 + 2 * HANGOVER],
        exc_cline[300 + 2 * HANGOVER], newP_cline[300 + 2 * HANGOVER];
    float	y1, y2, r, tpr, delta, pp, S[3], cds_len;

    /* Filename buffers */
    char *unb_file, *con_file, * new_file, *newP_file, *cg_file, *pub_file,
        *mod_file, *block_file, *BLOCK_file, *codpot_file, *Hits_file, *met_file,
        *kozak_file, *tata_file, *pali_file, *cap_file, *ccaa_file, *gcbox_file,
        *stop_file, *CG200_file, *read_file;


    /* by declaring them initially to be NULL, when realloc sees it, it
       will treat it as a malloc */
    int **unb=NULL, **con=NULL, **new=NULL, **newP=NULL, **cg=NULL, **pub=NULL,
        **modified=NULL, **block=NULL, **BLOCK=NULL, **codpot=NULL, *Hits=NULL,
        *met=NULL, *tata=NULL, *cap=NULL, *ccaa=NULL, *gcbox=NULL, *stop=NULL,
        *kozak=NULL, **pali=NULL, *X=NULL, *x=NULL;
    float **y=NULL,**Y=NULL, **reads= NULL;

    char *unb_str=NULL, *con_str=NULL, *new_str=NULL, *newP_str=NULL, *cg_str=NULL,
        *pub_str=NULL, *mod_str=NULL, *block_str=NULL, *BLOCK_str=NULL,
        *codpot_str=NULL, *codpot_col=NULL, *Hits_str=NULL, *Hits_col=NULL, *Hits_type=NULL,
        *met_str=NULL, *tata_str=NULL, *cap_str=NULL, *ccaa_str=NULL,
        *gcbox_str=NULL, *stop_str=NULL, *kozak_str=NULL;

    char **pub_name=NULL,**mod_name=NULL,**new_name=NULL,**newP_name=NULL;

    char longstr[200], *p, tatastr[20], *title1, *title2, ts[2], pos_string[15];
    FILE *input, *files;

    /**** Parse command line options ****/
    while (argi < argc && argv[argi][0] == '-') {
        opt = argv[argi];
        if(strcmp(opt, "-q") == 0) {
            quiet += 10;
            argi++;
        }
        else if(strcmp(opt, "-C") == 0) {
            /* Active the alternate color mode */
            ANOPIAS = 1;
            argi++;
            logmsg(0,"Using alternate color scheme.\n");
        }
        else if(strcmp(opt, "--stdin") == 0) {
            /* read Allplots.def from stdin instead of from file */
            use_stdin = 1;
            argi++;
        }
        else if(strcmp(opt, "-h") == 0) {
            printHelp();
            exit(0);
        }
        else {
            logmsg(20,"Got invalid optional argument %s", opt)
        }
    }

    if (argi+3 > argc) {
        printHelp();
        logmsg(10, "Not enough cmdline arguments.\n");
        exit(1);
    }

    tstart = atoi(argv[argi++]);
    bp_per_page = atoi(argv[argi++]);
    lines = atoi(argv[argi++]);
    line_range = bp_per_page / lines;
    delta = (float)(line_range);

    if(argi < argc)
        TIC_X = atof(argv[argi++]);
    else
        TIC_X = (float)(line_range/10);

    if(argi < argc)
        period = atoi(argv[argi++]);

    /* What's the final base we should evaluate through */
    if(argi < argc)
        tend = atoi(argv[argi++]);
    else
        tend = tstart + bp_per_page;
    logmsg(0, "Using end value %d\n", tend);


    /****** Parse definition file ******/
    logmsg(10, "Starting read of Allplots.def\n");
    files = use_stdin ? stdin : fopen("Allplots.def","r");
    if(!files || feof(files)) {
        printHelp();
        logmsg(20, "Missing or empty Allplots.def.\n");
        exit(1);
    }

    fgets(longstr,98,files);

    sscanf(longstr,"%s %d", NAME, &len);
    logmsg(10, "%s %d nt\n", NAME, len);

    fgets(NUCLEOTIDES,8,files);
    NUCLEOTIDES[strlen(NUCLEOTIDES)-1]='\0';
    logmsg(0, "%s\n", NUCLEOTIDES);

    title1 = np_getl(files);     logmsg(0,"\n%s",title1);
    title2 = np_getl(files);     logmsg(0,"%s\n",title2);
    /* File_of_unbiased_CDSs */
    unb_file = np_getl(files);   logmsg(0,"%s\n",unb_file);
    /* File_of_conserved_CDSs */
    con_file = np_getl(files);   logmsg(0,"%s\n",con_file);
    /* File_of_new_CDSs */
    new_file = np_getl(files);   logmsg(0,"%s\n",new_file);
    /* File_of_potential_new_CDSs */
    newP_file = np_getl(files);  logmsg(0,"%s\n",newP_file);
    /*File_of_stretches_where_CG_is_asymmetric */
    cg_file = np_getl(files);    logmsg(0,"%s\n",cg_file);
    /*File_of_published_accepted_CDSs */
    pub_file = np_getl(files);   logmsg(0,"%s\n",pub_file);
    /*File_of_published_rejected_CDSs */
    mod_file = np_getl(files);   logmsg(0,"%s\n",mod_file);
    /*File_of_blocks_from_new_ORFs_as_cds */
    block_file = np_getl(files); logmsg(0,"%s\n",block_file);
    /*File_of_blocks_from_annotated_genes_as_cds */
    BLOCK_file = np_getl(files); logmsg(0,"%s\n",BLOCK_file);
    /*File_of_GeneMark_regions */
    codpot_file = np_getl(files); logmsg(0,"%s\n",codpot_file);
    /*File_of_G+C_coding_potential_regions */
    Hits_file = np_getl(files); logmsg(0,"%s\n",Hits_file);
    /*File_of_met_positions (e.g.:D 432) */
    met_file = np_getl(files);   logmsg(0,"%s\n",met_file);
    /*File_of_stop_positions (e.g.:D 432) */
    stop_file = np_getl(files);  logmsg(0,"%s\n",stop_file);
    /*File_of_tatabox_positions (e.g.:105.73 D 432 TATAAAAG) */
    tata_file = np_getl(files);  logmsg(0,"%s\n",tata_file);
    /*File_of_capbox_positions */
    cap_file = np_getl(files);   logmsg(0,"%s\n",cap_file);
    /*File_of_ccaatbox_positions */
    ccaa_file = np_getl(files);  logmsg(0,"%s\n",ccaa_file);
    /*File_of_gcbox_positions */
    gcbox_file = np_getl(files); logmsg(0,"%s\n",gcbox_file);
    /*File_of_kozak_positions */
    kozak_file = np_getl(files); logmsg(0,"%s\n",kozak_file);
    /*File_of_palindrom_positions_and_size */
    pali_file = np_getl(files);  logmsg(0,"%s\n",pali_file);
    /*File_list_of_nucleotides_in_200bp windows. */
    CG200_file = np_getl(files); logmsg(0,"%s\n",CG200_file);
    /*File_list_of_read_numbers. */
    read_file = np_getl(files); logmsg(0,"%s\n",read_file);

    fclose(files);
    logmsg(20,"Done reading Allplots.def\n");


    logmsg(10, "\nData read from position %d to position %d.\n", tstart, tstart+line_range);
    logmsg(10, "Nucleotide frequencies computed with period %d\n", period);
    logmsg(10, "Plots of %d line per page, %d positions per line,tics every %.0f positions.\n",
           lines, line_range, TIC_X);


    /************ Start Working **********/

    for(k=0; k<lines && tstart + k*line_range<len; ++k) {
        n= N= nub= nc= nn= nnP= ncg= nm= nk= nt= ncap= ncca= ngcb= ns= np= ne= nb= nB= ncp= nScp= npali= 0;
        start= tstart + k*line_range;
        end= tstart + (k+1)*line_range;
        if(end > len) end = len;
        if(end > tend) end = tend;
        if(end < start) break;

// BLOCK_file and block_file excluded

/*
        if((input=fopen(BLOCK_file,"r"))) {
            logmsg(10, "Reading BLOCK_file %s\n", BLOCK_file);
            while(fgets(longstr,198,input) && !feof(input)) {
                BLOCK_str= (char *)realloc(BLOCK_str,(nB+1)*sizeof(char));
                BLOCK= (int **)realloc(BLOCK,(nB+1)*sizeof(int *));
                BLOCK[nB]= (int *)malloc(2*sizeof(int));
                p= strchr(longstr,'.');
                ge= atoi(p+2);
                if(longstr[0]=='c') { gs= atoi(longstr+11); BLOCK_str[nB]='C'; }
                else { gs= atoi(longstr); BLOCK_str[nB]='D'; }
                if(gs>=start-line_range/50 && gs<end && ge>start && ge<=end+line_range/50) { BLOCK[nB][0]= gs; BLOCK[nB][1]= ge; ++nB; }
                else if(gs>=start && gs<end && ge>end+line_range/50) { BLOCK[nB][0]= gs; BLOCK[nB][1]= end+line_range/50; ++nB; }
                else if(ge<=end+line_range/50 && ge>start && gs<start-line_range/50)
                { BLOCK[nB][0]= start-line_range/50; BLOCK[nB][0] += gs%period-BLOCK[nB][0]%period; BLOCK[nB][1]= ge; ++nB; }
                else if(gs<start-line_range/50 && ge>end+line_range/50)
                { BLOCK[nB][0]= start-line_range/50; BLOCK[nB][0] += gs%period-BLOCK[nB][0]%period; BLOCK[nB][1]= end+line_range/50; ++nB; }
            }
            fclose(input);
        }
        else logmsg(10, "File of blocks from annotated genes %s NOT read\n", BLOCK_file);

        if((input=fopen(block_file, "r"))) {
            logmsg(10, "Reading block_file %s\n", block_file);
            while(fgets(longstr,198,input) && !feof(input)) {
                block_str= (char *)realloc(block_str,(nb+1)*sizeof(char));
                block= (int **)realloc(block,(nb+1)*sizeof(int *));
                block[nb]= (int *)malloc(2*sizeof(int));
                p= strchr(longstr,'.');
                ge= atoi(p+2);
                if(longstr[0]=='c') { gs= atoi(longstr+11); block_str[nb]='C'; }
                else { gs= atoi(longstr); block_str[nb]='D'; }
                if(gs>=start-line_range/50 && gs<end && ge>start && ge<=end+line_range/50) { block[nb][0]= gs; block[nb][1]= ge; ++nb; }
                else if(gs>=start && gs<end && ge>end+line_range/50) { block[nb][0]= gs; block[nb][1]= end+line_range/50; ++nb; }
                else if(ge<=end+line_range/50 && ge>start && gs<start-line_range/50)
                { block[nb][0]= start-line_range/50; block[nb][0] += gs%period-block[nb][0]%period; block[nb][1]= ge; ++nb; }
                else if(gs<start-line_range/50 && ge>end+line_range/50)
                { block[nb][0]= start-line_range/50; block[nb][0] += gs%period-block[nb][0]%period; block[nb][1]= end+line_range/50; ++nb; }
            }
            fclose(input);
        }
        else logmsg(10, "File of new blocks %s NOT read\n", block_file);
*/

        if((input=fopen(codpot_file, "r"))) {
            logmsg(10, "Reading codpot_file %s\n", codpot_file);
            ++cpf;
            while(fgets(longstr,198,input) && !feof(input)) {
                codpot_str= (char *)realloc(codpot_str,(ncp+1)*sizeof(char));
                codpot_col= (char *)realloc(codpot_col,(ncp+1)*sizeof(char));
                codpot= (int **)realloc(codpot,(ncp+1)*sizeof(int *));
                codpot[ncp]= (int *)malloc(2*sizeof(int));
                p= strchr(longstr,'.');
                ge= atoi(p+2);
                if(longstr[0]=='c') { gs= atoi(longstr+11); codpot_str[ncp]='C'; codpot_col[ncp]= gs%period; }
                else { gs= atoi(longstr); codpot_str[ncp]='D'; codpot_col[ncp]= ge%period; }
                if(gs>=start && gs<end && ge>start && ge<=end) { codpot[ncp][0]= gs; codpot[ncp][1]= ge; ++ncp; }
                else if(gs>=start && gs<end && ge>end) { codpot[ncp][0]= gs; codpot[ncp][1]= end; ++ncp; }
                else if(ge<=end && ge>start && gs<start)
                { codpot[ncp][0]= start; codpot[ncp][0] += gs%period-codpot[ncp][0]%period; codpot[ncp][1]= ge; ++ncp; }
                else if(gs<start && ge>end)
                { codpot[ncp][0]= start; codpot[ncp][0] += gs%period-codpot[ncp][0]%period; codpot[ncp][1]= end; ++ncp; }
            }
            fclose(input);
        }
        else logmsg(10, "File of GeneMark coding potential %s NOT read\n", codpot_file);

        if((input=fopen(Hits_file,"r"))) {
            logmsg(10, "Reading Hits_file %s\n", Hits_file);
            ++Scpf;
            expand = 1;
            while(fgets(longstr,198,input) && !feof(input)) {
                if(expand) {
                    Hits_str= (char *)realloc(Hits_str, (nScp + 1) * sizeof(char));
                    Hits_col= (char *)realloc(Hits_col, (nScp + 1) * sizeof(char));
                    Hits_type= (char *)realloc(Hits_type, 2 * (nScp + 1) * sizeof(char));
                    Hits= (int *)realloc(Hits, 2 * (nScp + 1) * sizeof(int));
                }
                Hits_type[nScp * 2 + 0]= longstr[0];
                Hits_type[nScp * 2 + 1]= atoi(longstr + 1);
                p= strchr(longstr, '.');
                ge= atoi(p+2);
                if(longstr[3]=='c') { gs= atoi(longstr+14); Hits_str[nScp]='C'; Hits_col[nScp]= gs%period; }
                else if(longstr[3]=='r') { gs= atoi(longstr+10); Hits_str[nScp]='R'; }
                else { gs= atoi(longstr + 3); Hits_str[nScp]='D'; Hits_col[nScp]= ge%period; }
                expand = 1;   // Should we expand again next lop?
                if(gs>=start && gs<end && ge>start && ge<=end) { Hits[nScp * 2 + 0]= gs; Hits[nScp * 2 + 1]= ge; ++nScp; }
                else if(gs>=start && gs<end && ge>end) { Hits[nScp * 2 + 0]= gs; Hits[nScp * 2 + 1]= end; ++nScp; }
                else if(ge<=end && ge>start && gs<start)
                { Hits[nScp * 2 + 0]= start; Hits[nScp * 2 + 0] += gs % period - Hits[nScp * 2 + 0] % period; Hits[nScp * 2 + 1]= ge; ++nScp; }
                else if(gs<start && ge>end)
                { Hits[nScp * 2 + 0]= start; Hits[nScp * 2 + 0] += gs % period - Hits[nScp * 2 + 0] % period; Hits[nScp * 2 + 1]= end; ++nScp; }
                else {
                    expand = 0; // didn't use the current slot, don't need to expand
                }
            }
            fclose(input);
        }
        else logmsg(10, "File of G+C coding potential %s NOT read\n", Hits_file);

        if((input=fopen(met_file,"r"))) {
            logmsg(10, "Reading met_file %s\n", met_file);
            while(fgets(longstr,198,input) && !feof(input)) {
                gs= atoi(longstr+2);
                if(gs>=start && gs<=end)
                {
                    met_str= (char *)realloc(met_str,(nm+1)*sizeof(char));
                    met= (int *)realloc(met,(nm+1)*sizeof(int));
                    if(longstr[0]=='C') met_str[nm]='C';
                    else met_str[nm]='D';
                    met[nm]= gs;
                    ++nm;
                }
            }
            fclose(input);
        }
        else logmsg(10,"File of Met %s NOT read\n", met_file);

        if((input=fopen(stop_file,"r"))) {
            logmsg(10, "Reading stop_file %s\n", stop_file);
            while(fgets(longstr,198,input) && !feof(input)) {
                gs= atoi(longstr+2);
                if(gs>=start && gs<=end)
                {
                    stop_str= (char *)realloc(stop_str,(ns+1)*sizeof(char));
                    stop= (int *)realloc(stop,(ns+1)*sizeof(int));
                    if(longstr[0]=='C') stop_str[ns]='C';
                    else stop_str[ns]='D';
                    stop[ns]= gs;
                    ++ns;
                }
            }
            fclose(input);
        }
        else logmsg(10,"File of Stop %s NOT read\n", stop_file);


        if((input=fopen(tata_file,"r"))) {
            logmsg(10, "Reading tata_file %s\n", tata_file);
            while(fgets(longstr,198,input) && !feof(input)) {
                sscanf(longstr,"%f %s %d %s",&tpr,ts,&gs,tatastr);
                if(gs>=start && gs<=end) {
                    tata_str= (char *)realloc(tata_str,(nt+1)*sizeof(char));
                    tata= (int *)realloc(tata,(nt+1)*sizeof(int));
                    if(ts[0]=='C') tata_str[nt]='C';
                    else tata_str[nt]='D';
                    tata[nt]= gs;
                    ++nt;
                }
            }
            fclose(input);
        }
        else logmsg(10,"File of TATA box %s NOT read\n", tata_file);

        if((input=fopen(cap_file,"r"))) {
            logmsg(10, "Reading cap_file %s\n", cap_file);
            while(fgets(longstr,198,input) && !feof(input)) {
                sscanf(longstr,"%f %s %d %s",&tpr,ts,&gs,tatastr);
                if(gs>=start && gs<=end) {
                    cap_str= (char *)realloc(cap_str,(ncap+1)*sizeof(char));
                    cap= (int *)realloc(cap,(ncap+1)*sizeof(int));
                    if(ts[0]=='C') cap_str[ncap]='C';
                    else cap_str[ncap]='D';
                    cap[ncap]= gs;
                    ++ncap;
                }
            }
            fclose(input);
        }
        else logmsg(10,"File of CAP box %s NOT read\n", cap_file);

        if((input=fopen(ccaa_file,"r"))) {
            logmsg(10, "Reading ccaa_file %s\n", ccaa_file);
            while(fgets(longstr,198,input) && !feof(input)) {
                sscanf(longstr,"%f %s %d %s",&tpr,ts,&gs,tatastr);
                if(gs>=start && gs<=end) {
                    ccaa_str= (char *)realloc(ccaa_str,(ncca+1)*sizeof(char));
                    ccaa= (int *)realloc(ccaa,(ncca+1)*sizeof(int));
                    if(ts[0]=='C') ccaa_str[ncca]='C';
                    else ccaa_str[ncca]='D';
                    ccaa[ncca]= gs;
                    ++ncca;
                }
            }
            fclose(input);
        }
        else logmsg(10,"File of CCAAT box %s NOT read\n", ccaa_file);

        if((input=fopen(gcbox_file,"r"))) {
            logmsg(10, "Reading gcbox_file %s\n", gcbox_file);
            while(fgets(longstr,198,input) && !feof(input)) {
                sscanf(longstr,"%f %s %d %s",&tpr,ts,&gs,tatastr);
                if(gs>=start && gs<=end) {
                    gcbox_str= (char *)realloc(gcbox_str,(ngcb+1)*sizeof(char));
                    gcbox= (int *)realloc(gcbox,(ngcb+1)*sizeof(int));
                    if(ts[0]=='C') gcbox_str[ngcb]='C';
                    else gcbox_str[ngcb]='D';
                    gcbox[ngcb]= gs;
                    ++ngcb;
                }
            }
            fclose(input);
        }
        else logmsg(10,"File of GC box %s NOT read\n", gcbox_file);

        if((input=fopen(kozak_file,"r"))) {
            logmsg(10, "Reading kozak_file %s\n", kozak_file);
            while(fgets(longstr,198,input) && !feof(input)) {
                sscanf(longstr,"%f %s %d %s",&tpr,ts,&gs,tatastr);
                if(gs>=start && gs<=end) {
                    kozak_str= (char *)realloc(kozak_str,(nk+1)*sizeof(char));
                    kozak= (int *)realloc(kozak,(nk+1)*sizeof(int));
                    if(ts[0]=='C') kozak_str[nk]='C';
                    else kozak_str[nk]='D';
                    kozak[nk]= gs;
                    ++nk;
                }
            }
            fclose(input);
        }
        else logmsg(10,"File of Kozak sequences %s NOT read\n", kozak_file);


        if((input=fopen(pali_file,"r"))) {
            logmsg(10, "Reading pali_file %s\n", pali_file);
            while(fgets(longstr,198,input) && !feof(input)) {
                sscanf(longstr,"%f %d",&pp,&lp);
                if(pp>=start && pp<=end) {
                    pali= (int **)realloc(pali,(npali+1)*sizeof(int *));
                    pali[npali]= (int *)malloc(2*sizeof(int));
                    pali[npali][0]= (int)pp;
                    pali[npali][1]= lp;
                    ++npali;
                }
            }
            fclose(input);
        }
        else logmsg(10,"File of palindromes %s NOT read\n", pali_file);


        /* READS FILE OF UNBIASED ORFS */

        if((input=fopen(unb_file,"r"))) {
            logmsg(10, "Reading unb_file %s\n", unb_file);
            unbf= 1;
            while(fgets(longstr,198,input) && !feof(input)) {
                unb_str= (char *)realloc(unb_str,(nub+1)*sizeof(char));
                unb= (int **)realloc(unb,(nub+1)*sizeof(int *));
                unb[nub]= (int *)malloc(2*sizeof(int));
                p= strchr(longstr,'.'); p += 2;
			if(p[0] == '>' || p[0] == '<') ++p;
                ge= atoi(p);
                if(longstr[0]=='c')
		{
		p= longstr + 11;
			if(p[0] == '>' || p[0] == '<') ++p;
		gs= atoi(p);
		unb_str[nub]='C';
		}
                else
		{
		p= longstr;
			if(p[0] == '>' || p[0] == '<') ++p;
		gs= atoi(p);
		unb_str[nub]='D';
		}
                if(gs>=start-line_range/50 && gs<end && ge>start && ge<=end+line_range/50) { unb[nub][0]= gs; unb[nub][1]= ge; ++nub; }
                else if(gs>=start && gs<end && ge>end+line_range/50) { unb[nub][0]= gs; unb[nub][1]= end+line_range/50; ++nub; }
                else if(ge<=end+line_range/50 && ge>start && gs<start-line_range/50)
                { unb[nub][0]= start-line_range/50; unb[nub][0] += gs%period-unb[nub][0]%period; unb[nub][1]= ge; ++nub; }
                else if(gs<start-line_range/50 && ge>end+line_range/50)
                { unb[nub][0]= start-line_range/50; unb[nub][0] += gs%period-unb[nub][0]%period; unb[nub][1]= end+line_range/50; ++nub; }
            }
            fclose(input);
        }
        else logmsg(10,"Acc file NOT read\n") ;


        /* READS FILE OF CONSERVED ORFS */

        if((input= fopen(con_file,"r"))) {
            logmsg(10, "Reading con_file %s\n", con_file);
            conf= 1;
            while(fgets(longstr,198,input) && !feof(input)) {
                con_str= (char *)realloc(con_str,(nc+1)*sizeof(char));
                con= (int **)realloc(con,(nc+1)*sizeof(int *));
                con[nc]= (int *)malloc(2*sizeof(int));
                p= strchr(longstr,'.'); p += 2;
			if(p[0] == '>' || p[0] == '<') ++p;
                ge= atoi(p);
                if(longstr[0]=='c')
		{
		p= longstr + 11;
			if(p[0] == '>' || p[0] == '<') ++p;
		gs= atoi(p);
		con_str[nc]='C';
		}
                else
		{
		p= longstr;
			if(p[0] == '>' || p[0] == '<') ++p;
		gs= atoi(p);
		con_str[nc]='D';
		}
                if(gs>=start-line_range/50 && gs<end && ge<=end+line_range/50) { con[nc][0]= gs; con[nc][1]= ge; ++nc; }
                else if(gs>=start-line_range/50 && gs<end && ge>end+line_range/50) { con[nc][0]= gs; con[nc][1]= end+line_range/50; ++nc; }
                else if(ge<=end+line_range/50 && ge>start && gs<start-line_range/50)
                { con[nc][0]= start-line_range/50; con[nc][0] += gs%period-con[nc][0]%period; con[nc][1]= ge; ++nc; }
                else if(gs<start-line_range/50 && ge>end+line_range/50)
                { con[nc][0]= start-line_range/50; con[nc][0] += gs%period-con[nc][0]%period; con[nc][1]= end+line_range/50; ++nc; }
            }
            fclose(input);
        }
        else logmsg(10,"Conserved file NOT read\n") ;

        /* READS FILE OF NEW PROPOSED CODING REGIONS */

        if((input= fopen(new_file,"r"))) {
            logmsg(10, "Reading new_file %s\n", new_file);
            newf= 1;
            expand = 1;
            while(fgets(longstr,198,input) && !feof(input)) {
                if(expand) {
                    new_name= (char **)realloc(new_name, (nn + 1) * sizeof(char *));
                    new_name[nn]= (char *)malloc(50 * sizeof(char *));
                    new_str= (char *)realloc(new_str, (nn + 1) * sizeof(char));
                    new= (int **)realloc(new, (nn + 1)*sizeof(int *));
                    new[nn]= (int *)malloc(2 * sizeof(int));
                }
                strncpy(new_name[nn], longstr, 48);
                p= strchr(new_name[nn], ' ');
                p[0]= '\0';
                p= strrchr(longstr, ' '); ++p;
                p= strchr(p, '.'); p += 2;
                if(p[0] == '>' || p[0] == '<') ++p;
                ge= atoi(p);
                p= strrchr(longstr, ' '); ++p;
                if(p[0] == 'c') {
                    p += 11;
                    if(p[0] == '>' || p[0] == '<') ++p;
                    gs= atoi(p);
                    new_str[nn]='C';
                }
                else {
                    if(p[0] == '>' || p[0] == '<') ++p;
                    gs= atoi(p);
                    new_str[nn]='D';
                }
                expand = 1;   // Should we expand again next lop?
                if(gs >= start - line_range / 50 && gs < end && ge <= end + line_range / 50) { new[nn][0]= gs; new[nn][1]= ge; ++nn; }
                else if(gs >= start - line_range / 50 && gs < end && ge > end + line_range / 50) { new[nn][0]= gs; new[nn][1]= end + line_range / 50; ++nn; }
                else if(ge <= end + line_range / 50 && ge > start && gs < start - line_range / 50)
                { new[nn][0]= start - line_range / 50; new[nn][0] += gs % period - new[nn][0] % period; new[nn][1]= ge; ++nn; }
                else if(gs < start - line_range / 50 && ge > end + line_range / 50)
                { new[nn][0]= start - line_range / 50; new[nn][0] += gs % period - new[nn][0] % period; new[nn][1]= end + line_range / 50; ++nn; }
                else {
                    expand = 0; // didn't use the current slot, don't need to expand
                }
            }
            fclose(input);
        }
        else logmsg(10,"New file NOT read\n") ;

// READS FILE OF NEW PREDICTION CORRESPONDING TO ANNOTATED BUT WITH DIFFERENT PREDICTED START

        if((input= fopen(newP_file, "r"))) {
            logmsg(10, "Reading newP file %s\n", newP_file);
            newPf= 1;
            /* Header line not printed in current version of acgt_gamma */
//          fgets(longstr, 198, input);
            expand = 1;
            while(fgets(longstr,198,input) && !feof(input)) {
                if(expand) {
                    newP_name= (char **)realloc(newP_name, (nnP + 1) * sizeof(char *));
                    newP_name[nnP]= (char *)malloc(50 * sizeof(char *));
                    newP_str= (char *)realloc(newP_str, (nnP + 1) * sizeof(char));
                    newP= (int **)realloc(newP, (nnP + 1) * sizeof(int *));
                    newP[nnP]= (int *)malloc(2 * sizeof(int));
                }
                strncpy(newP_name[nnP], longstr, 48);
                p= strchr(newP_name[nnP], ' ');
                p[0]= '\0';
                p= strrchr(longstr, ' '); ++p;
                p= strchr(p, '.');
                ge= atoi(p + 2);
                p= strrchr(longstr, ' '); ++p;
                if(p[0] == 'c') {
                    p += 11;
                    gs= atoi(p);
                    newP_str[nnP]= 'C';
                }
                else {
                    gs= atoi(p);
                    newP_str[nnP]= 'D';
                }
                expand = 1;   // Should we expand again next lop?

                if(gs >= start - line_range / 50 && gs < end && ge <= end + line_range / 50) { newP[nnP][0]= gs; newP[nnP][1]= ge; ++nnP; }
                else if(gs >= start - line_range / 50 && gs < end && ge > end + line_range / 50) { newP[nnP][0]= gs; newP[nnP][1]= end + line_range / 50; ++nnP; }
                else if(ge <= end + line_range / 50 && ge > start && gs < start - line_range / 50)
                { newP[nnP][0]= start - line_range / 50; newP[nnP][0] += gs % period - newP[nnP][0] % period; newP[nnP][1]= ge; ++nnP; }
                else if(gs < start - line_range / 50 && ge > end + line_range / 50)
                { newP[nnP][0]= start - line_range / 50; newP[nnP][0] += gs % period - newP[nnP][0] % period; newP[nnP][1]= end + line_range / 50; ++nnP; }
                else {
                    expand = 0; // didn't use the current slot, don't need to expand
                }
            }
            fclose(input);
        }
        else logmsg(10,"newP file NOT read\n") ;



// READS FILE OF BLOCKS OF CONTRASTING S_PATTERNS

        if((input= fopen(cg_file,"r"))) { cgf= 1;
            logmsg(10, "Reading cg_file %s\n", cg_file);
            while(fgets(longstr,198,input) && !feof(input)) {
                cg_str= (char *)realloc(cg_str,(ncg+1)*sizeof(char));
                cg= (int **)realloc(cg,(ncg+1)*sizeof(int *));
                cg[ncg]= (int *)malloc(2*sizeof(int));
                p= strchr(longstr,'.');
                ge= atoi(p+2)+50;
                if(longstr[0]=='c') { gs= atoi(longstr+11)-50; cg_str[ncg]='C'; }
                else { gs= atoi(longstr)-50; cg_str[ncg]='D'; }
                if(gs>=start && gs<end && ge<=end) { cg[ncg][0]= gs; cg[ncg][1]= ge; ++ncg; }
                else if(gs>=start && gs<end && ge>end) { cg[ncg][0]= gs; cg[ncg][1]= end; ++ncg; }
                else if(ge<=end && ge>start && gs<start) { cg[ncg][0]= start; cg[ncg][1]= ge; ++ncg; }
                else if(gs<start && ge>end)
                { cg[ncg][0]= start; cg[ncg][1]= end; ++ncg; }
            }
            fclose(input);
        }
        else logmsg(10,"Files with blocks of asymmetric CG content NOT read\n") ;

        /* READS FILE OF ACCEPTED PUBLIC GENES */

        if((input= fopen(pub_file,"r"))) {
            logmsg(10, "Reading pub_file %s\n", pub_file);
            expand = 1;
            while(fgets(longstr, 198, input) && !feof(input)) {
                if(expand) {
                    pub_name= (char **)realloc(pub_name, (np + 1) * sizeof(char *));
                    pub_name[np]= (char *)malloc(50 * sizeof(char *));
                    pub_str= (char *)realloc(pub_str, (np + 1) * sizeof(char));
                    pub= (int **)realloc(pub, (np + 1) * sizeof(int *));
                    pub[np]= (int *)malloc(2 * sizeof(int));
                }
                strncpy(pub_name[np], longstr, 48);
                p= strchr(pub_name[np],' ');
                p[0]= '\0';
                p= strrchr(longstr, ' '); ++p;
                p= strchr(p, '.'); p += 2;
                if(p[0] == '>' || p[0] == '<') ++p;
                ge= atoi(p);
                p= strrchr(longstr, ' '); ++p;
                if(p[0]=='c') {
                    p += 11;
                    if(p[0] == '>' || p[0] == '<') ++p;
                    gs= atoi(p);
                    pub_str[np]='C';
                }
                else {
                    if(p[0] == '>' || p[0] == '<') ++p;
                    gs= atoi(p);
                    pub_str[np]='D';
                }

                expand = 1;   // Should we expand again next lop?
                if(gs >= start - line_range / 50 && gs < end && ge <= end + line_range / 50 && ge > start) { 
                    pub[np][0] = gs;
                    pub[np][1] = ge;
                    ++np;
                }
                else if(gs >= start - line_range / 50 && gs < end && ge > end + line_range / 50) {
                    pub[np][0] = gs;
                    pub[np][1] = end + line_range / 50;
                    ++np;
                }
                else if(ge <= end + line_range / 50 && ge > start && gs < start - line_range / 50) {
                    pub[np][0]  = start - line_range / 50;
                    pub[np][0] += gs % period - pub[np][0] % period;
                    pub[np][1]  = ge;
                    ++np;
                }
                else if(gs < start - line_range / 50 && ge > end + line_range / 50) {
                    pub[np][0]  = start - line_range / 50;
                    pub[np][0] += gs % period - pub[np][0] % period;
                    pub[np][1]  = end + line_range / 50;
                    ++np;
                }
                else {
                    expand = 0; // didn't use the current slot, don't need to expand
                }
            }
            fclose(input);
        }
        else logmsg(10,"Pub file NOT read\n") ;


// READS FILE OF MODIFIED PUBLIC GENES.
// FILE OF NEW PREDICTIONS MODIFYING START OF TRANSLATION OF ANNOTATED GENES
// IS READ INSTEAD INTO NEW POTENTIAL CODING REGIONS (newP* variables)

        if((input= fopen(mod_file,"r"))) {
            logmsg(10, "Reading modified-predictions file %s\n", mod_file);
            /* Skip passed the header line */
            fgets(longstr, 198, input);
            while(fgets(longstr, 198, input) && !feof(input)) {
                mod_name= (char **)realloc(mod_name, (ne + 1) * sizeof(char *));
                mod_name[ne]= (char *)malloc(50 * sizeof(char *));
                mod_str= (char *)realloc(mod_str, (ne + 1)*sizeof(char));
                modified= (int **)realloc(modified,(ne + 1)*sizeof(int *));
                modified[ne]= (int *)malloc(2 * sizeof(int));
                strncpy(mod_name[ne], longstr, 48);
                p= strchr(mod_name[ne],' ');
                p[0]= '\0';
                p= strrchr(longstr,' '); ++p;
                p= strchr(p, '.'); p += 2;
			if(p[0] == '<' || p[0] == '>') ++p;
                ge= atoi(p);
                p= strrchr(longstr,' '); ++p;
                if(p[0] =='c')
		{
		p += 11;
			if(p[0] == '>' || p[0] == '<') ++p;
		gs= atoi(p);
		mod_str[ne]='C';
		}
                else
		{
			if(p[0] == '>' || p[0] == '<') ++p;
		gs= atoi(p);
		mod_str[ne]='D';
		}

                if(gs >= start - line_range / 50 && gs < end && ge <= end + line_range / 50 && ge > start) { 
                    modified[ne][0] = gs;
                    modified[ne][1] = ge;
                    ++ne;
                }
                else if(gs >= start - line_range / 50 && gs < end && ge > end + line_range / 50) {
                    modified[ne][0] = gs;
                    modified[ne][1] = end + line_range / 50;
                    ++ne;
                }
                else if(ge <= end + line_range / 50 && ge > start && gs < start - line_range / 50) {
                    modified[ne][0]  = start - line_range / 50;
                    modified[ne][0] += gs % period - modified[ne][0] % period;
                    modified[ne][1]  = ge;
                    ++ne;
                }
                else if(gs < start - line_range / 50 && ge > end + line_range / 50) {
                    modified[ne][0]  = start - line_range / 50;
                    modified[ne][0] += gs % period - modified[ne][0] % period;
                    modified[ne][1]  = end + line_range / 50;
                    ++ne;
                }
            }
            fclose(input);
            logmsg(10,"Modified file %s read\n", mod_file);
        }
        else logmsg(10,"Modified file NOT read\n") ;

        if((input= fopen(CG200_file,"r"))) {
            logmsg(10, "Reading CG200_file %s\n", CG200_file);
            while(!feof(input)) {
                fscanf(input,"%d",&pos);
                for(j=0; j<period; ++j) fscanf(input,"%f",S+j);
                if(pos>=start && pos<=end) {
                    x= (int *)realloc(x,(n+1)*sizeof(int));
                    y= (float **)realloc(y,(n+1)*sizeof(float *));
                    y[n]= (float *)malloc((period+1)*sizeof(float));
                    y[n][period]= 0.0;
                    x[n]= pos;
                    for(j=0;j<period;++j) { y[n][j]= S[j]; y[n][period] += y[n][j]; }
                    y[n][period] /= (float)period;
                    ++n;
                    pos= -1;
                }
                else if(pos > end)
                    break;
                else
                    fgets(longstr,198,input);
            }

            fclose(input);
            logmsg(10,"Large window composition file %s read\n", CG200_file);
        }
        else 
            logmsg(10,"Large-window composition file NOT read\n");

// READS FILE OF RNA-SEQ READ NUMBERS

/*
        if((input= fopen(read_file,"r"))) {
            logmsg(10, "Reading read_file %s\n", read_file);
        fscanf(input,"%d",&wind);
          while(!feof(input))
          {
          fscanf(input,"%d",&pos);
            for(j= 0; j < 4; ++j) fscanf(input,"%f", S + j);
            if(pos >= start && pos < end)
            {
            X= (int *)realloc(X, (N + 1) * sizeof(int));
            Y= (float **)realloc(Y, (N + 1) * sizeof(float *));
            Y[N]= (float *)malloc(4 * sizeof(float));
            X[N]= pos;
              for(j= 0; j < 4; ++j) Y[N][j]= S[j];
            ++N;
            pos= -1;
            }
            else fgets(longstr, 198, input);
          }
          fclose(input);
	  swflag= 1;
          logmsg(10,"Read-numbers file %s read\n",  read_file); }
          else { logmsg(10, "Read-numbers file NOT read\n") ; }


        if(!k) {
            fprintf(stdout,"%%!PS-Adobe-2.0\n\n");
            fprintf(stdout,"gsave\n");
            fprintf(stdout,"%d dict begin\n",DICT_SIZE);
            fprintf(stdout,"/L025 { 0.25 setlinewidth } def\n");
            fprintf(stdout,"/L05 { 0.5 setlinewidth } def\n");
            fprintf(stdout,"/L1 { 1.0 setlinewidth } def\n");
            fprintf(stdout,"/L15 { 1.5 setlinewidth } def\n");
            fprintf(stdout,"/L2 { 2.0 setlinewidth } def\n");
            fprintf(stdout,"/L3 { 3.0 setlinewidth } def\n");
            fprintf(stdout,"/L4 { 4.0 setlinewidth } def\n");
            fprintf(stdout,"/M {moveto} def\n");
            fprintf(stdout,"/L {lineto} def\n");
            fprintf(stdout,"/RM {rmoveto} def\n");
            fprintf(stdout,"/RL {rlineto} def\n");
            fprintf(stdout,"/NumberFontSize {%1.f} def\n",NFS);
            fprintf(stdout,"/NamesFontSize {%1.f} def\n",NaFS);
            fprintf(stdout,"/TitleFontSize {%.1f} def\n",TFS);
            fprintf(stdout,"/Title2FontSize {%.1f} def\n",TFS);
            fprintf(stdout,"/LegendFontSize {%.1f} def\n",LFS);

            if(ANOPIAS)
	    {
            fprintf(stdout,"/R {0.835 0.369 0.000 setrgbcolor} def\n");  // Vermillion
            fprintf(stdout,"/G {0.800 0.475 0.655 setrgbcolor} def\n");  // Reddish purple
            fprintf(stdout,"/B {0.000 0.447 0.698 setrgbcolor} def\n");  // Blue
            fprintf(stdout,"/LR {1.000 0.592 0.027 setrgbcolor} def\n");  // Light Vermillion
            fprintf(stdout,"/LG {0.914 0.686 0.812 setrgbcolor} def\n");  // Light Reddish purple
            fprintf(stdout,"/LB {0.300 0.724 0.849 setrgbcolor} def\n");  // Light Blue
	    }
            else
	    {
            fprintf(stdout,"/R {1.0 0.0 0.0 setrgbcolor} def\n");
            fprintf(stdout,"/G {0.0 1.0 0.0 setrgbcolor} def\n");
            fprintf(stdout,"/B {0.0 0.0 1.0 setrgbcolor} def\n");
            fprintf(stdout,"/LR {1.0 0.5 0.5 setrgbcolor} def\n");
            fprintf(stdout,"/LG {0.5 1.0 0.5 setrgbcolor} def\n");
            fprintf(stdout,"/LB {0.5 0.5 1.0 setrgbcolor} def\n");
	    }

            fprintf(stdout,"/DarkGray {0.25 0.25 0.25 setrgbcolor} def\n");
            fprintf(stdout,"/Gray {0.5 0.5 0.5 setrgbcolor} def\n");
            fprintf(stdout,"/LightGray {0.75 0.75 0.75 setrgbcolor} def\n");
            fprintf(stdout,"/GrayShade {0.95 0.95 0.95 setrgbcolor} def\n");
            fprintf(stdout,"/Black {0.0 0.0 0.0 setrgbcolor} def\n");
            fprintf(stdout,"/Cshow { dup stringwidth pop -2 div 0 RM show } def\n");
            fprintf(stdout,"/Lshow { dup stringwidth pop -1 mul 0 RM show } def\n");
//          fprintf(stdout,"/Rarrow { 1 setlinejoin 1 setlinecap dup 2 lt { 0 2 RM dup -2 RL -1 mul -2 RL } { 0 2 RM 2 sub dup 0 RL 2 -2 RL -2 -2 RL -1 mul 0 RL } ifelse closepath stroke 0 setlinejoin 0 setlinecap } def\n");
//          fprintf(stdout,"/Larrow { 1 setlinejoin 1 setlinecap dup 2 lt { 2 RL 0 -4 RL } { 2 2 RL dup 2 sub dup 0 RL 0 -4 RL -1 mul 0 RL } ifelse closepath stroke 0 setlinejoin 0 setlinecap } def\n");
           fprintf(stdout,"/Rarrow { 1 setlinejoin 1 setlinecap dup 3 lt { 0 3 RM dup -3 RL -1 mul -3 RL } { 0 3 RM 3 sub dup 0 RL 3 -3 RL -3 -3 RL -1 mul 0 RL } ifelse closepath stroke 0 setlinejoin 0 setlinecap } def\n");
           fprintf(stdout,"/Larrow { 1 setlinejoin 1 setlinecap dup 3 lt { 3 RL 0 -6 RL } { 3 3 RL dup 3 sub dup 0 RL 0 -6 RL -1 mul 0 RL } ifelse closepath stroke 0 setlinejoin 0 setlinecap } def\n");
            fprintf(stdout,"/radius { %.1f } def\n",HAIRPIN_RADIUS);
            fprintf(stdout,"/thick { %.2f } def\n",STEM_SEPARATION);
            fprintf(stdout,"/hairpin { dup thick -2.0 div exch RL thick 0 RL -1 mul thick -2.0 div exch RL closepath } def\n");
            fprintf(stdout,"/Ucol { dup 0.0 exch RL %.3f 0.0 RL %.3f exch dup 0.0 exch -1 mul RL pop pop } def\n", (float)wind / delta * WIDTH, (float)wind / delta * WIDTH);
            fprintf(stdout,"/Dcol { dup 0.0 exch -1 mul RL %.3f 0.0 RL %.3f exch dup 0.0 exch RL pop pop } def\n", (float)wind / delta * WIDTH, (float)wind / delta * WIDTH);
	    fprintf(stdout,"/ShadeBox { dup -1 mul 3 1 roll 0 exch RL 0 RL 0 exch RL closepath fill } def\n");
            fprintf(stdout,"100 %d translate\n\n",TITLES_POSITION);


            if(tstart==0) {
                fprintf(stdout,"/Times-Bold findfont TitleFontSize scalefont setfont\n");
                fprintf(stdout,"%.3f %.3f M ",0.5*WIDTH-30,HIGHT-90.0);
                fprintf(stdout,"(%s-profiles of %s) Cshow\n", NUCLEOTIDES, title1);
            }
            else {
                fprintf(stdout,"/Times-Italic findfont Title2FontSize scalefont setfont\n");
                fprintf(stdout,"%.3f %.3f M ",0.0,HIGHT-90.0);
                fprintf(stdout,"(%s-profiles - %s) show\n", NUCLEOTIDES, title2);
            }
        }

        if(k) fprintf(stdout,"0 %d translate\n",TRANSLATE);
        else  fprintf(stdout,"0 %d translate\n",FIRST_TRANS);

        /*******************************************************/
        /* Prints boxes of unexplained regions of CG asymmetry */
        /*******************************************************/

        if(cgf) {
            fprintf(stdout,"L025\n");
            for(i=0;i<ncg;++i) {
                fprintf(stdout,"DarkGray\n");
                fprintf(stdout," %.1f %.2f M %.1f %.2f L %.1f %.1f RL %.2f %.2f L closepath fill\n",(cg[i][0]-(float)start)/delta*WIDTH,HIGHT,(cg[i][1]-(float)start)/delta*WIDTH,HIGHT,0.0,-HIGHT,(cg[i][0]-(float)start)/delta*WIDTH,0.0);
            }
        }

        /************************************************/
        /* End printing regions of nucleotide contrasts */
        /************************************************/

        fprintf(stdout,"/Helvetica-Roman findfont NumberFontSize scalefont setfont\n");

/*  Vertical grid: */

//      fprintf(stdout,"L025 Gray\n");
//      for(r= TIC_X / 5.0; r < end - start; r += TIC_X / 5.0)
//          fprintf(stdout,"%.3f 0.0 M 0 %.1f RL\n", r / delta * WIDTH, HIGHT);
//      fprintf(stdout,"stroke\n");

/*  Shade pattern: */

        fprintf(stdout,"GrayShade\n");
        for(r= TIC_X / 2.0; r < end - start - TIC_X / 2.0; r += 2.0 * TIC_X / 2.0)
            fprintf(stdout,"%.3f 0.0 M %.1f %.1f ShadeBox\n", r / delta * WIDTH, (TIC_X * WIDTH) / (2.0 * delta), HIGHT);
	if(r < end - start) fprintf(stdout,"%.3f 0.0 M %.1f %.1f ShadeBox\n", r / delta * WIDTH, (end - start - r) * WIDTH / delta, HIGHT);

// X-axis tics and labels

        fprintf(stdout, "L025 Black\n");
        for(r= 0; r <= end - start; r += TIC_X)
	{
//          fprintf(stdout,"%.3f %.3f M 0 %.3f RL 0 %.3f RM (%.0f) Cshow\n", r / delta * WIDTH, -TIC_W, TIC_W, -8.0 - TIC_W, r + start);
            fprintf(stdout,"%.3f %.3f M 0 %.3f RL 0 %.3f RM (", r / delta * WIDTH, -TIC_W, TIC_W, -8.0 - TIC_W);
	    sprintf(pos_string, "%.0f", r + start);

            for(i= 0; i < strlen(pos_string) - 1; ++i)
            {
            fprintf(stdout,"%c", pos_string[i]);
                if(!((strlen(pos_string) - 1 - i) % 3)) fprintf(stdout, ",");
            }
            fprintf(stdout,"%c) Cshow\n", pos_string[i]);
	}
	fprintf(stdout, "stroke\n");

  if(swflag)
  {
  TIC_Y= MAX / (MAXLOG * 2.0);

    for(i= 0; i <= MAXLOG; ++i)
    {
    fprintf(stdout, "%.3f %.3f M %.3f 0 RL %.3f %.3f RM ", -TIC_W, HIGHT / 2.0 + i * TIC_Y / MAX * HIGHT, TIC_W, -TIC_W - 5.0, -3.0);
      if(!(i % EVERY)) fprintf(stdout, "(%d) Lshow\n", i);
      else             fprintf(stdout, "\n");
    }
  fprintf(stdout, "stroke\n");
  }
  else
  {
  TIC_Y= 20.0;

  for(r= 0; r <= MAX; r += TIC_Y)
    fprintf(stdout, "%.3f %.3f M %.3f 0 RL %.3f %.3f RM (%.1f) Lshow\n", -TIC_W, r / MAX * HIGHT, TIC_W, -TIC_W - 3.0, -3.0, r);
  fprintf(stdout, "stroke\n");

/*  Horizontal grid: */

//  fprintf(stdout,"L025 Gray\n");

//  for(r= TIC_Y; r < MAX; r += TIC_Y)
//    fprintf(stdout, "0.0 %.3f M %.3f 0 RL\n", r / MAX * HIGHT, (end - start) / delta * WIDTH);
//  fprintf(stdout, "stroke\n");
  }

fprintf(stdout,"L05 Black\n");
fprintf(stdout,"newpath 0 0 M 0 %.3f RL %.3f 0 RL 0 %.3f RL %.3f 0 RL closepath stroke\n",HIGHT,(end-start)/delta*WIDTH,-HIGHT,(start-end)/delta*WIDTH);

fprintf(stdout,"Black %.3f %.3f M (Input file CDS) Lshow\n",-15.0,HIGHT+HIGHT_PUB - 2);

  if(swflag)
  {
  fprintf(stdout,"%.3f %.3f M (0.0) show\n", (end-start)/delta*WIDTH + 3.0, -2.0);
  fprintf(stdout,"%.3f %.3f M (50.0) show\n", (end-start)/delta*WIDTH + 3.0, HIGHT / 2.0 - 2.0);
  fprintf(stdout,"%.3f %.3f M (100.0) show\n", (end-start)/delta*WIDTH + 3.0, HIGHT - 2.0);
  fprintf(stdout,"90 rotate %.3f %.3f M (%% %s) Cshow -90 rotate\n", 0.5 * HIGHT, -((end-start)/delta*WIDTH) - 25.0, NUCLEOTIDES);
  fprintf(stdout," 0.0 0.0 M 90 rotate %.3f %.3f M (Log-count) Cshow -90 rotate\n", 0.5 * HIGHT, 20.0);
  }

        // Prints HSSs

        if(Scpf) {
            fprintf(stdout,"%.3f %.3f M (Hits) Lshow\n",-15.0, HIGHT + HIGHT_SCP - 2);
            fprintf(stdout,"L025 LightGray\n");
            fprintf(stdout,"-4 %.3f M %.3f 0 RL -3 +3 RL stroke\n",HIGHT + HIGHT_SCP + 2, (end - start) / delta * WIDTH + 8);
            fprintf(stdout,"-1 %.3f -3 add M -3 3 RL %.3f 0 RL stroke\n", HIGHT + HIGHT_SCP - 2,(end - start) / delta * WIDTH + 8);
            fprintf(stdout," Black\n");
        }

        if(conf) {
            fprintf(stdout,"%.3f %.3f M (Conserved) Lshow\n",-25.0,HIGHT+HIGHT_HOM-2);
            fprintf(stdout," L025 LightGray\n");
            fprintf(stdout,"0 %.3f M %.3f 0 RL stroke\n", HIGHT+HIGHT_HOM + 1.5, (end - start) / delta * WIDTH);
            fprintf(stdout,"0 %.3f M %.3f 0 RL stroke\n", HIGHT+HIGHT_HOM - 1.5, (end - start) / delta * WIDTH);
            fprintf(stdout," Black\n");
        }
        if(newf) fprintf(stdout,"%.3f %.3f M (Newly identified ORFs) Lshow\n", -15.0, HIGHT + HIGHT_NEW-2);

        fprintf(stdout,"/Helvetica-Roman findfont LegendFontSize scalefont setfont Black\n");
        if(!swflag) fprintf(stdout,"90 rotate %.3f %.3f M (%% %s) Cshow -90 rotate\n", 0.5 * HIGHT, 25.0, NUCLEOTIDES);
        if(k == (lines - 1) || (end-start) < line_range) fprintf(stdout,"%.3f %.3f M (Sequence position / nt) Cshow\n", 0.5 * WIDTH, -30.0);

        fprintf(stdout,"stroke DarkGray\n");

        fprintf(stdout,"%.1f %.2f M\n",(x[0]-(float)start)/delta*WIDTH,HIGHT*y[0][period]/MAX);
        for(i=1;i<n;++i)
            fprintf(stdout,"%.1f %.2f L\n",(x[i]-(float)start)/delta*WIDTH,HIGHT*y[i][period]/MAX);

      /**************************************/
      /* Prints histograms of RNA-seq reads */
      /**************************************/
      
        if(swflag)
        {
            logmsg(10, "\nPrints profiles of RNA-seq reads\n");
        fprintf(stdout,"Black\n");
      
      // Writes histogram of most frequent reads
      
        fprintf(stdout, "newpath %.1f %.2f M\n", (X[0] - (float)start) / delta * WIDTH, 0.5 * HIGHT);
          for(i= 0; i < N; ++i)
          {
              if(Y[i][0] + Y[i][2] > Y[i][1] + Y[i][3]) j= 0;
              else j= 1;
      
              if(Y[i][j] == 0.0) r= 0.0;
              else
              {
              r= (float)log10((double)Y[i][j]) / MAXLOG * (HIGHT / 2.0);
                if(r < 0.25) r= 0.25;
              }
            if(j) fprintf(stdout, "%.3f Dcol\n", r);
            else  fprintf(stdout, "%.3f Ucol\n", r);
          }
      
      // Writes histogram of less frequent reads
      
        fprintf(stdout, "closepath fill newpath Gray %.1f %.2f M\n", (X[0] - (float)start) / delta * WIDTH, 0.5 * HIGHT);
          for(i= 0; i < N; ++i)
          {
              if(Y[i][0] + Y[i][2] <= Y[i][1] + Y[i][3]) j= 0;
              else j= 1;
      
              if(Y[i][j] == 0.0) r= 0.0;
              else
              {
              r= (float)log10((double)Y[i][j]) / MAXLOG * (HIGHT / 2.0);
                if(r < 0.25) r= 0.25;
              }
            if(j) fprintf(stdout, "%.3f Dcol\n", r);
            else  fprintf(stdout, "%.3f Ucol\n", r);
          }
      
        fprintf(stdout,"closepath fill\n");
      
        fprintf(stdout,"L025 Black 0 setlinejoin 0 setlinecap\n");
      
      // Adds histogram of most frequent multiple reads
      
          for(i= 0; i < N; ++i)
          {
            if(Y[i][0] + Y[i][2] > Y[i][1] + Y[i][3]) j= 0;
            else j= 1;
      
            if(Y[i][j + 2])
            {
              if(Y[i][j])
              {
              y1= (float)log10((double)Y[i][j]) / MAXLOG * (HIGHT / 2.0);
                if(y1 < 0.25) y1= 0.25;
              }
              else y1= 0.0;
      
              if(Y[i][j] + Y[i][j + 2])
              {
              y2= (float)log10((double)(Y[i][j] + Y[i][j + 2])) / MAXLOG * (HIGHT / 2.0);
                if(y2 < 0.25) y2= 0.25;
              }
              else y2= 0.0;
      
              if(y2 > y1)
              {
                if(j) fprintf(stdout, "newpath %.3f %.3f %.3f Dbox\n", y2 - y1 - 0.125, (X[i] - (float)start) / delta * WIDTH, 0.5 * HIGHT - y1);
                else  fprintf(stdout, "newpath %.3f %.3f %.3f Ubox\n", y2 - y1 - 0.125, (X[i] - (float)start) / delta * WIDTH, 0.5 * HIGHT + y1);
              }
      
            fprintf(stdout, "stroke\n");
            }
          }
      
        fprintf(stdout,"Gray\n");
      
      // Adds histogram of least frequent multiple reads
      
          for(i= 0; i < N; ++i)
          {
            if(Y[i][0] + Y[i][2] <= Y[i][1] + Y[i][3]) j= 0;
            else j= 1;
      
            if(Y[i][j + 2])
            {
              if(Y[i][j])
              {
              y1= (float)log10((double)Y[i][j]) / MAXLOG * (HIGHT / 2.0);
                if(y1 < 0.25) y1= 0.25;
              }
              else y1= 0.0;
      
              if(Y[i][j] + Y[i][j + 2])
              {
              y2= (float)log10((double)(Y[i][j] + Y[i][j + 2])) / MAXLOG * (HIGHT / 2.0);
                if(y2 < 0.25) y2= 0.25;
              }
              else y2= 0.0;
      
              if(y2 > y1)
              {
                if(j) fprintf(stdout, "newpath %.3f %.3f %.3f Dbox\n", y2 - y1 - 0.125, (X[i] - (float)start) / delta * WIDTH, 0.5 * HIGHT - y1);
                else  fprintf(stdout, "newpath %.3f %.3f %.3f Ubox\n", y2 - y1 - 0.125, (X[i] - (float)start) / delta * WIDTH, 0.5 * HIGHT + y1);
              }
      
            fprintf(stdout, "stroke\n");
            }
          }
      
        }

        /****************************************/
        /* Prints S-profiles from large windows */
        /****************************************/

        fprintf(stdout,"\nstroke L05 R %.1f %.2f M 1 setlinejoin 1 setlinecap\n",(x[0]-(float)start)/delta*WIDTH,y[0][0]/MAX*HIGHT);
        for(i=1;i<n;++i)
            fprintf(stdout,"%.1f %.2f L\n",(x[i]-(float)start)/delta*WIDTH,y[i][0]/MAX*HIGHT);

        if(period>=2) {
            fprintf(stdout,"stroke L05 G\n");

            fprintf(stdout,"%.1f %.2f M\n",(x[0]-(float)start)/delta*WIDTH,y[0][1]/MAX*HIGHT);
            for(i=1;i<n;++i)
                fprintf(stdout,"%.1f %.2f L\n",(x[i]-(float)start)/delta*WIDTH,y[i][1]/MAX*HIGHT);

            if(period>=3) {
                fprintf(stdout,"stroke L05 B\n");

                fprintf(stdout,"%.1f %.2f M\n",(x[0]-(float)start)/delta*WIDTH,y[0][2]/MAX*HIGHT);
                for(i=1;i<n;++i)
                    fprintf(stdout,"%.1f %.2f L\n",(x[i]-(float)start)/delta*WIDTH,y[i][2]/MAX*HIGHT);
            }
        }

        fprintf(stdout,"stroke 0 setlinejoin 0 setlinecap\n");

        /****************************************/
        /* Prints new ORFs with low biases (DISABLED)

           if(unbf)
           {
           for(i=0;i<nub;++i)
           if(unb_str[i]=='D')
           {
           if(unb[i][0]%period==1) fprintf(stdout,"B");
           else if(unb[i][0]%period==2) fprintf(stdout,"R");
           else if(unb[i][0]%period==0) fprintf(stdout,"G");
           fprintf(stdout," %.1f %.2f M 0 2 RM 0 -4 RL 0 2 RM %.1f %.2f L -2 2 RM 2 -2 RL -2 -2 RL stroke\n",(unb[i][0]-(float)start)/delta*WIDTH,HIGHT+HIGHT_BIA+2.0,(unb[i][1]-(float)start)/delta*WIDTH,HIGHT+HIGHT_BIA+2.0);
           }
           else
           {
           if(unb[i][0]%period==1) fprintf(stdout,"R");
           else if(unb[i][0]%period==2) fprintf(stdout,"G");
           else if(unb[i][0]%period==0) fprintf(stdout,"B");
           fprintf(stdout," %.1f %.2f M 2 2 RL 0 -4 RM -2 2 RL %.1f %.2f L 0 2 RM 0 -4 RL stroke\n",(unb[i][0]-(float)start)/delta*WIDTH,HIGHT+HIGHT_BIA-2,(unb[i][1]-(float)start)/delta*WIDTH,HIGHT+HIGHT_BIA-2.0);
           }
           }

        */

        /*********************************/
        /* Prints new ORFs with homology */
        /*********************************/

        if(conf) {
            for(i=0;i<nc;++i) {
                if(con_str[i]=='D') {
                    if(con[i][0]%period==1) fprintf(stdout,"B");
                    else if(con[i][0]%period==2) fprintf(stdout,"R");
                    else if(con[i][0]%period==0) fprintf(stdout,"G");
                    fprintf(stdout," %.1f %.2f M 0 2 RM 0 -4 RL 0 2 RM %.1f %.2f L -2 2 RM 2 -2 RL -2 -2 RL stroke\n",(con[i][0]-(float)start)/delta*WIDTH,HIGHT+HIGHT_HOM+1.5,(con[i][1]-(float)start)/delta*WIDTH,HIGHT+HIGHT_HOM+1.5);
                }
                else {
                    if(con[i][0]%period==1) fprintf(stdout,"R");
                    else if(con[i][0]%period==2) fprintf(stdout,"G");
                    else if(con[i][0]%period==0) fprintf(stdout,"B");
                    fprintf(stdout," %.1f %.2f M 2 2 RL 0 -4 RM -2 2 RL %.1f %.2f L 0 2 RM 0 -4 RL stroke\n",(con[i][0]-(float)start)/delta*WIDTH,HIGHT+HIGHT_HOM-1.5,(con[i][1]-(float)start)/delta*WIDTH,HIGHT+HIGHT_HOM-1.5);
                }
            }
        }

        /****************************************/
        /* Prints S-profiles from small windows */
        /****************************************/

/*
        if(swflag) {
            fprintf(stdout,"stroke 1 setlinejoin 1 setlinecap L05 LR\n");

            fprintf(stdout,"%.1f %.2f M\n",(X[0]-(float)start)/delta*WIDTH,Y[0][0]/MAX*HIGHT);
            for(i=1;i<N;++i)
                fprintf(stdout,"%.1f %.2f L\n",(X[i]-(float)start)/delta*WIDTH,Y[i][0]/MAX*HIGHT);

            if(period>=2) {
                fprintf(stdout,"stroke L05 LG\n");

                fprintf(stdout,"%.1f %.2f M\n",(X[0]-(float)start)/delta*WIDTH,Y[0][1]/MAX*HIGHT);
                for(i=1;i<N;++i)
                    fprintf(stdout,"%.1f %.2f L\n",(X[i]-(float)start)/delta*WIDTH,Y[i][1]/MAX*HIGHT);

                if(period>=3) {
                    fprintf(stdout,"stroke L05 LB\n");

                    fprintf(stdout,"%.1f %.2f M\n",(X[0]-(float)start)/delta*WIDTH,Y[0][2]/MAX*HIGHT);
                    for(i=1;i<N;++i)
                        fprintf(stdout,"%.1f %.2f L\n",(X[i]-(float)start)/delta*WIDTH,Y[i][2]/MAX*HIGHT);
                }
            }
        }

        fprintf(stdout,"stroke 1 setlinejoin 1 setlinecap L05 R\n");
*/

        /****************************************/
        /* Prints S-profiles from large windows */
        /****************************************/

/*
        fprintf(stdout,"%.1f %.2f M\n",(x[0]-(float)start)/delta*WIDTH,y[0][0]/MAX*HIGHT);
        for(i=1;i<n;++i)
            fprintf(stdout,"%.1f %.2f L\n",(x[i]-(float)start)/delta*WIDTH,y[i][0]/MAX*HIGHT);

        if(period>=2) {
            fprintf(stdout,"stroke L05 G\n");

            fprintf(stdout,"%.1f %.2f M\n",(x[0]-(float)start)/delta*WIDTH,y[0][1]/MAX*HIGHT);
            for(i=1;i<n;++i)
                fprintf(stdout,"%.1f %.2f L\n",(x[i]-(float)start)/delta*WIDTH,y[i][1]/MAX*HIGHT);

            if(period>=3) {
                fprintf(stdout,"stroke L05 B\n");

                fprintf(stdout,"%.1f %.2f M\n",(x[0]-(float)start)/delta*WIDTH,y[0][2]/MAX*HIGHT);
                for(i=1;i<n;++i)
                    fprintf(stdout,"%.1f %.2f L\n",(x[i]-(float)start)/delta*WIDTH,y[i][2]/MAX*HIGHT);
            }
        }

        fprintf(stdout,"stroke\n");
*/

        /****************************************/
        /* Prints new ORFs with low biases (DISABLED)

fprintf(stdout,"\n0 setlinejoin 0 setlinecap\n");

           if(unbf)
           {
           for(i=0;i<nub;++i)
           if(unb_str[i]=='D')
           {
           if(unb[i][0]%period==1) fprintf(stdout,"B");
           else if(unb[i][0]%period==2) fprintf(stdout,"R");
           else if(unb[i][0]%period==0) fprintf(stdout,"G");
           fprintf(stdout," %.1f %.2f M 0 2 RM 0 -4 RL 0 2 RM %.1f %.2f L -2 2 RM 2 -2 RL -2 -2 RL stroke\n",(unb[i][0]-(float)start)/delta*WIDTH,HIGHT+HIGHT_BIA+2.0,(unb[i][1]-(float)start)/delta*WIDTH,HIGHT+HIGHT_BIA+2.0);
           }
           else
           {
           if(unb[i][0]%period==1) fprintf(stdout,"R");
           else if(unb[i][0]%period==2) fprintf(stdout,"G");
           else if(unb[i][0]%period==0) fprintf(stdout,"B");
           fprintf(stdout," %.1f %.2f M 2 2 RL 0 -4 RM -2 2 RL %.1f %.2f L 0 2 RM 0 -4 RL stroke\n",(unb[i][0]-(float)start)/delta*WIDTH,HIGHT+HIGHT_BIA-2,(unb[i][1]-(float)start)/delta*WIDTH,HIGHT+HIGHT_BIA-2.0);
           }
           }

        */

        /*********************************/
        /* Prints new ORFs with homology */
        /*********************************/

        if(conf) {
            for(i=0;i<nc;++i) {
                if(con_str[i]=='D') {
                    if(con[i][0]%period==1) fprintf(stdout,"B");
                    else if(con[i][0]%period==2) fprintf(stdout,"R");
                    else if(con[i][0]%period==0) fprintf(stdout,"G");
                    fprintf(stdout," %.1f %.2f M 0 2 RM 0 -4 RL 0 2 RM %.1f %.2f L -2 2 RM 2 -2 RL -2 -2 RL stroke\n",(con[i][0]-(float)start)/delta*WIDTH,HIGHT+HIGHT_HOM+1.5,(con[i][1]-(float)start)/delta*WIDTH,HIGHT+HIGHT_HOM+1.5);
                }
                else {
                    if(con[i][0]%period==1) fprintf(stdout,"R");
                    else if(con[i][0]%period==2) fprintf(stdout,"G");
                    else if(con[i][0]%period==0) fprintf(stdout,"B");
                    fprintf(stdout," %.1f %.2f M 2 2 RL 0 -4 RM -2 2 RL %.1f %.2f L 0 2 RM 0 -4 RL stroke\n",(con[i][0]-(float)start)/delta*WIDTH,HIGHT+HIGHT_HOM-1.5,(con[i][1]-(float)start)/delta*WIDTH,HIGHT+HIGHT_HOM-1.5);
                }
            }
        }

// Prints stops

        /*

          for(i=0;i<ns;++i)
          if(stop_str[i]=='D')
          {
          if(stop[i]%period==1) fprintf(stdout,"B");
          else if(stop[i]%period==2) fprintf(stdout,"R");
          else if(stop[i]%period==0) fprintf(stdout,"G");
          fprintf(stdout," %.1f %.2f M 0 2 RM 0 -4 RL stroke\n",(stop[i]-(float)start)/delta*WIDTH,HIGHT+HIGHT_MET + 2.0);
          }
          else
          {
          if(stop[i]%period==1) fprintf(stdout,"G");
          else if(stop[i]%period==2) fprintf(stdout,"B");
          else if(stop[i]%period==0) fprintf(stdout,"R");
          fprintf(stdout," %.1f %.2f M 0 2 RM 0 -4 RL stroke\n",(stop[i]-(float)start)/delta*WIDTH,HIGHT+HIGHT_MET - 2.0);
          }

        */

// Prints methionines

        /*
          fprintf(stdout,"\nL2 1 setlinecap\n");

          for(i=0;i<nm;++i)
          if(met_str[i]=='D')
          {
          if(met[i]%period==1) fprintf(stdout,"B");
          else if(met[i]%period==2) fprintf(stdout,"R");
          else if(met[i]%period==0) fprintf(stdout,"G");
          fprintf(stdout," %.1f %.2f M 0 0 RL stroke\n",(met[i]-(float)start)/delta*WIDTH,HIGHT+HIGHT_MET + 2.0);
          }
          else
          {
          if(met[i]%period==1) fprintf(stdout,"G");
          else if(met[i]%period==2) fprintf(stdout,"B");
          else if(met[i]%period==0) fprintf(stdout,"R");
          fprintf(stdout," %.1f %.2f M 0 0 RL stroke\n",(met[i]-(float)start)/delta*WIDTH,HIGHT+HIGHT_MET - 2.0);
          }
        */

// Prints Kozak sequence positions

/*
        fprintf(stdout,"\nL05 0 setlinecap\n");

        for(i=0;i<nk;++i) {
            if(kozak_str[i]=='D') {
                if(kozak[i]%period==1) fprintf(stdout,"B");
                else if(kozak[i]%period==2) fprintf(stdout,"R");
                else if(kozak[i]%period==0) fprintf(stdout,"G");
                fprintf(stdout," %.1f %.2f M 0 3 RL 2 0 RL stroke\n",(kozak[i]-(float)start)/delta*WIDTH,HIGHT);
            }
            else {
                if(kozak[i]%period==1) fprintf(stdout,"G");
                else if(kozak[i]%period==2) fprintf(stdout,"B");
                else if(kozak[i]%period==0) fprintf(stdout,"R");
                fprintf(stdout," %.1f %.2f M 0 3 RL -2 0 RL stroke\n",(kozak[i]-(float)start)/delta*WIDTH,HIGHT);
            }
        }
*/

// Prints PALINDROMES

/*
        fprintf(stdout,"\nL05 0 setlinecap Black\n");

        for(i=0;i<npali;++i)
            fprintf(stdout," %.1f %.2f M %d hairpin stroke\n",((float)pali[i][0]-(float)start)/delta*WIDTH,HIGHT,pali[i][1]);
*/

// Prints TATABOX

/*
        fprintf(stdout,"\nL025 0 setlinecap Black\n");

        for(i=0;i<nt;++i)
            if(tata_str[i]=='D')
                fprintf(stdout," %.1f %.2f M 0 5 RL 2 0 RL stroke\n",(tata[i]-(float)start)/delta*WIDTH,0.0);
            else
                fprintf(stdout," %.1f %.2f M 0 5 RL -2 0 RL stroke\n",(tata[i]-(float)start)/delta*WIDTH,0.0);
*/

// Prints CCAATBOX

/*
        fprintf(stdout,"\nDarkGray\n");

        for(i=0;i<ncca;++i)
            if(ccaa_str[i]=='D')
                fprintf(stdout," %.1f %.2f M 0 3 RL 2 0 RL stroke\n",(ccaa[i]-(float)start)/delta*WIDTH,0.0);
            else
                fprintf(stdout," %.1f %.2f M 0 3 RL -2 0 RL stroke\n",(ccaa[i]-(float)start)/delta*WIDTH,0.0);
*/

// Prints CAPBOX

/*
        fprintf(stdout,"\nGray\n");

        for(i=0;i<ncap;++i)
            if(cap_str[i]=='D')
                fprintf(stdout," %.1f %.2f M 0 4 RL 2 0 RL stroke\n",(cap[i]-(float)start)/delta*WIDTH,0.0);
            else
                fprintf(stdout," %.1f %.2f M 0 4 RL -2 0 RL stroke\n",(cap[i]-(float)start)/delta*WIDTH,0.0);
*/

// Prints GCBOX

/*
        fprintf(stdout,"\nLightGray\n");

        for(i=0;i<ngcb;++i)
            if(gcbox_str[i]=='D')
                fprintf(stdout," %.1f %.2f M 0 2 RL 2 0 RL stroke\n",(gcbox[i]-(float)start)/delta*WIDTH,0.0);
            else
                fprintf(stdout," %.1f %.2f M 0 2 RL -2 0 RL stroke\n",(gcbox[i]-(float)start)/delta*WIDTH,0.0);
*/


        /*****************************************************/
        /* Prints blocks of homology between annotated genes */
        /*****************************************************/

/*
        fprintf(stdout,"\nL2\n");

        for(i=0;i<nB;++i) {
            if(BLOCK_str[i]=='D') {
                if(BLOCK[i][0]%period==1) fprintf(stdout,"B");
                else if(BLOCK[i][0]%period==2) fprintf(stdout,"R");
                else if(BLOCK[i][0]%period==0) fprintf(stdout,"G");
                fprintf(stdout," %.1f %.2f M %.1f %.2f L stroke\n",(BLOCK[i][0]-(float)start)/delta*WIDTH,HIGHT+HIGHT_PUB+2,(BLOCK[i][1]-(float)start)/delta*WIDTH,HIGHT+HIGHT_PUB+2);
            }
            else {
                if(BLOCK[i][1]%period==1) fprintf(stdout,"G");
                else if(BLOCK[i][1]%period==2) fprintf(stdout,"B");
                else if(BLOCK[i][1]%period==0) fprintf(stdout,"R");
                fprintf(stdout," %.1f %.2f M %.1f %.2f L stroke\n",(BLOCK[i][0]-(float)start)/delta*WIDTH,HIGHT+HIGHT_PUB-2,(BLOCK[i][1]-(float)start)/delta*WIDTH,HIGHT+HIGHT_PUB-2);
            }
        }
*/


        /*********************************************************/
        /* Prints blocks of homology between potential new genes */
        /*********************************************************/

/*
        for(i=0;i<nb;++i) {
            if(block_str[i]=='D') {
                if(block[i][0]%period==1) fprintf(stdout,"B");
                else if(block[i][0]%period==2) fprintf(stdout,"R");
                else if(block[i][0]%period==0) fprintf(stdout,"G");
                fprintf(stdout," %.1f %.2f M %.1f %.2f L stroke\n",(block[i][0]-(float)start)/delta*WIDTH,HIGHT+HIGHT_HOM+1.5,(block[i][1]-(float)start)/delta*WIDTH,HIGHT+HIGHT_HOM+1.5);
            }
            else {
                if(block[i][1]%period==1) fprintf(stdout,"G");
                else if(block[i][1]%period==2) fprintf(stdout,"B");
                else if(block[i][1]%period==0) fprintf(stdout,"R");
                fprintf(stdout," %.1f %.2f M %.1f %.2f L stroke\n",(block[i][0]-(float)start)/delta*WIDTH,HIGHT+HIGHT_HOM-1.5,(block[i][1]-(float)start)/delta*WIDTH,HIGHT+HIGHT_HOM-1.5);
            }
        }
*/


        /***************************************************/
        /* Prints blocks of high GeneMark coding potential */
        /***************************************************/

/*
        fprintf(stdout,"\nL1\n");

        for(i=0;i<ncp;++i) {
            if(codpot_col[i]==1) fprintf(stdout,"R");
            else if(codpot_col[i]==2) fprintf(stdout,"G");
            else if(codpot_col[i]==0) fprintf(stdout,"B");
            if(codpot_str[i]=='D') fprintf(stdout," %.1f %.2f M %.1f %.2f L stroke\n",(codpot[i][0]-(float)start)/delta*WIDTH,HIGHT+HIGHT_CP+1.0,(codpot[i][1]-(float)start)/delta*WIDTH,HIGHT+HIGHT_CP+1.0);
            else fprintf(stdout," %.1f %.2f M %.1f %.2f L stroke\n",(codpot[i][0]-(float)start)/delta*WIDTH,HIGHT+HIGHT_CP-1.0,(codpot[i][1]-(float)start)/delta*WIDTH,HIGHT+HIGHT_CP-1.0);
        }
*/


        /***************/
        /* Prints Hits */
        /***************/

    if(nScp)
    {
    fprintf(stdout,"\n");

        for(i= 0; i < nScp; ++i)
	{
            if(Hits_str[i] != 'R')
	    {
// 1. Line width proportional to significance level.
// 2. H-hits in full color; G-hits in light color.

	    fprintf(stdout, "L%d ", Hits_type[i * 2 + 1] + 1);

		if(Hits_type[i * 2 + 0] == 'G') fprintf(stdout, "L");

                if(Hits_col[i] == 1) fprintf(stdout, "R");
                else if(Hits_col[i] == 2) fprintf(stdout, "G");
                else if(Hits_col[i] == 0) fprintf(stdout, "B");
            }

            if(Hits_str[i]=='D')
                fprintf(stdout, " %.1f %.2f M %.1f %.2f L stroke\n", (Hits[i * 2 + 0] - (float)start) / delta * WIDTH, HIGHT + HIGHT_SCP + 2.0, (Hits[i * 2 + 1] - (float)start) / delta * WIDTH, HIGHT + HIGHT_SCP + 2.0);
            else if(Hits_str[i] == 'C')
                fprintf(stdout, " %.1f %.2f M %.1f %.2f L stroke\n", (Hits[i * 2 + 0] - (float)start) / delta * WIDTH, HIGHT + HIGHT_SCP - 2.0, (Hits[i * 2 + 1] - (float)start) / delta * WIDTH, HIGHT + HIGHT_SCP - 2.0);
            else if(Hits_str[i] == 'R') {
                fprintf(stdout, "L1 Gray");
                fprintf(stdout, " %.1f %.2f M %.1f %.2f L stroke\n", (Hits[i * 2 + 0] - (float)start) / delta * WIDTH, HIGHT + HIGHT_SCP, (Hits[i * 2 + 1] - (float)start) / delta * WIDTH, HIGHT + HIGHT_SCP);
            }
        }
    free(Hits_str); Hits_str = NULL;
    free(Hits_col); Hits_col = NULL;
    free(Hits_type); Hits_type = NULL;
    free(Hits); Hits = NULL;
    }

        fprintf(stdout,"\nL05\n");

        if(np || ne)
            fprintf(stdout,"/Helvetica-Roman findfont NamesFontSize scalefont setfont\n");


        /****************************************************/
        /* Prints accepted annotated genes with short names */
        /****************************************************/
	
        for(i=0;i<np;++i)
	{
	name_len= 2 * strlen(pub_name[i]);
	cds_len= (pub[i][1] - pub[i][0]) / delta * WIDTH;
	  if(name_len * 1.4 < cds_len)
	  {
	  name_pos= (int)(((pub[i][0] + pub[i][1]) / 2 - (float)start) * 300 / delta) - name_len + HANGOVER;

            if(pub_str[i]=='D') {
                if(period%3) fprintf(stdout,"Black");
                else if(pub[i][0]%period==1) fprintf(stdout,"B");
                else if(pub[i][0]%period==2) fprintf(stdout,"R");
                else if(pub[i][0]%period==0) fprintf(stdout,"G");
                fprintf(stdout," %.1f %.2f M %.2f Rarrow\n", (pub[i][0] - (float)start) / delta * WIDTH, HIGHT + HIGHT_PUB + 2.0, (pub[i][1] - pub[i][0]) / delta * WIDTH);
                fprintf(stdout,"Black %.1f %.2f M (%s) Cshow\n", ((pub[i][1] + pub[i][0]) / 2 - (float)start) / delta * WIDTH, HIGHT + HIGHT_PUB + (DSHORT_OFFSET), pub_name[i]);
            }
            else {
                if(period%3) fprintf(stdout,"Black");
                else if(pub[i][0]%period==1) fprintf(stdout,"R");
                else if(pub[i][0]%period==2) fprintf(stdout,"G");
                else if(pub[i][0]%period==0) fprintf(stdout,"B");
                fprintf(stdout, " %.1f %.2f M %.2f Larrow\n", (pub[i][0] - (float)start) / delta * WIDTH, HIGHT + HIGHT_PUB - 2.0, (pub[i][1] - pub[i][0]) / delta * WIDTH);
                fprintf(stdout,"Black %.1f %.2f M (%s) Cshow\n", ((pub[i][1] + pub[i][0]) / 2 - (float)start) / delta * WIDTH, HIGHT + HIGHT_PUB + (CSHORT_OFFSET), pub_name[i]);
            }
	  }
        }

        /****************************************************/
        /* Prints modified annotated genes with short names */
        /****************************************************/

        for(i= 0; i < ne; ++i)
	{
	name_len= 2 *strlen(mod_name[i]);
	cds_len= (modified[i][1] - modified[i][0]) / delta * WIDTH;
	  if(name_len * 1.4 < cds_len)
	  {
	  name_pos= (int)(((modified[i][0] + modified[i][1]) / 2 - (float)start) * 300 / delta) - name_len + HANGOVER;

            if(mod_str[i] == 'D') {
                if(period % 3) fprintf(stdout, "Black");
                else if(modified[i][0] % period == 1) fprintf(stdout, "LB");
                else if(modified[i][0] % period == 2) fprintf(stdout, "LR");
                else if(modified[i][0] % period==0) fprintf(stdout, "LG");
                fprintf(stdout," %.1f %.2f M %.2f Rarrow\n", (modified[i][0] - (float)start) / delta * WIDTH, HIGHT + HIGHT_MOD + 2.0, (modified[i][1] - modified[i][0]) / delta * WIDTH);
                fprintf(stdout,"Gray %.1f %.2f M (%s) Cshow\n", ((modified[i][1] + modified[i][0]) / 2 - (float)start) / delta * WIDTH, HIGHT + HIGHT_MOD + (DSHORT_OFFSET), mod_name[i]);
            }
            else {
                if(period % 3) fprintf(stdout, "Black");
                else if(modified[i][0] % period == 1) fprintf(stdout, "LR");
                else if(modified[i][0] % period == 2) fprintf(stdout, "LG");
                else if(modified[i][0] % period == 0) fprintf(stdout, "LB");
                fprintf(stdout," %.1f %.2f M %.2f Larrow\n", (modified[i][0] - (float)start) / delta * WIDTH, HIGHT + HIGHT_MOD - 2.0, (modified[i][1] - modified[i][0]) / delta * WIDTH);
                fprintf(stdout,"Gray %.1f %.2f M (%s) Cshow\n", ((modified[i][1] + modified[i][0]) / 2 - (float)start) / delta * WIDTH, HIGHT + HIGHT_MOD + (CSHORT_OFFSET), mod_name[i]);
            }
	  }
        }

        /**********************************************/
        /* Prints new potential ORFs with short names */
        /**********************************************/

        if(newPf)
	{
            for(i= 0; i < nnP; ++i)
	    {
	    name_len= 2 * strlen(newP_name[i]);
	    cds_len=  (newP[i][1] - newP[i][0]) / delta * WIDTH;
	      if(name_len * 1.4 < cds_len)
	      {
	      name_pos= (int)(((newP[i][0] + newP[i][1]) / 2 - (float)start) * 300 / delta) - name_len + HANGOVER;

                if(newP_str[i] == 'D') {
                    if(period % 3) fprintf(stdout,"Gray");
                    else if(newP[i][0] % period == 1) fprintf(stdout,"LB");
                    else if(newP[i][0] % period == 2) fprintf(stdout,"LR");
                    else if(newP[i][0] % period == 0) fprintf(stdout,"LG");
                    fprintf(stdout," %.1f %.2f M %.2f Rarrow\n", (newP[i][0] - (float)start) / delta * WIDTH, HIGHT + HIGHT_NEW + 2.0, (newP[i][1] - newP[i][0]) / delta * WIDTH);
                    fprintf(stdout,"DarkGray %.1f %.2f M (%s) Cshow\n", ((newP[i][1] + newP[i][0]) / 2 - (float)start) / delta * WIDTH, HIGHT + HIGHT_NEW + (DSHORT_OFFSET), newP_name[i]);
                }
                else {
                    if(period % 3) fprintf(stdout, "Gray");
                    else if(newP[i][0] % period == 1) fprintf(stdout, "LR");
                    else if(newP[i][0] % period == 2) fprintf(stdout, "LG");
                    else if(newP[i][0] % period == 0) fprintf(stdout, "LB");
                    fprintf(stdout, " %.1f %.2f M %.2f Larrow\n", (newP[i][0] - (float)start) / delta * WIDTH, HIGHT + HIGHT_NEW - 2.0, (newP[i][1] - newP[i][0]) / delta * WIDTH);
                    fprintf(stdout, "DarkGray %.1f %.2f M (%s) Cshow\n", ((newP[i][1] + newP[i][0]) / 2 - (float)start) / delta * WIDTH, HIGHT + HIGHT_NEW + (CSHORT_OFFSET), newP_name[i]);
                }
	      }
            }
        }

        /************************************/
        /* Prints new ORFs with short names */
        /************************************/

        if(newf)
	{
            fprintf(stdout,"\nL1 ");
            fprintf(stdout,"/Helvetica-Bold findfont NamesFontSize scalefont setfont\n");
            for(i= 0; i < nn; ++i)
	    {
            name_len= 2 * strlen(new_name[i]);
	    cds_len= (new[i][1] - new[i][0]) / delta * WIDTH;
	      if(name_len * 1.4 < cds_len)
	      {
              name_pos= (int)(((new[i][0] + new[i][1]) / 2 - (float)start) * 300 / delta) - name_len + HANGOVER;

                if(new_str[i] == 'D') {
                    if(period % 3) fprintf(stdout, "Black");
                    else if(new[i][0] % period == 1) fprintf(stdout, "B");
                    else if(new[i][0] % period == 2) fprintf(stdout, "R");
                    else if(new[i][0] % period == 0) fprintf(stdout, "G");
                    fprintf(stdout, " %.1f %.2f M %.2f Rarrow\n", (new[i][0] - (float)start) / delta * WIDTH, HIGHT + HIGHT_NEW + 2.0, (new[i][1] - new[i][0]) / delta * WIDTH);
                    fprintf(stdout, "Black %.1f %.2f M (%s) Cshow\n", ((new[i][1] + new[i][0]) / 2 - (float)start) / delta * WIDTH, HIGHT + HIGHT_NEW + (DSHORT_OFFSET), new_name[i]);
                }
                else {
                    if(period % 3) fprintf(stdout, "Black");
                    else if(new[i][0] % period == 1) fprintf(stdout, "R");
                    else if(new[i][0] % period == 2) fprintf(stdout, "G");
                    else if(new[i][0] % period == 0) fprintf(stdout, "B");
                    fprintf(stdout, " %.1f %.2f M %.2f Larrow\n", (new[i][0] - (float)start) / delta * WIDTH, HIGHT + HIGHT_NEW - 2.0, (new[i][1] - new[i][0]) / delta * WIDTH);
                    fprintf(stdout, "Black %.1f %.2f M (%s) Cshow\n", ((new[i][1] + new[i][0]) / 2 - (float)start) / delta * WIDTH, HIGHT + HIGHT_NEW + (CSHORT_OFFSET), new_name[i]);
                }
	      }
            }
        }

        fprintf(stdout,"\nL05\n");

        if(np || ne)
            fprintf(stdout,"/Helvetica-Roman findfont NamesFontSize scalefont setfont\n");


        /***************************************************/
        /* Prints accepted annotated genes with long names */
        /***************************************************/

	for(i= 0; i < 300 + 2 * HANGOVER; ++i) pub_line[i]= pub_cline[i]= 0;

        for(i= 0; i < np; ++i)
	{
	name_len= 2 * strlen(pub_name[i]);
	cds_len= (pub[i][1] - pub[i][0]) / delta * WIDTH;
	  if(name_len * 1.4 >= cds_len)
	  {
	  name_pos= (int)(((pub[i][0] + pub[i][1]) / 2 - (float)start) * 300 / delta) - name_len + HANGOVER;

	  d= (int)NaFS * name_yoffset(pub_line, pub_cline, name_pos, name_len, pub_str[i]);

            if(pub_str[i]=='D') {
                if(period%3) fprintf(stdout,"Black");
                else if(pub[i][0]%period==1) fprintf(stdout,"B");
                else if(pub[i][0]%period==2) fprintf(stdout,"R");
                else if(pub[i][0]%period==0) fprintf(stdout,"G");
                fprintf(stdout," %.1f %.2f M %.2f Rarrow\n", (pub[i][0] - (float)start) / delta * WIDTH, HIGHT + HIGHT_PUB + 2.0, (pub[i][1] - pub[i][0]) / delta * WIDTH);
                fprintf(stdout,"Black %.1f %.2f M (%s) Cshow\n", ((pub[i][1] + pub[i][0]) / 2 - (float)start) / delta * WIDTH, HIGHT + HIGHT_PUB + (DLONG_OFFSET) + (float)d, pub_name[i]);
            }
            else {
                if(period%3) fprintf(stdout,"Black");
                else if(pub[i][0]%period==1) fprintf(stdout,"R");
                else if(pub[i][0]%period==2) fprintf(stdout,"G");
                else if(pub[i][0]%period==0) fprintf(stdout,"B");
                fprintf(stdout, " %.1f %.2f M %.2f Larrow\n", (pub[i][0] - (float)start) / delta * WIDTH, HIGHT + HIGHT_PUB - 2.0, (pub[i][1] - pub[i][0]) / delta * WIDTH);
                fprintf(stdout,"Black %.1f %.2f M (%s) Cshow\n", ((pub[i][1] + pub[i][0]) / 2 - (float)start) / delta * WIDTH, HIGHT + HIGHT_PUB + (CLONG_OFFSET) - (float)d, pub_name[i]);
            }
	  }
        }

        /***************************************************/
        /* Prints modified annotated genes with long names */
        /***************************************************/

	for(i= 0; i < 300 + 2 * HANGOVER; ++i) new_line[i]= new_cline[i]= 0;

        for(i= 0; i < ne; ++i)
	{
	name_len= 2 * strlen(mod_name[i]);
	cds_len= (modified[i][1] - modified[i][0]) / delta * WIDTH;
	  if(name_len * 1.4 >= cds_len)
	  {
	  name_pos= (int)(((modified[i][0] + modified[i][1]) / 2 - (float)start) * 300 / delta) - name_len + HANGOVER;

	  d= (int)NaFS * name_yoffset(new_line, new_cline, name_pos, name_len, mod_str[i]);

            if(mod_str[i] == 'D') {
                if(period % 3) fprintf(stdout, "Black");
                else if(modified[i][0] % period == 1) fprintf(stdout, "LB");
                else if(modified[i][0] % period == 2) fprintf(stdout, "LR");
                else if(modified[i][0] % period==0) fprintf(stdout, "LG");
                fprintf(stdout," %.1f %.2f M %.2f Rarrow\n", (modified[i][0] - (float)start) / delta * WIDTH, HIGHT + HIGHT_MOD + 2.0, (modified[i][1] - modified[i][0]) / delta * WIDTH);
                fprintf(stdout,"Gray %.1f %.2f M (%s) Cshow\n", ((modified[i][1] + modified[i][0]) / 2 - (float)start) / delta * WIDTH, HIGHT + HIGHT_MOD + (DLONG_OFFSET) + (float)d, mod_name[i]);
            }
            else {
                if(period % 3) fprintf(stdout, "Black");
                else if(modified[i][0] % period == 1) fprintf(stdout, "LR");
                else if(modified[i][0] % period == 2) fprintf(stdout, "LG");
                else if(modified[i][0] % period == 0) fprintf(stdout, "LB");
                fprintf(stdout," %.1f %.2f M %.2f Larrow\n", (modified[i][0] - (float)start) / delta * WIDTH, HIGHT + HIGHT_MOD - 2.0, (modified[i][1] - modified[i][0]) / delta * WIDTH);
                fprintf(stdout,"Gray %.1f %.2f M (%s) Cshow\n", ((modified[i][1] + modified[i][0]) / 2 - (float)start) / delta * WIDTH, HIGHT + HIGHT_MOD + (CLONG_OFFSET) - (float)d, mod_name[i]);
            }
	  }
        }


        /*********************************************/
        /* Prints new potential ORFs with long names */
        /*********************************************/

        if(newPf)
	{
            for(i=0;i<nnP;++i)
	    {
	    name_len= 2 * strlen(newP_name[i]);
	    cds_len= (newP[i][1] - newP[i][0]) / delta * WIDTH;
	      if(name_len * 1.4 >= cds_len)
	      {
	      name_pos= (int)(((newP[i][0] + newP[i][1]) / 2 - (float)start) * 300 / delta) - name_len + HANGOVER;

	      d= (int)NaFS * name_yoffset(new_line, new_cline, name_pos, name_len, newP_str[i]);

                if(newP_str[i] == 'D') {
                    if(period % 3) fprintf(stdout,"Gray");
                    else if(newP[i][0] % period == 1) fprintf(stdout,"LB");
                    else if(newP[i][0] % period == 2) fprintf(stdout,"LR");
                    else if(newP[i][0] % period == 0) fprintf(stdout,"LG");
                    fprintf(stdout," %.1f %.2f M %.2f Rarrow\n", (newP[i][0] - (float)start) / delta * WIDTH, HIGHT + HIGHT_NEW + 2.0, (newP[i][1] - newP[i][0]) / delta * WIDTH);
                    fprintf(stdout,"DarkGray %.1f %.2f M (%s) Cshow\n", ((newP[i][1] + newP[i][0]) / 2 - (float)start) / delta * WIDTH, HIGHT + HIGHT_NEW + (DLONG_OFFSET) + (float)d, newP_name[i]);
                }
                else {
                    if(period % 3) fprintf(stdout, "Gray");
                    else if(newP[i][0] % period == 1) fprintf(stdout, "LR");
                    else if(newP[i][0] % period == 2) fprintf(stdout, "LG");
                    else if(newP[i][0] % period == 0) fprintf(stdout, "LB");
                    fprintf(stdout, " %.1f %.2f M %.2f Larrow\n", (newP[i][0] - (float)start) / delta * WIDTH, HIGHT + HIGHT_NEW - 2.0, (newP[i][1] - newP[i][0]) / delta * WIDTH);
                    fprintf(stdout, "DarkGray %.1f %.2f M (%s) Cshow\n", ((newP[i][1] + newP[i][0]) / 2 - (float)start) / delta * WIDTH, HIGHT + HIGHT_NEW + (CLONG_OFFSET) - (float)d, newP_name[i]);
                }
	      }
            }
        }

        /***********************************/
        /* Prints new ORFs with long names */
        /***********************************/

        if(newf)
	{
            fprintf(stdout,"\nL1 ");
            fprintf(stdout,"/Helvetica-Bold findfont NamesFontSize scalefont setfont\n");
            for(i= 0; i < nn; ++i)
	    {
            name_len= 2 * strlen(new_name[i]);
	    cds_len= (new[i][1] - new[i][0]) / delta * WIDTH;
	      if(name_len * 1.4 >= cds_len)
	      {
              name_pos= (int)(((new[i][0] + new[i][1]) / 2 - (float)start) * 300 / delta) - name_len + HANGOVER;

              d= (int)NaFS * name_yoffset(new_line, new_cline, name_pos, name_len, new_str[i]);

                if(new_str[i] == 'D') {
                    if(period % 3) fprintf(stdout, "Black");
                    else if(new[i][0] % period == 1) fprintf(stdout, "B");
                    else if(new[i][0] % period == 2) fprintf(stdout, "R");
                    else if(new[i][0] % period == 0) fprintf(stdout, "G");
                    fprintf(stdout, " %.1f %.2f M %.2f Rarrow\n", (new[i][0] - (float)start) / delta * WIDTH, HIGHT + HIGHT_NEW + 2.0, (new[i][1] - new[i][0]) / delta * WIDTH);
                    fprintf(stdout, "Black %.1f %.2f M (%s) Cshow\n", ((new[i][1] + new[i][0]) / 2 - (float)start) / delta * WIDTH, HIGHT + HIGHT_NEW + (DLONG_OFFSET) + (float)d, new_name[i]);
                }
                else {
                    if(period % 3) fprintf(stdout, "Black");
                    else if(new[i][0] % period == 1) fprintf(stdout, "R");
                    else if(new[i][0] % period == 2) fprintf(stdout, "G");
                    else if(new[i][0] % period == 0) fprintf(stdout, "B");
                    fprintf(stdout, " %.1f %.2f M %.2f Larrow\n", (new[i][0] - (float)start) / delta * WIDTH, HIGHT + HIGHT_NEW - 2.0, (new[i][1] - new[i][0]) / delta * WIDTH);
                    fprintf(stdout, "Black %.1f %.2f M (%s) Cshow\n", ((new[i][1] + new[i][0]) / 2 - (float)start) / delta * WIDTH, HIGHT + HIGHT_NEW + (CLONG_OFFSET) - (float)d, new_name[i]);
                }
	      }
            }
        }

    }

    /*************** File footer ************************************/

    fprintf(stdout,"grestore\n");
    fprintf(stdout,"end\n");
    fprintf(stdout,"showpage\n");

    fclose(stdout);
    exit(0);
}
