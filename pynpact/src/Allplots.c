/*  -*- c-file-style:"linux" c-basic-offset:4 tab-width:4 -*-  */
# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# define DICT_SIZE 500
# define WIDTH 450.0
# define TITLES_POSITION 780
# define FIRST_TRANS -150
# define TRANSLATE -130
# define HIGHT 50.0
# define HIGHT_PUB 22.0  // Published annotation ORFs
# define HIGHT_MOD 40.0  // Modified ORFs
# define HIGHT_NEW 40.0  // New ORFs (bold-face) and Potential new ORFs (roman)
# define HIGHT_BIA 22.0
# define HIGHT_CP 8.0
# define HIGHT_SCP 8.0
# define HIGHT_HOM 52.5
# define HIGHT_MET 10.0
# define HAIRPIN_RADIUS 2.0
# define STEM_SEPARATION 1.0

# define ANOPIAS 1

# define VIEWER "showps"
# define TIC_W 5.0
# define NFS 7.0
# define LFS 12.0
# define TFS 18.0
# define TFS2 16.0

float	MAX= 100.00;
float	TIC_X;
float	TIC_Y= 20.0;

char NAME[100];
char NUCLEOTIDES[10];

void printHelp() {
    fprintf(stderr,"\nUsage:\n\n  Allplots start interval lines [x-tics period_of_frames]\n");
    fprintf(stderr,"\nstart                 Genome interval first base.");
    fprintf(stderr,"\ninterval              Number of bases.\n");
    fprintf(stderr,"\nlines                 Number of lines on page (one page).\n");
    fprintf(stderr,"\nx-tics                Number of subdivisions.\n");
    fprintf(stderr,"\nperiod_of_frame       Number of frames.\n");

    fprintf(stderr,"\nNeeds imput file \"Allplots.def\" of the form:\n\n");
    fprintf(stderr,"Plot_title\n");
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
    fprintf(stderr,"File_list_of_nucleotides_in_100bp windows.\n");
}

int quiet = 0;
#define logmsg(level,...) if(level >= quiet) { fprintf(stderr,__VA_ARGS__); }

/*
 * Read a line of arbitrary length from a file, returning a malloced
 * null-terminated string with the \n already stripped.
 */
char * ap_getl(FILE * f) {
    size_t size = 0,last= 0, len = 0;
    char * buf  = NULL;
    if(feof(f)) 
        return NULL;
    do {
        size += BUFSIZ; /* BUFSIZ is defined as "the optimal read size for this platform" */
        buf = (char*) realloc(buf,size); /* realloc(NULL,n) is the same as malloc(n) */
        /* Actually do the read. Note that fgets puts a terminal '\0' on the
           end of the string, so we make sure we overwrite this */
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


main(int argc, char *argv[]) {
    int	i,j,k,n,N,lines,nub,nc,nn,nnP,np,ne,nb,ns,ncg,nB,ncp,nScp,nm,nk,nt,ncap,ncca,ngcb,gs,ge,lp,len,period=3,swflag=1,
        start,end,pos,tstart,line_range,unbf=0,conf=0,newf=0,newPf=0,cgf=0,cpf= 0,Scpf= 0, npali= 0;
    float	r,tpr,delta,pp,S[3];

    /* Filename buffers */
    char *unb_file, *con_file, * new_file, *newP_file, *cg_file, *pub_file,
        *mod_file, *block_file,*BLOCK_file,*codpot_file,*Scodpot_file,*met_file,
        *kozak_file,*tata_file,*pali_file,*cap_file,*ccaa_file,*gcbox_file,
        *stop_file,*CG200_file,*CG100_file;


    /* by declaring them initially to be NULL, when realloc sees it, it
       will treat it as a malloc */
    int **unb=NULL, **con=NULL, **new=NULL, **newP=NULL, **cg=NULL, **pub=NULL,
        **modified=NULL, **block=NULL, **BLOCK=NULL, **codpot=NULL, **Scodpot=NULL,
        *met=NULL, *tata=NULL, *cap=NULL, *ccaa=NULL, *gcbox=NULL, *stop=NULL,
        *kozak=NULL, **pali=NULL, *X=NULL, *x=NULL;
    float **y=NULL,**Y=NULL;

    char *unb_str=NULL, *con_str=NULL, *new_str=NULL, *newP_str=NULL, *cg_str=NULL,
        *pub_str=NULL, *mod_str=NULL, *block_str=NULL, *BLOCK_str=NULL,
        *codpot_str=NULL, *codpot_col=NULL, *Scodpot_str=NULL, *Scodpot_col=NULL,
        *met_str=NULL, *tata_str=NULL, *cap_str=NULL, *ccaa_str=NULL,
        *gcbox_str=NULL, *stop_str=NULL, *kozak_str=NULL;

    char **pub_name=NULL,**mod_name=NULL,**new_name=NULL,**newP_name=NULL;

    char	longstr[200],*p,tatastr[20],*title1,*title2,ts[2];
    FILE	*input,*files;

    if(argc==1)
    {
        printHelp();
        exit(1);
    }

    fprintf(stderr,"Starting read of Allplots.def\n");
    files= fopen("Allplots.def","r");
    if(!files) {
        printHelp();
        fprintf(stderr,"Couldn't find Allplots.def in current directory.\n");
        exit(1);
    }

    fgets(longstr,98,files);

    sscanf(longstr,"%s %d",NAME,&len);
    fprintf(stderr,"\n%s %d nt",NAME,len);
    fgets(NUCLEOTIDES,8,files);
    NUCLEOTIDES[strlen(NUCLEOTIDES)-1]='\0';
    fprintf(stderr,"\n%s",NUCLEOTIDES);
    title1 = ap_getl(files);     fprintf(stderr,"\n%s",title1);
    title2 = ap_getl(files);     fprintf(stderr,"\n%s",title2);
    /* File_of_unbiased_CDSs */
    unb_file = ap_getl(files);   fprintf(stderr,"\n%s",unb_file);
    /* File_of_conserved_CDSs */
    con_file = ap_getl(files);   fprintf(stderr,"\n%s",con_file);
    /* File_of_new_CDSs */
    new_file = ap_getl(files);   fprintf(stderr,"\n%s",new_file);
    /* File_of_potential_new_CDSs */
    newP_file = ap_getl(files);  fprintf(stderr,"\n%s",newP_file);
    /*File_of_stretches_where_CG_is_asymmetric */
    cg_file = ap_getl(files);    fprintf(stderr,"\n%s",cg_file);
    /*File_of_published_accepted_CDSs */
    pub_file = ap_getl(files);   fprintf(stderr,"\n%s",pub_file);
    /*File_of_published_rejected_CDSs */
    mod_file = ap_getl(files);   fprintf(stderr,"\n%s",mod_file);
    /*File_of_blocks_from_new_ORFs_as_cds */
    block_file = ap_getl(files); fprintf(stderr,"\n%s",block_file);
    /*File_of_blocks_from_annotated_genes_as_cds */
    BLOCK_file = ap_getl(files); fprintf(stderr,"\n%s",BLOCK_file);
    /*File_of_GeneMark_regions */
    codpot_file = ap_getl(files); fprintf(stderr,"\n%s",codpot_file);
    /*File_of_G+C_coding_potential_regions */
    Scodpot_file = ap_getl(files); fprintf(stderr,"\n%s",Scodpot_file);
    /*File_of_met_positions (e.g.:D 432) */
    met_file = ap_getl(files);   fprintf(stderr,"\n%s",met_file);
    /*File_of_stop_positions (e.g.:D 432) */
    stop_file = ap_getl(files);  fprintf(stderr,"\n%s",stop_file);
    /*File_of_tatabox_positions (e.g.:105.73 D 432 TATAAAAG) */
    tata_file = ap_getl(files);  fprintf(stderr,"\n%s",tata_file);
    /*File_of_capbox_positions */
    cap_file = ap_getl(files);   fprintf(stderr,"\n%s",cap_file);
    /*File_of_ccaatbox_positions */
    ccaa_file = ap_getl(files);  fprintf(stderr,"\n%s",ccaa_file);
    /*File_of_gcbox_positions */
    gcbox_file = ap_getl(files); fprintf(stderr,"\n%s",gcbox_file);
    /*File_of_kozak_positions */
    kozak_file = ap_getl(files); fprintf(stderr,"\n%s",kozak_file);
    /*File_of_palindrom_positions_and_size */
    pali_file = ap_getl(files);  fprintf(stderr,"\n%s",pali_file);
    /*File_list_of_nucleotides_in_200bp windows. */
    CG200_file = ap_getl(files); fprintf(stderr,"\n%s",CG200_file);
    /*File_list_of_nucleotides_in_100bp windows. */
    CG100_file = ap_getl(files); fprintf(stderr,"\n%s",CG100_file);

    fclose(files);
    fprintf(stderr,"Done reading Allplots.def\n");

    tstart= atoi(argv[1]);
    line_range= atoi(argv[2]);
    lines= atoi(argv[3]);
    line_range /= lines;
    delta= (float)(line_range);

    if(argc==5) TIC_X= atof(argv[4]);
    else        TIC_X= (float)(line_range/10);

    if(argc==6) period= atoi(argv[5]);

    fprintf(stderr,"\nData read from position %d to position %d.\n",tstart,tstart+line_range);
    fprintf(stderr,"Nucleotide frequencies computed with period %d\n",period);
    fprintf(stderr,"Plots of %d line per page, %d positions per line,tics every %.0f positions.\n",lines,line_range,TIC_X);

    for(k=0; k<lines && tstart + k*line_range<len; ++k) {
        n= N= nub= nc= nn= nnP= ncg= nm= nk= nt= ncap= ncca= ngcb= ns= np= ne= nb= nB= ncp= nScp= npali= 0;
        start= tstart + k*line_range;
        end= tstart + (k+1)*line_range;
        if(end > len) end= len;

        if(input=fopen(BLOCK_file,"r")) {
            while(fgets(longstr,198,input) && !feof(input)) {
                BLOCK_str= (char *)realloc(BLOCK_str,(nB+1)*sizeof(char));
                BLOCK= (int **)realloc(BLOCK,(nB+1)*sizeof(int *));
                BLOCK[nB]= (int *)malloc(2*sizeof(int));
                p= strchr(longstr,'.');
                ge= atoi(p+2);
                if(longstr[0]=='c') { gs= atoi(longstr+11); BLOCK_str[nB]='C'; }
                else { gs= atoi(longstr); BLOCK_str[nB]=' '; }
                if(gs>=start-line_range/50 && gs<end && ge>start && ge<=end+line_range/50) { BLOCK[nB][0]= gs; BLOCK[nB][1]= ge; ++nB; }
                else if(gs>=start && gs<end && ge>end+line_range/50) { BLOCK[nB][0]= gs; BLOCK[nB][1]= end+line_range/50; ++nB; }
                else if(ge<=end+line_range/50 && ge>start && gs<start-line_range/50)
                { BLOCK[nB][0]= start-line_range/50; BLOCK[nB][0] += gs%period-BLOCK[nB][0]%period; BLOCK[nB][1]= ge; ++nB; }
                else if(gs<start-line_range/50 && ge>end+line_range/50)
                { BLOCK[nB][0]= start-line_range/50; BLOCK[nB][0] += gs%period-BLOCK[nB][0]%period; BLOCK[nB][1]= end+line_range/50; ++nB; }
            }
            fclose(input);
        }
        else fprintf(stderr,"\nFile of blocks from annotated genes %s NOT read",BLOCK_file);

        if(input=fopen(block_file,"r")) {
            while(fgets(longstr,198,input) && !feof(input)) {
                block_str= (char *)realloc(block_str,(nb+1)*sizeof(char));
                block= (int **)realloc(block,(nb+1)*sizeof(int *));
                block[nb]= (int *)malloc(2*sizeof(int));
                p= strchr(longstr,'.');
                ge= atoi(p+2);
                if(longstr[0]=='c') { gs= atoi(longstr+11); block_str[nb]='C'; }
                else { gs= atoi(longstr); block_str[nb]=' '; }
                if(gs>=start-line_range/50 && gs<end && ge>start && ge<=end+line_range/50) { block[nb][0]= gs; block[nb][1]= ge; ++nb; }
                else if(gs>=start && gs<end && ge>end+line_range/50) { block[nb][0]= gs; block[nb][1]= end+line_range/50; ++nb; }
                else if(ge<=end+line_range/50 && ge>start && gs<start-line_range/50)
                { block[nb][0]= start-line_range/50; block[nb][0] += gs%period-block[nb][0]%period; block[nb][1]= ge; ++nb; }
                else if(gs<start-line_range/50 && ge>end+line_range/50)
                { block[nb][0]= start-line_range/50; block[nb][0] += gs%period-block[nb][0]%period; block[nb][1]= end+line_range/50; ++nb; }
            }
            fclose(input);
        }
        else fprintf(stderr,"\nFile of new blocks %s NOT read",block_file);

        if(input=fopen(codpot_file,"r")) {
            ++cpf;
            while(fgets(longstr,198,input) && !feof(input)) {
                codpot_str= (char *)realloc(codpot_str,(ncp+1)*sizeof(char));
                codpot_col= (char *)realloc(codpot_col,(ncp+1)*sizeof(char));
                codpot= (int **)realloc(codpot,(ncp+1)*sizeof(int *));
                codpot[ncp]= (int *)malloc(2*sizeof(int));
                p= strchr(longstr,'.');
                ge= atoi(p+2);
                if(longstr[0]=='c') { gs= atoi(longstr+11); codpot_str[ncp]='C'; codpot_col[ncp]= gs%period; }
                else { gs= atoi(longstr); codpot_str[ncp]=' '; codpot_col[ncp]= ge%period; }
                if(gs>=start && gs<end && ge>start && ge<=end) { codpot[ncp][0]= gs; codpot[ncp][1]= ge; ++ncp; }
                else if(gs>=start && gs<end && ge>end) { codpot[ncp][0]= gs; codpot[ncp][1]= end; ++ncp; }
                else if(ge<=end && ge>start && gs<start)
                { codpot[ncp][0]= start; codpot[ncp][0] += gs%period-codpot[ncp][0]%period; codpot[ncp][1]= ge; ++ncp; }
                else if(gs<start && ge>end)
                { codpot[ncp][0]= start; codpot[ncp][0] += gs%period-codpot[ncp][0]%period; codpot[ncp][1]= end; ++ncp; }
            }
            fclose(input);
        }
        else fprintf(stderr,"\nFile of GeneMark coding potential %s NOT read",codpot_file);

        if(input=fopen(Scodpot_file,"r")) {
            ++Scpf;
            while(fgets(longstr,198,input) && !feof(input)) {
                Scodpot_str= (char *)realloc(Scodpot_str,(nScp+1)*sizeof(char));
                Scodpot_col= (char *)realloc(Scodpot_col,(nScp+1)*sizeof(char));
                Scodpot= (int **)realloc(Scodpot,(nScp+1)*sizeof(int *));
                Scodpot[nScp]= (int *)malloc(2*sizeof(int));
                p= strchr(longstr,'.');
                ge= atoi(p+2);
                if(longstr[0]=='c') { gs= atoi(longstr+11); Scodpot_str[nScp]='C'; Scodpot_col[nScp]= gs%period; }
                else if(longstr[0]=='r') { gs= atoi(longstr+7); Scodpot_str[nScp]='R'; }
                else { gs= atoi(longstr); Scodpot_str[nScp]=' '; Scodpot_col[nScp]= ge%period; }
                if(gs>=start && gs<end && ge>start && ge<=end) { Scodpot[nScp][0]= gs; Scodpot[nScp][1]= ge; ++nScp; }
                else if(gs>=start && gs<end && ge>end) { Scodpot[nScp][0]= gs; Scodpot[nScp][1]= end; ++nScp; }
                else if(ge<=end && ge>start && gs<start)
                { Scodpot[nScp][0]= start; Scodpot[nScp][0] += gs%period-Scodpot[nScp][0]%period; Scodpot[nScp][1]= ge; ++nScp; }
                else if(gs<start && ge>end)
                { Scodpot[nScp][0]= start; Scodpot[nScp][0] += gs%period-Scodpot[nScp][0]%period; Scodpot[nScp][1]= end; ++nScp; }
            }
            fclose(input);
        }
        else fprintf(stderr,"\nFile of G+C coding potential %s NOT read",Scodpot_file);

        if(input=fopen(met_file,"r")) {
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
        else fprintf(stderr,"\nFile of Met %s NOT read",met_file);

        if(input=fopen(stop_file,"r")) {
            while(fgets(longstr,198,input) && !feof(input)) {
                gs= atoi(longstr+2);
                if(gs>=start && gs<=end)
                {
                    stop_str= (char *)realloc(stop_str,(ns+1)*sizeof(char));
                    stop= (int *)realloc(stop,(ns+1)*sizeof(int));
                    if(longstr[0]=='C') stop_str[ns]='C';
                    else stop_str[ns]=' ';
                    stop[ns]= gs;
                    ++ns;
                }
            }
            fclose(input);
        }
        else fprintf(stderr,"\nFile of Stop %s NOT read",stop_file);


        if(input=fopen(tata_file,"r")) {
            while(fgets(longstr,198,input) && !feof(input)) {
                sscanf(longstr,"%f %s %d %s",&tpr,ts,&gs,tatastr);
                if(gs>=start && gs<=end) {
                    tata_str= (char *)realloc(tata_str,(nt+1)*sizeof(char));
                    tata= (int *)realloc(tata,(nt+1)*sizeof(int));
                    if(ts[0]=='C') tata_str[nt]='C';
                    else tata_str[nt]=' ';
                    tata[nt]= gs;
                    ++nt;
                }
            }
            fclose(input);
        }
        else fprintf(stderr,"\nFile of TATA box %s NOT read",tata_file);

        if(input=fopen(cap_file,"r")) {
            while(fgets(longstr,198,input) && !feof(input)) {
                sscanf(longstr,"%f %s %d %s",&tpr,ts,&gs,tatastr);
                if(gs>=start && gs<=end) {
                    cap_str= (char *)realloc(cap_str,(ncap+1)*sizeof(char));
                    cap= (int *)realloc(cap,(ncap+1)*sizeof(int));
                    if(ts[0]=='C') cap_str[ncap]='C';
                    else cap_str[ncap]=' ';
                    cap[ncap]= gs;
                    ++ncap;
                }
            }
            fclose(input);
        }
        else fprintf(stderr,"\nFile of CAP box %s NOT read",cap_file);

        if(input=fopen(ccaa_file,"r")) {
            while(fgets(longstr,198,input) && !feof(input)) {
                sscanf(longstr,"%f %s %d %s",&tpr,ts,&gs,tatastr);
                if(gs>=start && gs<=end) {
                    ccaa_str= (char *)realloc(ccaa_str,(ncca+1)*sizeof(char));
                    ccaa= (int *)realloc(ccaa,(ncca+1)*sizeof(int));
                    if(ts[0]=='C') ccaa_str[ncca]='C';
                    else ccaa_str[ncca]=' ';
                    ccaa[ncca]= gs;
                    ++ncca;
                }
            }
            fclose(input);
        }
        else fprintf(stderr,"\nFile of CCAAT box %s NOT read",ccaa_file);

        if(input=fopen(gcbox_file,"r")) {
            while(fgets(longstr,198,input) && !feof(input)) {
                sscanf(longstr,"%f %s %d %s",&tpr,ts,&gs,tatastr);
                if(gs>=start && gs<=end) {
                    gcbox_str= (char *)realloc(gcbox_str,(ngcb+1)*sizeof(char));
                    gcbox= (int *)realloc(gcbox,(ngcb+1)*sizeof(int));
                    if(ts[0]=='C') gcbox_str[ngcb]='C';
                    else gcbox_str[ngcb]=' ';
                    gcbox[ngcb]= gs;
                    ++ngcb;
                }
            }
            fclose(input);
        }
        else fprintf(stderr,"\nFile of GC box %s NOT read",gcbox_file);

        if(input=fopen(kozak_file,"r")) {
            while(fgets(longstr,198,input) && !feof(input)) {
                sscanf(longstr,"%f %s %d %s",&tpr,ts,&gs,tatastr);
                if(gs>=start && gs<=end) {
                    kozak_str= (char *)realloc(kozak_str,(nk+1)*sizeof(char));
                    kozak= (int *)realloc(kozak,(nk+1)*sizeof(int));
                    if(ts[0]=='C') kozak_str[nk]='C';
                    else kozak_str[nk]=' ';
                    kozak[nk]= gs;
                    ++nk;
                }
            }
            fclose(input);
        }
        else fprintf(stderr,"\nFile of Kozak sequences %s NOT read",kozak_file);


        if(input=fopen(pali_file,"r")) {
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
        else fprintf(stderr,"\nFile of palindromes %s NOT read",pali_file);


        /* READS FILE OF UNBIASED ORFS */

        if(input=fopen(unb_file,"r")) {
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
		unb_str[nub]=' ';
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
        else fprintf(stderr,"\nAcc file NOT read");


        /* READS FILE OF CONSERVED ORFS */

        if(input= fopen(con_file,"r")) {
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
		con_str[nc]=' ';
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
        else fprintf(stderr,"\nConserved file NOT read");

        /* READS FILE OF NEW PROPOSED CODING REGIONS */

        if(input= fopen(new_file,"r")) {
            newf= 1;
            while(fgets(longstr,198,input) && !feof(input)) {
                new_name= (char **)realloc(new_name, (nn + 1) * sizeof(char *));
                new_name[nn]= (char *)malloc(50 * sizeof(char *));
                new_str= (char *)realloc(new_str, (nn + 1) * sizeof(char));
                new= (int **)realloc(new, (nn + 1)*sizeof(int *));
                new[nn]= (int *)malloc(2 * sizeof(int));
                strncpy(new_name[nn], longstr, 48);
                p= strchr(new_name[nn], ' ');
                p[0]= '\0';
		p= strrchr(longstr, ' '); ++p;
                p= strchr(p, '.'); p += 2;
			if(p[0] == '>' || p[0] == '<') ++p;
                ge= atoi(p);
		p= strrchr(longstr, ' '); ++p;
                if(p[0] == 'c')
		{
		p += 11;
			if(p[0] == '>' || p[0] == '<') ++p;
		gs= atoi(p);
		new_str[nn]='C';
		}
                else
		{
			if(p[0] == '>' || p[0] == '<') ++p;
		gs= atoi(p);
		new_str[nn]=' ';
		}
                if(gs >= start - line_range / 50 && gs < end && ge <= end + line_range / 50) { new[nn][0]= gs; new[nn][1]= ge; ++nn; }
                else if(gs >= start - line_range / 50 && gs < end && ge > end + line_range / 50) { new[nn][0]= gs; new[nn][1]= end + line_range / 50; ++nn; }
                else if(ge <= end + line_range / 50 && ge > start && gs < start - line_range / 50)
                { new[nn][0]= start - line_range / 50; new[nn][0] += gs % period - new[nn][0] % period; new[nn][1]= ge; ++nn; }
                else if(gs < start - line_range / 50 && ge > end + line_range / 50)
                { new[nn][0]= start - line_range / 50; new[nn][0] += gs % period - new[nn][0] % period; new[nn][1]= end + line_range / 50; ++nn; }
            }
            fclose(input);
        }
        else fprintf(stderr,"\nNew file NOT read");

        /* READS FILE OF NEW POTENTIAL CODING REGIONS */

        if(input= fopen(newP_file,"r")) { newPf= 1;
            while(fgets(longstr,198,input) && !feof(input)) {
                newP_name= (char **)realloc(newP_name, (nnP + 1) * sizeof(char *));
                newP_name[nnP]= (char *)malloc(50 * sizeof(char *));
                newP_str= (char *)realloc(newP_str, (nnP + 1) * sizeof(char));
                newP= (int **)realloc(newP, (nnP + 1) * sizeof(int *));
                newP[nnP]= (int *)malloc(2 * sizeof(int));
                strncpy(newP_name[nnP], longstr, 48);
                p= strchr(newP_name[nnP], ' ');
                p[0]= '\0';
		p= strrchr(longstr, ' '); ++p;
                p= strchr(p, '.');
                ge= atoi(p + 2);
		p= strrchr(longstr, ' '); ++p;
                if(p[0] == 'c')
		{
		p += 11;
		gs= atoi(p);
		newP_str[nnP]= 'C';
		}
                else
		{
		gs= atoi(p);
		newP_str[nnP]= ' ';
		}
                if(gs >= start - line_range / 50 && gs < end && ge <= end + line_range / 50) { newP[nnP][0]= gs; newP[nnP][1]= ge; ++nnP; }
                else if(gs >= start - line_range / 50 && gs < end && ge > end + line_range / 50) { newP[nnP][0]= gs; newP[nnP][1]= end + line_range / 50; ++nnP; }
                else if(ge <= end + line_range / 50 && ge > start && gs < start - line_range / 50)
                { newP[nnP][0]= start - line_range / 50; newP[nnP][0] += gs % period - newP[nnP][0] % period; newP[nnP][1]= ge; ++nnP; }
                else if(gs < start - line_range / 50 && ge > end + line_range / 50)
                { newP[nnP][0]= start - line_range / 50; newP[nnP][0] += gs % period - newP[nnP][0] % period; newP[nnP][1]= end + line_range / 50; ++nnP; }
            }
            fclose(input);
        }
        else fprintf(stderr,"\nNew file NOT read");



        /* READS FILE OF BLOCKS OF CONTRASTING S_PATTERNS */

        if(input= fopen(cg_file,"r")) { cgf= 1;
            while(fgets(longstr,198,input) && !feof(input)) {
                cg_str= (char *)realloc(cg_str,(ncg+1)*sizeof(char));
                cg= (int **)realloc(cg,(ncg+1)*sizeof(int *));
                cg[ncg]= (int *)malloc(2*sizeof(int));
                p= strchr(longstr,'.');
                ge= atoi(p+2)+50;
                if(longstr[0]=='c') { gs= atoi(longstr+11)-50; cg_str[ncg]='C'; }
                else { gs= atoi(longstr)-50; cg_str[ncg]=' '; }
                if(gs>=start && gs<end && ge<=end) { cg[ncg][0]= gs; cg[ncg][1]= ge; ++ncg; }
                else if(gs>=start && gs<end && ge>end) { cg[ncg][0]= gs; cg[ncg][1]= end; ++ncg; }
                else if(ge<=end && ge>start && gs<start) { cg[ncg][0]= start; cg[ncg][1]= ge; ++ncg; }
                else if(gs<start && ge>end)
                { cg[ncg][0]= start; cg[ncg][1]= end; ++ncg; }
            }
            fclose(input);
        }
        else fprintf(stderr,"\nFiles with blocks of asymmetric CG content NOT read");

        /* READS FILE OF ACCEPTED PUBLIC GENES */

        if(input= fopen(pub_file,"r")) {
            logmsg(10, "Reading Pub file.\n");
            while(fgets(longstr, 198, input) && !feof(input)) {
                pub_name= (char **)realloc(pub_name, (np + 1) * sizeof(char *));
                pub_name[np]= (char *)malloc(50 * sizeof(char *));
                pub_str= (char *)realloc(pub_str, (np + 1) * sizeof(char));
                pub= (int **)realloc(pub, (np + 1) * sizeof(int *));
                pub[np]= (int *)malloc(2 * sizeof(int));
                strncpy(pub_name[np], longstr, 48);
                p= strchr(pub_name[np],' ');
                p[0]= '\0';
		p= strrchr(longstr, ' '); ++p;
                p= strchr(p, '.'); p += 2;
			if(p[0] == '>' || p[0] == '<') ++p;
                ge= atoi(p);
		p= strrchr(longstr, ' '); ++p;
                if(p[0]=='c')
		{
		p += 11;
			if(p[0] == '>' || p[0] == '<') ++p;
		gs= atoi(p);
		pub_str[np]='C';
		}
                else
		{
			if(p[0] == '>' || p[0] == '<') ++p;
		gs= atoi(p);
		pub_str[np]=' ';
		}
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
            }
            fclose(input);
        }
        else fprintf(stderr,"\nPub file NOT read");


        /* READS FILE OF MODIFIED PUBLIC GENES */

        if(input= fopen(mod_file,"r")) {
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
		mod_str[ne]=' ';
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
            fprintf(stderr,"\nModified file %s read",mod_file);
        }
        else fprintf(stderr,"\nModified file NOT read");

        if(input= fopen(CG200_file,"r")) {
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
            fprintf(stderr,"\nLarge window composition file %s read",CG200_file); }
        else fprintf(stderr,"\nLarge-window composition file NOT read");

        if(input= fopen(CG100_file,"r")) {
            while(!feof(input)) {
                fscanf(input,"%d",&pos);
                for(j=0;j<period;++j) fscanf(input,"%f",S+j);
                if(pos>=start && pos<=end) {
                    X= (int *)realloc(X,(N+1)*sizeof(int));
                    Y= (float **)realloc(Y,(N+1)*sizeof(float *));
                    Y[N]= (float *)malloc(period*sizeof(float));
                    X[N]= pos;
                    for(j=0;j<period;++j) Y[N][j]= S[j];
                    ++N;
                    pos= -1;
                }
                else fgets(longstr,198,input);
            }
            fclose(input);
            fprintf(stderr,"\nSmall-window composition file %s read",CG100_file); }
        else {
            fprintf(stderr,"\nSmall window composition file NOT read"); swflag= 0;
        }

        if(!k) {
            fprintf(stdout,"%%!PS-Adobe-2.0\n\n");
            fprintf(stdout,"gsave\n");
            fprintf(stdout,"%d dict begin\n",DICT_SIZE);
            fprintf(stdout,"/L2 { 2.0 setlinewidth } def\n");
            fprintf(stdout,"/L15 { 1.5 setlinewidth } def\n");
            fprintf(stdout,"/L1 { 1.0 setlinewidth } def\n");
            fprintf(stdout,"/L05 { 0.5 setlinewidth } def\n");
            fprintf(stdout,"/L025 { 0.25 setlinewidth } def\n");
            fprintf(stdout,"/M {moveto} def\n");
            fprintf(stdout,"/L {lineto} def\n");
            fprintf(stdout,"/RM {rmoveto} def\n");
            fprintf(stdout,"/RL {rlineto} def\n");
            fprintf(stdout,"/NumberFontSize {%1.f} def\n",NFS);
            fprintf(stdout,"/TitleFontSize {%.1f} def\n",TFS);
            fprintf(stdout,"/Title2FontSize {%.1f} def\n",TFS);
            fprintf(stdout,"/LegendFontSize {%.1f} def\n",LFS);

            if(ANOPIAS)
	    {
            fprintf(stdout,"/R {0.835 0.369 0.000 setrgbcolor} def\n");  // Vermillion
            fprintf(stdout,"/G {0.800 0.475 0.655 setrgbcolor} def\n");  // Reddish purple
            fprintf(stdout,"/B {0.000 0.447 0.698 setrgbcolor} def\n");  // Blue
            fprintf(stdout,"/LR {0.918 0.685 0.300 setrgbcolor} def\n");  // Light Vermillion
            fprintf(stdout,"/LG {0.900 0.738 0.300 setrgbcolor} def\n");  // Light Reddish purple
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
            fprintf(stdout,"/Black {0.0 0.0 0.0 setrgbcolor} def\n");
            fprintf(stdout,"/Cshow { dup stringwidth pop -2 div 0 RM show } def\n");
            fprintf(stdout,"/Lshow { dup stringwidth pop -1 mul 0 RM show } def\n");
            fprintf(stdout,"/Rarrow { 1 setlinejoin 1 setlinecap dup 4 lt { 0 3 RM dup -3 RL -1 mul -3 RL } { 0 1 RM 4 sub dup 0 RL 0 2 RL 4 -3 RL -4 -3 RL 0 2 RL -1 mul 0 RL } ifelse closepath stroke 0 setlinejoin 0 setlinecap } def\n");
            fprintf(stdout,"/Larrow { 1 setlinejoin 1 setlinecap dup 4 lt { 3 RL 0 -6 RL } { 4 3 RL 0 -2 RL dup 4 sub dup 0 RL 0 -2 RL -1 mul 0 RL 0 -2 RL } ifelse closepath stroke 0 setlinejoin 0 setlinecap } def\n");
            fprintf(stdout,"/Rhssline { dup 4 lt { 0 3 RM dup -3 RL -1 mul -3 RL } { 0 1 RM 4 sub dup 0 RL 0 2 RL 4 -3 RL -4 -3 RL 0 2 RL -1 mul 0 RL } ifelse closepath } def\n");
            fprintf(stdout,"/Lhssline { dup 4 lt { 3 RL 0 -6 RL } { 4 3 RL 0 -2 RL dup 4 sub dup 0 RL 0 -2 RL -1 mul 0 RL 0 -2 RL } ifelse closepath } def\n");
            fprintf(stdout,"/radius { %.1f } def\n",HAIRPIN_RADIUS);
            fprintf(stdout,"/thick { %.2f } def\n",STEM_SEPARATION);
            fprintf(stdout,"/hairpin { dup thick -2.0 div exch RL thick 0 RL -1 mul thick -2.0 div exch RL closepath } def\n");
            fprintf(stdout,"100 %d translate\n\n",TITLES_POSITION);


            if(atoi(argv[1])==0) {
                fprintf(stdout,"/Times-Bold findfont TitleFontSize scalefont setfont\n");
                fprintf(stdout,"%.3f %.3f M ",0.5*WIDTH-30,HIGHT-90.0);
                fprintf(stdout,"(S-profiles of %s) Cshow\n",title1);
            }
            else {
                fprintf(stdout,"/Times-Italic findfont Title2FontSize scalefont setfont\n");
                fprintf(stdout,"%.3f %.3f M ",0.0,HIGHT-90.0);
                fprintf(stdout,"(S-profiles - %s) show\n",title2);
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

        /****************************************/
        /* End printing regions of CG asymmetry */
        /****************************************/

        fprintf(stdout,"/Times-Roman findfont NumberFontSize scalefont setfont\n");
        fprintf(stdout,"L025 Gray\n");

        for(r=TIC_X/5;r<end-start;r+=TIC_X/5)
            fprintf(stdout,"%.3f 0.0 M 0 %.1f RL\n",r/delta*WIDTH,HIGHT);

        for(r=0;r<=MAX;r+=TIC_Y)
            fprintf(stdout,"0.0 %.3f M %.1f 0 RL\n",r/MAX*HIGHT,(end-start)/delta*WIDTH);

        fprintf(stdout,"stroke\nL05 Black\n");
        fprintf(stdout,"newpath 0 0 M 0 %.3f RL %.3f 0 RL 0 %.3f RL %.3f 0 RL closepath stroke\n",HIGHT,(end-start)/delta*WIDTH,-HIGHT,(start-end)/delta*WIDTH);

        for(r=0;r<=end-start;r+=TIC_X)
            fprintf(stdout,"%.3f %.3f M 0 %.3f RL 0 %.3f RM (%.0f) Cshow\n",r/delta*WIDTH,-TIC_W,TIC_W,-12.0,r+start);

        for(r=0;r<=MAX;r+=TIC_Y)
            fprintf(stdout,"%.3f %.3f M %.3f 0 RL %.3f %.3f RM (%.1f) Lshow\n",-TIC_W,r/MAX*HIGHT,TIC_W,-TIC_W-5.0,-3.0,r);


        fprintf(stdout,"stroke %.3f %.3f M (Annotation) Lshow\n",-15.0,HIGHT+HIGHT_PUB-2);
        /*
          fprintf(stdout,"stroke %.3f %.3f M (Annotation) Lshow\n",-15.0,HIGHT+HIGHT_MOD-2);
        */

        /* Prints S-profile and Genemark coding potentials (DISABLED)

           if(cpf && Scpf)
           {
           fprintf(stdout,"%.3f %.3f M (Coding Potential) Lshow\n",-25.0,HIGHT+(HIGHT_CP+HIGHT_SCP)/2-2);
           fprintf(stdout,"%.3f %.3f M (S) Lshow\n",-15.0,HIGHT+HIGHT_SCP-2);
           fprintf(stdout,"%.3f %.3f M (G) Lshow\n",-15.0,HIGHT+HIGHT_CP-2);
           fprintf(stdout,"L025 LightGray\n");
           fprintf(stdout,"0 %.3f M %.3f 0 RL stroke\n",HIGHT+HIGHT_SCP+1,(end-start)/delta*WIDTH);
           fprintf(stdout,"0 %.3f M %.3f 0 RL stroke\n",HIGHT+HIGHT_SCP-1,(end-start)/delta*WIDTH);
           fprintf(stdout,"0 %.3f M %.3f 0 RL stroke\n",HIGHT+HIGHT_CP+1,(end-start)/delta*WIDTH);
           fprintf(stdout,"0 %.3f M %.3f 0 RL stroke\n",HIGHT+HIGHT_CP-1,(end-start)/delta*WIDTH);
           fprintf(stdout," Black\n");
           }
        */

        // Prints HSSs

        if(Scpf) {
            fprintf(stdout,"%.3f %.3f M (HSSs) Lshow\n",-15.0,HIGHT+HIGHT_SCP-2);
            fprintf(stdout,"L025 LightGray\n");
            fprintf(stdout,"-4 %.3f M %.3f 0 RL -3 +3 RL stroke\n",HIGHT+HIGHT_SCP+1,(end-start)/delta*WIDTH+8);
            fprintf(stdout,"-1 %.3f -3 add M -3 3 RL %.3f 0 RL stroke\n",HIGHT+HIGHT_SCP-1,(end-start)/delta*WIDTH+8);
            fprintf(stdout," Black\n");
        }

        if(conf) {
            fprintf(stdout,"%.3f %.3f M (Conserved) Lshow\n",-25.0,HIGHT+HIGHT_HOM-2);
            fprintf(stdout," L025 LightGray\n");
            fprintf(stdout,"0 %.3f M %.3f 0 RL stroke\n",HIGHT+HIGHT_HOM+1.5,(end-start)/delta*WIDTH);
            fprintf(stdout,"0 %.3f M %.3f 0 RL stroke\n",HIGHT+HIGHT_HOM-1.5,(end-start)/delta*WIDTH);
            fprintf(stdout," Black\n");
        }
        if(newf) fprintf(stdout,"%.3f %.3f M (Newly annotated) Lshow\n",-15.0,HIGHT+HIGHT_NEW-2);

        fprintf(stdout,"/Times-Bold findfont LegendFontSize scalefont setfont Black\n");
        /*  if(k==lines-1) fprintf(stdout,"%.3f %.3f M (Genome position) Cshow\n",0.5*WIDTH,-40.0); */
        fprintf(stdout,"90 rotate %.3f %.3f M (%% %s) Cshow -90 rotate\n",0.5*HIGHT,40.0,NUCLEOTIDES);

        fprintf(stdout,"stroke DarkGray\n");

        fprintf(stdout,"%.1f %.2f M\n",(x[0]-(float)start)/delta*WIDTH,HIGHT*y[0][period]/MAX);
        for(i=1;i<n;++i)
            fprintf(stdout,"%.1f %.2f L\n",(x[i]-(float)start)/delta*WIDTH,HIGHT*y[i][period]/MAX);

        /****************************************/
        /* Prints S-profiles from small windows */
        /****************************************/

        if(swflag) {
            fprintf(stdout,"stroke L05 LR 1 setlinejoin 1 setlinecap\n");

            fprintf(stdout,"%.1f %.2f M\n",(X[0]-(float)start)/delta*WIDTH,Y[0][0]/MAX*HIGHT);
            for(i=1;i<N;++i) fprintf(stdout,"%.1f %.2f L\n",(X[i]-(float)start)/delta*WIDTH,Y[i][0]/MAX*HIGHT);

            if(period>=2) {
                fprintf(stdout,"stroke L05 LG\n");

                fprintf(stdout,"%.1f %.2f M\n",(X[0]-(float)start)/delta*WIDTH,Y[0][1]/MAX*HIGHT);
                for(i=1;i<N;++i) fprintf(stdout,"%.1f %.2f L\n",(X[i]-(float)start)/delta*WIDTH,Y[i][1]/MAX*HIGHT);

                if(period>=3) {
                    fprintf(stdout,"stroke L05 LB\n");

                    fprintf(stdout,"%.1f %.2f M\n",(X[0]-(float)start)/delta*WIDTH,Y[0][2]/MAX*HIGHT);
                    for(i=1;i<N;++i)
                        fprintf(stdout,"%.1f %.2f L\n",(X[i]-(float)start)/delta*WIDTH,Y[i][2]/MAX*HIGHT);
                }
            }
        }

        fprintf(stdout,"stroke L05 R 0 setlinejoin 0 setlinecap\n");

        /****************************************/
        /* Prints S-profiles from large windows */
        /****************************************/

        fprintf(stdout,"%.1f %.2f M 1 setlinejoin 1 setlinecap\n",(x[0]-(float)start)/delta*WIDTH,y[0][0]/MAX*HIGHT);
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
           if(unb_str[i]==' ')
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
                if(con_str[i]==' ') {
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

        /* Prints stops */

        /*

          for(i=0;i<ns;++i)
          if(stop_str[i]==' ')
          {
          if(conf)
          {
          fprintf(stdout,"%.3f %.3f M (Conserved) Lshow\n",-15.0,HIGHT+HIGHT_HOM-2);
          fprintf(stdout," L025 LightGray\n");
          fprintf(stdout,"0 %.3f M %.3f 0 RL stroke\n",HIGHT+HIGHT_HOM+1.5,(end-start)/delta*WIDTH);
          fprintf(stdout,"0 %.3f M %.3f 0 RL stroke\n",HIGHT+HIGHT_HOM-1.5,(end-start)/delta*WIDTH);
          fprintf(stdout," Black\n");
          }
          if(newf) fprintf(stdout,"%.3f %.3f M (Mrazek-Karlin annotation) Lshow\n",-15.0,HIGHT+HIGHT_NEW-2);

          fprintf(stdout,"/Times-Bold findfont LegendFontSize scalefont setfont Black\n");
          if(k==lines-1) fprintf(stdout,"%.3f %.3f M (Genome position) Cshow\n",0.5*WIDTH,-40.0);
          fprintf(stdout,"90 rotate %.3f %.3f M (%% %s) Cshow -90 rotate\n",0.5*HIGHT,40.0,NUCLEOTIDES);


          fprintf(stdout,"stroke DarkGray\n");

          fprintf(stdout,"%.1f %.2f M\n",(x[0]-(float)start)/delta*WIDTH,HIGHT*y[0][period]/MAX);
          for(i=1;i<n;++i)
          fprintf(stdout,"%.1f %.2f L\n",(x[i]-(float)start)/delta*WIDTH,HIGHT*y[i][period]/MAX);
        */

        /****************************************/
        /* Prints S-profiles from small windows */
        /****************************************/

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

        /****************************************/
        /* Prints S-profiles from large windows */
        /****************************************/

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

        fprintf(stdout,"stroke 0 setlinejoin 0 setlinecap\n");

        /****************************************/
        /* Prints new ORFs with low biases (DISABLED)

           if(unbf)
           {
           for(i=0;i<nub;++i)
           if(unb_str[i]==' ')
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
                if(con_str[i]==' ') {
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

        /* Prints stops */

        /*

          for(i=0;i<ns;++i)
          if(stop_str[i]==' ')
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

        /* Prints methionines */

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

        /* Prints Kozak sequence positions */

        fprintf(stdout,"\nL05 0 setlinecap\n");

        for(i=0;i<nk;++i) {
            if(kozak_str[i]==' ') {
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

        /* Prints PALINDROMES */

        fprintf(stdout,"\nL05 0 setlinecap Black\n");

        for(i=0;i<npali;++i)
            fprintf(stdout," %.1f %.2f M %d hairpin stroke\n",((float)pali[i][0]-(float)start)/delta*WIDTH,HIGHT,pali[i][1]);

        /* Prints TATABOX */

        fprintf(stdout,"\nL025 0 setlinecap Black\n");

        for(i=0;i<nt;++i)
            if(tata_str[i]==' ')
                fprintf(stdout," %.1f %.2f M 0 5 RL 2 0 RL stroke\n",(tata[i]-(float)start)/delta*WIDTH,0.0);
            else
                fprintf(stdout," %.1f %.2f M 0 5 RL -2 0 RL stroke\n",(tata[i]-(float)start)/delta*WIDTH,0.0);

        /* Prints CCAATBOX */

        fprintf(stdout,"\nDarkGray\n");

        for(i=0;i<ncca;++i)
            if(ccaa_str[i]==' ')
                fprintf(stdout," %.1f %.2f M 0 3 RL 2 0 RL stroke\n",(ccaa[i]-(float)start)/delta*WIDTH,0.0);
            else
                fprintf(stdout," %.1f %.2f M 0 3 RL -2 0 RL stroke\n",(ccaa[i]-(float)start)/delta*WIDTH,0.0);

        /* Prints CAPBOX */

        fprintf(stdout,"\nGray\n");

        for(i=0;i<ncap;++i)
            if(cap_str[i]==' ')
                fprintf(stdout," %.1f %.2f M 0 4 RL 2 0 RL stroke\n",(cap[i]-(float)start)/delta*WIDTH,0.0);
            else
                fprintf(stdout," %.1f %.2f M 0 4 RL -2 0 RL stroke\n",(cap[i]-(float)start)/delta*WIDTH,0.0);

        /* Prints GCBOX */

        fprintf(stdout,"\nLightGray\n");

        for(i=0;i<ngcb;++i)
            if(gcbox_str[i]==' ')
                fprintf(stdout," %.1f %.2f M 0 2 RL 2 0 RL stroke\n",(gcbox[i]-(float)start)/delta*WIDTH,0.0);
            else
                fprintf(stdout," %.1f %.2f M 0 2 RL -2 0 RL stroke\n",(gcbox[i]-(float)start)/delta*WIDTH,0.0);


        /*****************************************************/
        /* Prints blocks of homology between annotated genes */
        /*****************************************************/

        fprintf(stdout,"\nL2\n");

        for(i=0;i<nB;++i) {
            if(BLOCK_str[i]==' ') {
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


        /***************************************************************/
        /* Prints blocks of homology between potential new genes genes */
        /***************************************************************/

        for(i=0;i<nb;++i) {
            if(block_str[i]==' ') {
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


        /***************************************************/
        /* Prints blocks of high GeneMark coding potential */
        /***************************************************/

        fprintf(stdout,"\nL1\n");

        for(i=0;i<ncp;++i) {
            if(codpot_col[i]==1) fprintf(stdout,"R");
            else if(codpot_col[i]==2) fprintf(stdout,"G");
            else if(codpot_col[i]==0) fprintf(stdout,"B");
            if(codpot_str[i]==' ') fprintf(stdout," %.1f %.2f M %.1f %.2f L stroke\n",(codpot[i][0]-(float)start)/delta*WIDTH,HIGHT+HIGHT_CP+1.0,(codpot[i][1]-(float)start)/delta*WIDTH,HIGHT+HIGHT_CP+1.0);
            else fprintf(stdout," %.1f %.2f M %.1f %.2f L stroke\n",(codpot[i][0]-(float)start)/delta*WIDTH,HIGHT+HIGHT_CP-1.0,(codpot[i][1]-(float)start)/delta*WIDTH,HIGHT+HIGHT_CP-1.0);
        }


        /***************************************/
        /* Prints blocks of S coding potential */
        /***************************************/

        fprintf(stdout,"\nL1\n");

        for(i=0;i<nScp;++i) {
            if(Scodpot_str[i]!='R') {
                if(Scodpot_col[i]==1) fprintf(stdout,"R");
                else if(Scodpot_col[i]==2) fprintf(stdout,"G");
                else if(Scodpot_col[i]==0) fprintf(stdout,"B");
            }
            if(Scodpot_str[i]==' ')
                fprintf(stdout," %.1f %.2f M %.1f %.2f L stroke\n",(Scodpot[i][0]-(float)start)/delta*WIDTH,HIGHT+HIGHT_SCP+1.0,(Scodpot[i][1]-(float)start)/delta*WIDTH,HIGHT+HIGHT_SCP+1.0);
            else if(Scodpot_str[i]=='C')
                fprintf(stdout," %.1f %.2f M %.1f %.2f L stroke\n",(Scodpot[i][0]-(float)start)/delta*WIDTH,HIGHT+HIGHT_SCP-1.0,(Scodpot[i][1]-(float)start)/delta*WIDTH,HIGHT+HIGHT_SCP-1.0);
            else if(Scodpot_str[i]=='R') {
                fprintf(stdout,"Gray");
                fprintf(stdout," %.1f %.2f M %.1f %.2f L stroke\n",(Scodpot[i][0]-(float)start)/delta*WIDTH,HIGHT+HIGHT_SCP,(Scodpot[i][1]-(float)start)/delta*WIDTH,HIGHT+HIGHT_SCP);
            }
        }

        fprintf(stdout,"\nL05\n");

        if(np || ne)
            fprintf(stdout,"/Times-Roman findfont NumberFontSize scalefont setfont\n");


        /***********************************/
        /* Prints modified annotated genes */
        /***********************************/

        for(i= 0; i < ne; ++i) {
            if(mod_str[i] == ' ') {
                if(period % 3) fprintf(stdout, "Black");
                else if(modified[i][0] % period == 1) fprintf(stdout, "LB");
                else if(modified[i][0] % period == 2) fprintf(stdout, "LR");
                else if(modified[i][0] % period==0) fprintf(stdout, "LG");
                fprintf(stdout," %.1f %.2f M %.2f Rarrow\n", (modified[i][0] - (float)start) / delta * WIDTH, HIGHT + HIGHT_MOD + 2.0, (modified[i][1] - modified[i][0]) / delta * WIDTH);
                fprintf(stdout,"Gray %.1f %.2f M (%s) Cshow\n", ((modified[i][1] + modified[i][0]) / 2 - (float)start) / delta * WIDTH, HIGHT + HIGHT_MOD + 5.0, mod_name[i]);
            }
            else {
                if(period % 3) fprintf(stdout, "Black");
                else if(modified[i][0] % period == 1) fprintf(stdout, "LR");
                else if(modified[i][0] % period == 2) fprintf(stdout, "LG");
                else if(modified[i][0] % period == 0) fprintf(stdout, "LB");
                fprintf(stdout," %.1f %.2f M %.2f Larrow\n", (modified[i][0] - (float)start) / delta * WIDTH, HIGHT + HIGHT_MOD - 2.0, (modified[i][1] - modified[i][0]) / delta * WIDTH);
                fprintf(stdout,"Gray %.1f %.2f M (%s) Cshow\n", ((modified[i][1] + modified[i][0]) / 2 - (float)start) / delta * WIDTH, HIGHT + HIGHT_MOD - 10.0, mod_name[i]);
            }
        }


        /***********************************/
        /* Prints accepted annotated genes */
        /***********************************/

        for(i=0;i<np;++i) {
            if(pub_str[i]==' ') {
                if(period%3) fprintf(stdout,"Black");
                else if(pub[i][0]%period==1) fprintf(stdout,"B");
                else if(pub[i][0]%period==2) fprintf(stdout,"R");
                else if(pub[i][0]%period==0) fprintf(stdout,"G");
                fprintf(stdout," %.1f %.2f M %.2f Rarrow\n",(pub[i][0]-(float)start)/delta*WIDTH,HIGHT+HIGHT_PUB+2.0,(pub[i][1]-pub[i][0])/delta*WIDTH);
                fprintf(stdout,"Black %.1f %.2f M (%s) Cshow\n",((pub[i][1]+pub[i][0])/2-(float)start)/delta*WIDTH,HIGHT+HIGHT_PUB+5.0,pub_name[i]);
            }
            else {
                if(period%3) fprintf(stdout,"Black");
                else if(pub[i][0]%period==1) fprintf(stdout,"R");
                else if(pub[i][0]%period==2) fprintf(stdout,"G");
                else if(pub[i][0]%period==0) fprintf(stdout,"B");
                fprintf(stdout," %.1f %.2f M %.2f Larrow\n",(pub[i][0]-(float)start)/delta*WIDTH,HIGHT+HIGHT_PUB-2.0,(pub[i][1]-pub[i][0])/delta*WIDTH);
                fprintf(stdout,"Black %.1f %.2f M (%s) Cshow\n",((pub[i][1]+pub[i][0])/2-(float)start)/delta*WIDTH,HIGHT+HIGHT_PUB-10.0,pub_name[i]);
            }
        }

        /*****************************/
        /* Prints new potential ORFs */
        /*****************************/

        if(newPf) {
            for(i=0;i<nnP;++i) {
                if(newP_str[i]==' ') {
                    if(period%3) fprintf(stdout,"Gray");
                    else if(newP[i][0]%period==1) fprintf(stdout,"LB");
                    else if(newP[i][0]%period==2) fprintf(stdout,"LR");
                    else if(newP[i][0]%period==0) fprintf(stdout,"LG");
                    fprintf(stdout," %.1f %.2f M %.2f Rarrow\n",(newP[i][0]-(float)start)/delta*WIDTH,HIGHT+HIGHT_NEW+2.0,(newP[i][1]-newP[i][0])/delta*WIDTH);
                    fprintf(stdout,"DarkGray %.1f %.2f M (%s) Cshow\n",((newP[i][1]+newP[i][0])/2-(float)start)/delta*WIDTH,HIGHT+HIGHT_NEW+5.0,newP_name[i]);
                }
                else {
                    if(period%3) fprintf(stdout,"Gray");
                    else if(newP[i][0]%period==1) fprintf(stdout,"LR");
                    else if(newP[i][0]%period==2) fprintf(stdout,"LG");
                    else if(newP[i][0]%period==0) fprintf(stdout,"LB");
                    fprintf(stdout," %.1f %.2f M %.2f Larrow\n",(newP[i][0]-(float)start)/delta*WIDTH,HIGHT+HIGHT_NEW-2.0,(newP[i][1]-newP[i][0])/delta*WIDTH);
                    fprintf(stdout,"DarkGray %.1f %.2f M (%s) Cshow\n",((newP[i][1]+newP[i][0])/2-(float)start)/delta*WIDTH,HIGHT+HIGHT_NEW-10.0,newP_name[i]);
                }
            }
        }


        /*******************/
        /* Prints new ORFs */
        /*******************/

        if(newf) {
            fprintf(stdout,"\nL15 ");
            fprintf(stdout,"/Times-Bold findfont NumberFontSize scalefont setfont\n");
            for(i=0;i<nn;++i) {
                if(new_str[i]==' ') {
                    if(period%3) fprintf(stdout,"Black");
                    else if(new[i][0]%period==1) fprintf(stdout,"B");
                    else if(new[i][0]%period==2) fprintf(stdout,"R");
                    else if(new[i][0]%period==0) fprintf(stdout,"G");
                    fprintf(stdout," %.1f %.2f M %.2f Rarrow\n",(new[i][0]-(float)start)/delta*WIDTH,HIGHT+HIGHT_NEW+2.0,(new[i][1]-new[i][0])/delta*WIDTH);
                    fprintf(stdout,"Black %.1f %.2f M (%s) Cshow\n",((new[i][1]+new[i][0])/2-(float)start)/delta*WIDTH,HIGHT+HIGHT_NEW+5.0,new_name[i]);
                }
                else {
                    if(period%3) fprintf(stdout,"Black");
                    else if(new[i][0]%period==1) fprintf(stdout,"R");
                    else if(new[i][0]%period==2) fprintf(stdout,"G");
                    else if(new[i][0]%period==0) fprintf(stdout,"B");
                    fprintf(stdout," %.1f %.2f M %.2f Larrow\n",(new[i][0]-(float)start)/delta*WIDTH,HIGHT+HIGHT_NEW-2.0,(new[i][1]-new[i][0])/delta*WIDTH);
                    fprintf(stdout,"Black %.1f %.2f M (%s) Cshow\n",((new[i][1]+new[i][0])/2-(float)start)/delta*WIDTH,HIGHT+HIGHT_NEW-10.0,new_name[i]);
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
