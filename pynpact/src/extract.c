/*  -*- c-file-style:"linux" c-basic-offset:3 tab-width:3 -*-  */
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <string.h>

int SKIP= 0;
int SKIP2= 0;
char DESCRIPTION[20];
char DESCRIPTION2[20];

main(argc,argv)
int	argc;
char	*argv[];
{
   char	longstr[2000],shortstr[200],extrastr[200],*p,*q, exon[100][30],name[50],str;
   int	i,ex,stop=0;
   FILE	*input;

   if(argc==1)
   {
      fprintf(stderr,"\nUsage: extract GenBank_file [ first_n_characters_of_gene_name_to_skip [Gene_descriptor_string1 [default: gene] [ characters_of_name2_to_skip Gene_descriptor_string2] ] ]\n");
      exit(1);
   }

   if(argc>=3) SKIP = atoi(argv[2]);
   if(argc==4) strcpy(DESCRIPTION+1,argv[3]);
   else strcpy(DESCRIPTION+1,"gene");
   if(argc==5) { strcpy(DESCRIPTION2+1,argv[4]); SKIP2= SKIP; }
   else if(argc==6)
   {
      SKIP2= atoi(argv[4]);
      strcpy(DESCRIPTION2+1,argv[5]);
   }

   DESCRIPTION[0]='/';
   strcat(DESCRIPTION,"=");

   DESCRIPTION2[0]='/';
   strcat(DESCRIPTION2,"=");

   if((input=fopen(argv[1],"r"))==NULL) {
      fprintf(stderr,"\nInput file %s not found\n",argv[1]);
      exit(1);
   }

   fprintf(stderr,"\nCDS described in line \"%s\" skipping %d+1+%d\n", DESCRIPTION, (int)strlen(DESCRIPTION), SKIP);
   fprintf(stderr,"\nCDS described in line \"%s\" skipping %d+1+%d\n", DESCRIPTION2, (int)strlen(DESCRIPTION2), SKIP2);

   while(fgets(longstr,198,input) && strncmp(longstr, "ORIGIN", 6) && !feof(input))
   {
      if(!strncmp(longstr,"     CDS",8))
      {
         ex= 0;
         if(strstr(longstr,"join"))
         {
            while(longstr[strlen(longstr)-2]==',')
            {
               fgets(extrastr,98,input);
               strcpy(longstr+strlen(longstr)-1,extrastr+21);
            }
            if(!strncmp(longstr+21,"comp",4)) { q= longstr+37; str='C'; }
            else { q= longstr+26; str=' '; }
            while(p= strchr(q,','))
            {
               p[0]='\0';
               strcpy(exon[ex],q);
               q= p+1; ++ex;
            }
            p= strchr(q,')');
            p[0]='\0';
            strcpy(exon[ex],q);
            ++ex;
         }
         strcpy(shortstr,longstr+21);
         while(fgets(longstr,198,input) && 
               !strstr(longstr,DESCRIPTION) && 
               !strstr(longstr,DESCRIPTION2) && !stop)
         {
            if(strstr(longstr,"/translation=")) {
               fprintf(stderr,"\nDescriptor line not found for CDS %s",shortstr);
               stop= 1;
            }
         }
         stop= 0;
         if(strstr(longstr,DESCRIPTION))
         {
            p= strchr(longstr,DESCRIPTION[0]);
            p += strlen(DESCRIPTION)+1+SKIP;
            if(q=strchr(p,'"')) q[0]='\0';
            else p[strlen(p)-2]='\0';

            if(q= strchr(p,' ')) q[0]='\0';
            if(q= strchr(p,';')) q[0]='\0';
            strcpy(name,p);
         }
         else if(strstr(longstr,DESCRIPTION2))
         {
            p= strchr(longstr,DESCRIPTION2[0]);
            p += strlen(DESCRIPTION2)+1+SKIP2;
            if(q=strchr(p,'"')) q[0]='\0';
            else p[strlen(p)-2]='\0';

            if(q= strchr(p,' ')) q[0]='\0';
            if(q= strchr(p,';')) q[0]='\0';
            strcpy(name,p);
         }
         else strcpy(name,"?");

         if(!ex) fprintf(stdout,"%s %s",name,shortstr);
         else
            for(i=0;i<ex;++i)
            {
               strcpy(shortstr,name);
               if(str=='C') sprintf(shortstr+strlen(shortstr),"_E%d",ex-i);
               else sprintf(shortstr+strlen(shortstr),"_E%d",i+1);
               fprintf(stdout,"%s ",shortstr);
               if(str=='C') fprintf(stdout,"complement(");
               fprintf(stdout,"%s",exon[i]);
               if(str=='C') fprintf(stdout,")\n");
               else fprintf(stdout,"\n");
            }
         fflush(stdout);
      }
   }
   exit(0);
}
