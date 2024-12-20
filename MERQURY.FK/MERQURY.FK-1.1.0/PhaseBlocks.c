/********************************************************************************************
 *
 *  Refactoring of Merqury CN-spectra script as a command line tool using FastK
 *
 *  Author:  Gene Myers
 *  Date  :  March, 2021
 *
 *  Nancy Hansen created PhaseBlocks.c from MerquryFK.c to enable running phase block
 *   portion of the code without needing plotting software, etc. (October 2024)
 *
 *********************************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <math.h>

#define DEBUG

#include "libfastk.h"

  //  Phase Block parameters

static int ANCHOR_MARK    = 5;
static int ANCHOR_LENGTH  = 20000;

  //  Usage

static char *Usage[5] = { " [-w<double(6.0)>] [-h<double(4.5)>]",
                          " [-[xX]<number(x2.1)>] [-[yY]<number(y1.1)>]",
                          " [-vk] [-lfs] [-pdf] [-z] [-T<int(4)>] [-P<dir(/tmp)>]",
                          " <read>[.ktab] [ <mat>[.hap[.ktab]] <pat>[.hap[.ktab]] ]",
                          " <asm1:dna> [<asm2:dna>] <out>"
                        };

//  Expected inputs from FastK ...
//    MAT.hap.ktab   HAPmaker MAT PAT
//    PAT.hap.ktab
//
//  TRIO Outputs:
//    OUT.ASM[i].phased_block.hapmers.bed
//    OUT.ASM[i].phased_block.bed
//    OUT.ASM[i].phased_block.stats

static char  template[20] = "._MQY_GEN.XXXXXX";

static char  templateA1[20] = "._MQY_A1.XXXXXX";
static char  templateA2[20] = "._MQY_A2.XXXXXX";
static char *templateA[2]   = { templateA1, templateA2 };

static char  templateR1[20] = "._MQY_R1.XXXXXX";
static char  templateR2[20] = "._MQY_R2.XXXXXX";
static char *templateR[2]   = { templateR1, templateR2 };

static char  templateM1[20] = "._MQY_M1.XXXXXX";
static char  templateM2[20] = "._MQY_M2.XXXXXX";
static char *templateM[2]   = { templateM1, templateM2 };

static char  templateP1[20] = "._MQY_P1.XXXXXX";
static char  templateP2[20] = "._MQY_P2.XXXXXX";
static char *templateP[2]   = { templateP1, templateP2 };

static int   VERBOSE;
static int   KMER;

  //  Phase block codes: mark record

typedef struct
  { int beg;     //  interval [beg,end]
    int end;
    int mrk;     //  sum of majority marks  (+ = hap 1, - = hap 2)
    int opp;     //  sum of minority marks  (- = hap 1, + = hap 2)
  } Mark;

  //  mark[0..mtop) is an array of *pure* (opp=0) blocks, partition and merge blocks
  //    so that no block is unreliable save for the 1st and last, doing so in a way
  //    that minimized the impurity of the blocks.

static int merge_blocks(int mtop, Mark *mark)
{ int i, j, k;
  int last, first;

  //  First merge adjacent reliable blocks of the same polarity.

  first = -1;
  last = -1;
  j = 0;
  for (i = 0; i < mtop; i++, j++)
    { mark[j] = mark[i];
      if (mark[i].end-mark[i].beg >= ANCHOR_LENGTH || abs(mark[i].mrk) >= ANCHOR_MARK)
        { if (last >= 0 && (mark[i].mrk < 0) == (mark[last].mrk < 0))
            { for (k = last+1; k < j; k += 2)
                { mark[last].opp += mark[k].mrk;
                  mark[last].mrk += mark[k+1].mrk;
                }
              mark[last].end  = mark[j].end;
              j = last;
            }
          else
            { if (last < 0)
                first = j;
              last = j;
            }
        }
    }
  mtop = j;

  //  Now decide how to divide any unreliable blocks between opposite polarity reliable
  //    blocks to arrive at a final partitioning

  if ((mtop-last) % 2 == 0)
    { mark[mtop].beg = mark[mtop].end = mark[mtop-1].end;
      mark[mtop].mrk = mark[mtop].opp = 0;
      mtop += 1;
    }
  mark[mtop].mrk = mark[mtop].opp = 0;
  mark[mtop].end = mark[mtop-1].end;
  mark[mtop].beg = mark[mtop].end - ANCHOR_LENGTH;

  if (first % 2 != 0)
    { mark[-1].beg = mark[-1].end = mark[0].beg;
      mark[-1].mrk = mark[-1].opp = 0;
      last = -2;
    }
  else
    last = -1;

  j = 0;
  for (i = 0; i <= mtop; i++, j++)
    { mark[j] = mark[i];
      if (mark[i].end-mark[i].beg >= ANCHOR_LENGTH || abs(mark[i].mrk) >= ANCHOR_MARK)
        { int score, best, bidx;

          score = 0;
          for (k = j-1; k > last; k -= 2)
            score += mark[k].mrk;

          best = abs(score);
          bidx = k;
          for (k += 2 ; k < j; k += 2)
            { score -= (mark[k-1].mrk + mark[k].mrk);
              if (abs(score) < best)
                { best = abs(score);
                  bidx = k;
                }
              else if (abs(score) == best && last >= 0 && abs(mark[last].mrk) > abs(mark[j].mrk))
                { best = abs(score);
                  bidx = k;
                }
            }

          if (last >= 0)
            { for (k = last+2; k <= bidx; k += 2) 
                { mark[last].mrk += mark[k].mrk;
                  mark[last].opp += mark[k-1].mrk;
                }
              mark[last++].end = mark[bidx].end;
            }
          else if (bidx >= 0)
            { int tmrk, topp;

              tmrk = topp = 0;
              for (k = last+2; k <= bidx; k += 2) 
                { tmrk += mark[k].mrk;
                  topp += mark[k-1].mrk;
                }
              mark[0].mrk = tmrk;
              mark[0].opp = topp;
              mark[0].end = mark[bidx].end;
              last = 1;
            }
          else
            last = 0;
          
          for (k = bidx+2; k < j; k += 2)
            { mark[j].mrk += mark[k-1].mrk;
              mark[j].opp += mark[k].mrk;
            }
          mark[j].beg = mark[bidx+1].beg;

          mark[last] = mark[j];
          j = last;
        }
    }

  if (mark[j-1].mrk == 0)
    j -= 1;

  return (j);
}

static int64 BTOT;   //  block length sum, #, min & max
static int   NBLK;   //  set by phase_blocks, consumed by block_stats
static int   BMIN;
static int   BMAX;

  //  In a scan of an assembly's profile and relative hapmer profiles:
  //    Compute a phased-block partitioning of the contigs of the assembly
  //    and output them to <out>.<asmb>.phased_block.bed.
  //    Output temporary files containing the sizes of blocks, contigs,
  //    and scaffolds in <troot>.[blk+ctg+scf].un (these are later sorted
  //    and used to produce N-curves.

static int phase_blocks(char *aroot, char *mroot, char *proot,
                        char *aprf, char *mprf, char *pprf,
                        char *out, char *troot)
{ Profile_Index *AP, *MP, *PP;
  FILE   *bed, *fscf, *fctg, *fblk, *fmap;
  uint16 *aprof, *mprof, *pprof;
  int64   pmax, plen;
  Mark   *mark;
  int     nctg, nscf;
  int64   stot, mtot;
  int     p, i, x;

  fprintf(stderr,"\n About to open assembly profile\n");
  AP = Open_Profiles(aprf);
  fprintf(stderr,"\n About to open mat hapmer profile\n");
  MP = Open_Profiles(mprf);
  fprintf(stderr,"\n About to open pat hapmer profile\n");
  PP = Open_Profiles(pprf);
  if (AP == NULL)
    { fprintf(stderr,"\n%s: Cannot open/find FastK self-profile of %s\n",Prog_Name,aroot);
      exit (1);
    }
  if (MP == NULL)
    { fprintf(stderr,"\n%s: Cannot open/find FastK profile of %s relative to %s\n",
                     Prog_Name,aroot,mroot);
      exit (1);
    }
  if (PP == NULL)
    { fprintf(stderr,"\n%s: Cannot open/find FastK profile of %s relative to %s\n",
                     Prog_Name,aroot,proot);
      exit (1);
    }

  fprintf(stderr,"\n Assembly and hapmer profiles successfully opened?\n");
  pmax  = 20000;
  aprof = Malloc((3*pmax+2)*sizeof(uint16) + (pmax/KMER+4)*sizeof(Mark),"Profile array");
  if (aprof == NULL)
    exit (1);
  mprof = aprof + (pmax+2);
  pprof = mprof + pmax;
  mark  = ((Mark *) (pprof + pmax)) + 1;

  fscf = fopen(Catenate(troot,".scf.un","",""),"w");
  fctg = fopen(Catenate(troot,".ctg.un","",""),"w");
  fblk = fopen(Catenate(troot,".blk.un","",""),"w");

  bed  = fopen(Catenate(out,".",aroot,".phased_block.bed"),"w");
  fmap = fopen(Catenate(out,".",aroot,".hapmers.bed"),"w");
  fprintf(bed,"Scaffold\tStart\tEnd\tPhase\tPurity\tSwitches\tMarkers\n");


  stot = 0;
  mtot = 0;
  BTOT = 0;
  NBLK = 0;
  BMIN = 0x7fffffff;
  BMAX = 0;
  nctg = 0;
  nscf = AP->nreads;
  for (p = 0; p < nscf; p++)
    { int beg, end, frst;
      int mrk, sum;
      int eoc, mtop;

      sum = 0;
      mrk = 0;

      plen = Fetch_Profile(AP,p,pmax,aprof);
      if (plen > pmax)
        { pmax  = 1.2*plen + 1000;
          aprof = Realloc(aprof,(3*pmax+2)*sizeof(uint16) + (pmax/KMER+4)*sizeof(Mark),
                                "Profile array");
          if (aprof == NULL)
            exit (1);
          mprof = aprof + (pmax+2);
          pprof = mprof + pmax;
          mark  = ((Mark *) (pprof + pmax)) + 1;
          Fetch_Profile(AP,p,pmax,aprof);
        }
      Fetch_Profile(MP,p,pmax,mprof);
      Fetch_Profile(PP,p,pmax,pprof);
      aprof[plen] = 0;
      aprof[plen+1] = 1;
      fprintf(fscf,"scaffold\t%s\t%lld\n",aroot,plen);

      //  For each contig, build stack of hap sites (= overlapping hap-mers), merging adjacent
      //    sites of the same polarity and then call process_blocks to finish the partitioning
      //    into putative phased blocks

      eoc  = 0;
      end  = 0;
      frst = 0;
      mtop = 0;
      for (x = 0; x <= plen; x++)
        { int d;

          if (aprof[x] == 0)
            eoc = 1;
          else
            { if (mprof[x] > 0)
		{ d = 1;
		  fprintf(fmap,"%d\t%d\t%d\t%s.mat\n",p,x,x+KMER,aroot);
		}
              else if (pprof[x] > 0)
		{ d = -1;
		  fprintf(fmap,"%d\t%d\t%d\t%s.pat\n",p,x,x+KMER,aroot);
		}
              else
                continue;
              if (x < end)
                { mrk += d;
                  end  = x+KMER;
                  sum += 1;
                  continue;
                }
            }
          if (end > 0 && abs(mrk) >= (sum+3)/4)
            { if (mtop > 0 && (mrk < 0) == (mark[mtop-1].mrk < 0))
                { mark[mtop-1].end = end;
                  if (mrk < 0)
                    mark[mtop-1].mrk += (sum+mrk)/2 - (end-beg);
                  else
                    mark[mtop-1].mrk += (end-beg) - (sum-mrk)/2;
                }
              else
                { mark[mtop].beg = beg;
                  mark[mtop].end = end;
                  mark[mtop].opp = 0;
                  if (mrk < 0)
                    mark[mtop].mrk = (sum+mrk)/2 - (end-beg);
                  else
                    mark[mtop].mrk = (end-beg) - (sum-mrk)/2;
                  mtop += 1;
                }
            }
          if (eoc)
            { int ns, to, len;

              if (mtop > 0)
                { nctg += 1;
                  fprintf(fctg,"contig\t%s\t%d\n",aroot,(x-frst)+(KMER-1));

                  mtop = merge_blocks(mtop,mark);

                  mark[0].beg      = frst;
                  mark[mtop-1].end = x + (KMER-1);
                  for (i = 0; i < mtop; i++)
                    { len = mark[i].end - mark[i].beg;
                      if (mark[i].mrk > 0)
                        { ns = -mark[i].opp;
                          to = ns + mark[i].mrk;
                          fprintf(bed,"%d\t%d\t%d\t%s\t%.3f\t%d\t%d\n",
                                      p,mark[i].beg,mark[i].end,mroot,(100.*ns)/to,ns,to);
                          fprintf(fblk,"block\t%s\t%d\n",mroot,len);
                        }
                      else
                        { ns = mark[i].opp;
                          to = ns - mark[i].mrk;
                          fprintf(bed,"%d\t%d\t%d\t%s\t%.3f\t%d\t%d\n",
                                      nscf,mark[i].beg,mark[i].end,proot,(100.*ns)/to,ns,to);
                          fprintf(fblk,"block\t%s\t%d\n",proot,len);
                        }
                      stot += ns;
                      mtot += to;
                      BTOT += len;
                      if (len > BMAX)
                        BMAX = len;
                      else if (len < BMIN)
                        BMIN = len;
                    }
                  NBLK += mtop;
                }
               
              while (aprof[x+1] == 0)
                x += 1;
              eoc  = 0;
              end  = 0;
              mtop = 0;
              frst = x+1;
            }
          else
            { beg = x;
              end = x+KMER;
              mrk = d;
              sum = 1;
            }
        }
    }
  fprintf(bed,"total\t-t-\t-\t%.3f\t%lld\t%lld\n",
              (100.*stot)/mtot,stot,mtot);

  fclose(bed);
  fclose(fblk);
  fclose(fctg);
  fclose(fscf);

  Free_Profiles(MP);
  Free_Profiles(PP);
  Free_Profiles(AP);

  return (nctg != nscf);
}

static void block_stats(char *aroot, char *out, char *troot)
{ FILE *sfile, *nfile;
  int   bn50;
  int64 sum, thr;

  nfile = fopen(Catenate(troot,".block.sizes","",""),"r");
  sum = 0;
  thr = BTOT/2;
  while (fscanf(nfile,"block\t%*s\t%d\n",&bn50) == 1)
    { sum += bn50;
      if (sum >= thr)
        break;
    }
  fclose(nfile);
  
  sfile = fopen(Catenate(out,".",aroot,".phased_block.stats"),"w");
  fprintf(sfile,"# Blocks\tSum\tMin.\tAvg.\tN50\tMax.\n");
  fprintf(sfile,"%d\t%lld\t%d\t%lld\t%d\t%d\n",NBLK,BTOT,BMIN,BTOT/NBLK,bn50,BMAX);
  fclose(sfile);
}


/****************************************************************************************
 *
 *  Main Routine
 *
 *****************************************************************************************/

static int check_table(char *name, int lmer)
{ int   kmer;
  FILE *f;

  f = fopen(name,"r");
  if (f == NULL)
    { fprintf(stderr,"\n%s: Cannot find FastK table %s\n",Prog_Name,name);
      exit (1);
    }
  else
    { fread(&kmer,sizeof(int),1,f);
      if (lmer != 0 && kmer != lmer)
        { fprintf(stderr,"\n%s: Kmer (%d) of table %s != %d\n",Prog_Name,kmer,name,lmer);
          exit (1);
        }
      fclose(f);
      return (kmer);
    }
}

int main(int argc, char *argv[])
{ char  *ASM[2], *MAT, *PAT, *OUT;
  char  *AROOT[2], *MROOT, *PROOT;
  char  *ATAB[2], *APRF[2], *MPRF[2], *PPRF[2], *troot;
  int    LINE, FILL, STACK;
  int    NTHREADS;
  char  *SORT_PATH;

  char   command[5000];
  
  //  Command line processing

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    (void) eptr;

    ARG_INIT("PhaseBlocks");

    NTHREADS = 4;
    SORT_PATH = "/tmp";

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vlfszk")
            break;
          case 'P':
            SORT_PATH = argv[i]+2;
            break;
          case 'T':
	    ARG_POSITIVE(NTHREADS,"Number of threads")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];
    LINE    = flags['l'];
    FILL    = flags['f'];
    STACK   = flags['s'];

    if (argc < 4 || argc > 7)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
	fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[2]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[3]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[4]);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: verbose output to stderr\n");
        fprintf(stderr,"      -T: number of threads to use\n");
        fprintf(stderr,"      -P: Place all temporary files in directory -P.\n");
        exit (1);
      }

    if (LINE+FILL+STACK == 0)
      LINE = FILL = STACK = 1;

    MAT = NULL;
    PAT = NULL;
    switch (argc)
    { case 5:
        if (VERBOSE)
          fprintf(stderr,"\n Single diploid assembly with trio data\n");
        MAT    = argv[1];
        PAT    = argv[2];
        ASM[0] = argv[3];
        ASM[1] = NULL;
        break;
      case 6:
        if (VERBOSE)
          fprintf(stderr,"\n Two haploid assemblies with trio data\n");
        MAT    = argv[1];
        PAT    = argv[2];
        ASM[0] = argv[3];
        ASM[1] = argv[4];
        break;
    }
    OUT = argv[argc-1];

    //  Remove any suffixes from argument names

    { char *suffix[10] = { ".gz", ".fa", ".fq", ".fasta", ".fastq", ".db",
                           ".dam", ".sam", ".bam", ".cram" };

      troot = mktemp(template);

      char *x;

      x   = PathnRoot(MAT,".ktab");
      MAT = PathnRoot(x,".hap");
      free(x);
      x   = PathnRoot(PAT,".ktab");
      PAT = PathnRoot(x,".hap");
      free(x);

      MROOT = Root(MAT,"");
      PROOT = Root(PAT,"");
      if (strcmp(MROOT,PROOT) == 0)
        { fprintf(stderr,"%s: Parent haplotype tables have the same root name %s\n",
                         Prog_Name,MROOT);
          exit (1);
        }

      for (i = 0; i < 2; i++)
        { char *A = ASM[i];
          int   len;

          if (A == NULL)
            { AROOT[i] = NULL;
              continue;
            }

          for (j = 0; j < 10; j++)
            { len = strlen(A) - strlen(suffix[j]);
              if (strcmp(A+len,suffix[j]) == 0)
                A[len] = '\0';
            }
          AROOT[i] = Root(A,"");
        }
      if (ASM[1] != NULL && strcmp(AROOT[0],AROOT[1]) == 0)
        { fprintf(stderr,"%s: Two assemblies have the same root name %s\n",Prog_Name,AROOT[0]);
          exit (1);
        }

      KMER = check_table(Catenate(MAT,".hap",".ktab",""),0);
      KMER = check_table(Catenate(PAT,".hap",".ktab",""),KMER);

      if (VERBOSE)
        fprintf(stderr,"\n Kmer size is %d\n",KMER);

      ANCHOR_MARK *= KMER;
    }
  }

  { int    i;
    for (i = 0; i < 2; i++)
      { char  *A = ASM[i];
  
        if (A == NULL)
          continue;
  
        //  Create assembly tables and profiles
  
        ATAB[i] = mktemp(templateA[i]);
        APRF[i] = mktemp(templateR[i]);
  
        sprintf(command,"FastK -k%d -T%d -P%s -t1 -p %s -N%s",
                        KMER,NTHREADS,SORT_PATH,A,ATAB[i]); 
        fprintf(stderr,"\n Running command %s\n",command);
        SystemX(command);
  
        sprintf(command,"FastK -k%d -T%d -P%s -p:%s %s -N%s",
                        KMER,NTHREADS,SORT_PATH,ATAB[i],A,APRF[i]); 
        fprintf(stderr,"\n Running command %s\n",command);
        SystemX(command);
      }
  }

  //  Trio actions ...

  { int        i;

    //  For each assembly, determine the block phasing 

    for (i = 0; i < 2; i++)
      { char *A = ASM[i];
        char *R = AROOT[i];

        if (A == NULL)
          continue;

        if (VERBOSE)
          fprintf(stderr,"\n Producing relative profiles for phasing block calculation on %s\n",R);

        MPRF[i] = mktemp(templateM[i]);
        PPRF[i] = mktemp(templateP[i]);

        sprintf(command,"FastK -T%d -P%s -k%d -p:%s.hap %s -N%s",
                        NTHREADS,SORT_PATH,KMER,MAT,A,MPRF[i]);
        fprintf(stderr,"\n Running command %s\n",command);
        SystemX(command);

        sprintf(command,"FastK -T%d -P%s -k%d -p:%s.hap %s -N%s",
                        NTHREADS,SORT_PATH,KMER,PAT,A,PPRF[i]);
        fprintf(stderr,"\n Running command %s\n",command);
        SystemX(command);

        if (VERBOSE)
          fprintf(stderr,"\n Computing phasing blocks for assembly %s\n",R);

        phase_blocks(R, MROOT, PROOT, APRF[i], MPRF[i], PPRF[i], OUT, troot);

        if (VERBOSE)
          fprintf(stderr,"\n Sorting phase block and assembly sizes for %s\n",R);

        sprintf(command,"sort -nr -k3 %s.scf.un >%s.scaff.sizes",troot,troot);
        SystemX(command);

        sprintf(command,"sort -nr -k3 %s.ctg.un >%s.contig.sizes",troot,troot);
        SystemX(command);

        sprintf(command,"sort -nr -k3 %s.blk.un >%s.block.sizes",troot,troot);
        SystemX(command);

        //  Output block partition stats now that have block.sizes needed to compute N50

        block_stats(R,OUT,troot);

        sprintf(command,"rm %s.scf.un %s.scaff.sizes",troot,troot);
        SystemX(command);

        sprintf(command,"rm  %s.ctg.un %s.contig.sizes",troot,troot);
        SystemX(command);

        sprintf(command,"rm %s.blk.un %s.block.sizes",troot,troot);
        SystemX(command);
      }

    sprintf(command,"rm -f %s.N.R",troot);
    SystemX(command);

  }

  { int  i;

    for (i = 0; i < 2; i++)
      { char *A = ASM[i];
        char *R = AROOT[i];
 
        if (A == NULL)
          continue;

        if (MAT != NULL)
          sprintf(command,"Fastrm %s %s %s %s",ATAB[i],APRF[i],MPRF[i],PPRF[i]);
        else
          sprintf(command,"Fastrm %s %s",ATAB[i],APRF[i]);
        SystemX(command);

        free(R);
      }

    if (MAT != NULL)
      { free(PROOT);
        free(MROOT);
        free(PAT);
        free(MAT);
      }
  }

  if (VERBOSE)
    printf("\n");

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}

