/********************************************************************/
/*                                                                  */
/*  parse (...) :                                                   */
/*       1. Reads shortest path problem in extended DIMACS format.  */
/*       2. Prepares internal data representation #1.               */
/*                                                                  */
/********************************************************************/
       
/* files to be included: */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "nodearc.h"

/* ----------------------------------------------------------------- */
int parse_gr(int *n_ad, int *m_ad, Node **nodes_ad, Arc **arcs_ad, 
	   int *node_min_ad, char *gName )

/* all parameters are output */
// long    *n_ad;                 /* address of the number of nodes */
// long    *m_ad;                 /* address of the number of arcs */
// Node    **nodes_ad;            /* address of the array of nodes */
// Arc     **arcs_ad;             /* address of the array of arcs */
// long    *node_min_ad;          /* address of the minimal node */
// char    *gName;                /* pointer to the string with file name */

{

#define MAXLINE       100	/* max line length in the input file */
#define ARC_FIELDS      3	/* no of fields in arc input  */
#define P_FIELDS        3       /* no of fields in problem line */
#define PROBLEM_TYPE "sp"       /* name of problem type*/
#define DEFAULT_NAME "unknown"  /* default name of the problem */


int    n,                      /* internal number of nodes */
        node_min=1,               /* minimal no of node  */
        node_max=0,               /* maximal no of nodes */
       *arc_first=NULL,              /* internal array for holding
                                     - node degree
                                     - position of the first outgoing arc */
       *arc_tail=NULL,               /* internal array: tails of the arcs */

        /* temporary variables carrying no of nodes */
        head, tail, i;

int    m,                      /* internal number of arcs */
        /* temporary variable carrying no of arcs */
        last, arc_num, arc_new_num;

Node    *nodes=NULL;                 /* pointer to the node structure */

Arc     *arcs=NULL,                  /* pointer to the arc structure */
        *arc_current=NULL,
        *arc_new;

long long length=0;               /* length of the current arc */

int    no_lines=0,             /* no of current input line */
        no_plines=0,            /* no of problem-lines */
        no_alines=0;            /* no of arc-lines */

char    in_line[MAXLINE],       /* for reading input line */
        pr_type[3];             /* for reading type of the problem */

 FILE *gFile;

int        err_no;                 /* no of detected error */

/* -------------- error numbers & error messages ---------------- */
#define EN1   0
#define EN2   1
#define EN3   2
#define EN4   3
#define EN6   4
#define EN10  5
#define EN7   6
#define EN8   7
#define EN9   8
#define EN11  9
#define EN12 10
#define EN13 11
#define EN14 12
#define EN16 13
#define EN15 14
#define EN17 15
#define EN18 16
#define EN21 17
#define EN19 18
#define EN20 19
#define EN22 20

static const char *err_message[] = 
  { 
/* 0*/    "more than one problem line.",
/* 1*/    "wrong number of parameters in the problem line.",
/* 2*/    "it is not a Shortest Path problem line.",
/* 3*/    "bad value of a parameter in the problem line.",
/* 4*/    "can't obtain enough memory to solve this problem.",
/* 5*/    "more than one line with the problem name.",
/* 6*/    "can't read problem name.",
/* 7*/    "problem description must preceed source/sink description.",
/* 8*/    "this parser doesn't support multiple sources/sinks.",
/* 9*/    "wrong number of parameters in the source/sink line.",
/*10*/    "wrong value of parameters in the source/sink line.",
/*11*/    "this parser doesn't support destination description.",
/*12*/    "source/sink description must be before arc descriptions.",
/*13*/    "too many arcs in the input.",
/*14*/    "wrong number of parameters in the arc line.",
/*15*/    "wrong value of parameters in the arc line.",
/*16*/    "unknown line type in the input.",
/*17*/    "reading error.",
/*18*/    "not enough arcs in the input.",
/*20*/    "can't read anything from the input file."
  };
/* --------------------------------------------------------------- */

/* The main loop:
        -  reads the line of the input,
        -  analises its type,
        -  checks correctness of parameters,
        -  puts data to the arrays,
        -  does service functions
*/

 gFile = fopen(gName, "r");
 if (gFile == NULL) {
   fprintf(stderr, "ERROR: file %s not found\n", gName);
   exit(1);
 }

while (fgets(in_line, MAXLINE, gFile) != NULL)
  {
  no_lines ++;


  switch (in_line[0])
    {
      case 'c':                  /* skip lines with comments */
      case '\n':                 /* skip empty lines   */
      case '\0':                 /* skip empty lines at the end of file */
                break;

      case 'p':                  /* problem description      */
                if ( no_plines > 0 )
                   /* more than one problem line */
                   { err_no = EN1 ; goto error; }

                no_plines = 1;
   
                if (
        /* reading problem line: type of problem, no of nodes, no of arcs */
                    sscanf ( in_line, "%*c %2s %ld %ld", pr_type, &n, &m )
                != P_FIELDS
                   )
		    /*wrong number of parameters in the problem line*/
		    { err_no = EN2; goto error; }

                if ( strcmp ( pr_type, PROBLEM_TYPE ) )
		    /*wrong problem type*/
		    { err_no = EN3; goto error; }

                if ( n <= 0  || m <= 0 )
		    /*wrong value of no of arcs or nodes*/
		    { err_no = EN4; goto error; }

        /* allocating memory for  'nodes', 'arcs'  and internal arrays */
                nodes    = (Node*) calloc ( n+2, sizeof(Node) );
		arcs     = (Arc*)  calloc ( m+1, sizeof(Arc) );
	        arc_tail = (int*) calloc ( m,   sizeof(long) ); 
		arc_first= (int*) calloc ( n+2, sizeof(long) );
                /* arc_first [ 0 .. n+1 ] = 0 - initialized by calloc */

                if ( nodes == NULL || arcs == NULL || 
                     arc_first == NULL || arc_tail == NULL )
                    /* memory is not allocated */
		    { 
		      printf("Need %lld bytes for data and %lld bytes temp. data\n", 
			     ((long long) (n+2))*((long long) sizeof(Node))+
			     ((long long) (m+1))*((long long) sizeof(Arc)),
			     ((long long) (n+m+2))*((long long) sizeof(long)));
		      err_no = EN6; goto error; 
		    }
		     
		/* setting pointer to the current arc */
		arc_current = arcs;
                break;

      case 'a':                    /* arc description */

		if ( no_alines >= m )
                  /*too many arcs on input*/
                  { err_no = EN16; goto error; }

		if (
                    /* reading an arc description */
                    sscanf ( in_line,"%*c %ld %ld %lld",
                                      &tail, &head, &length )
                    != 3 
                   ) 
                    /* arc description is not correct */
                    { err_no = EN15; goto error; }
		if ( tail < 0  ||  tail > n  ||
                     head < 0  ||  head > n  
		   )
                    /* wrong value of nodes */
		    { err_no = EN17; goto error; }

		arc_first[tail + 1] ++; /* no of arcs outgoing from tail
                                           is stored in arc_first[tail+1] */

                /* storing information about the arc */
		arc_tail[no_alines] = tail;
		arc_current -> head = nodes + head;
		arc_current -> len  = length;

		/* searching minimumu and maximum node */
                if ( head < node_min ) node_min = head;
                if ( tail < node_min ) node_min = tail;
                if ( head > node_max ) node_max = head;
                if ( tail > node_max ) node_max = tail;

		no_alines ++;
		arc_current ++;
		break;

	default:
		/* unknown type of line */
		err_no = EN18; goto error;
		break;

    } /* end of switch */
}     /* end of input loop */

/* ----- all is red  or  error while reading ----- */ 

if ( feof (gFile) == 0 ) /* reading error */
  { err_no=EN21; goto error; } 

if ( no_lines == 0 ) /* empty input */
  { err_no = EN22; goto error; } 

if ( no_alines < m ) /* not enough arcs */
  { err_no = EN19; goto error; } 

  
/********** ordering arcs - linear time algorithm ***********/

/* first arc from the first node */
( nodes + node_min ) -> first = arcs;

/* before below loop arc_first[i+1] is the number of arcs outgoing from i;
   after this loop arc_first[i] is the position of the first 
   outgoing from node i arcs after they would be ordered;
   this value is transformed to pointer and written to node.first[i]
   */
 
for ( i = node_min + 1; i <= node_max + 1; i ++ ) 
  {
    arc_first[i]          += arc_first[i-1];
    ( nodes + i ) -> first = arcs + arc_first[i];
  }


for ( i = node_min; i < node_max; i ++ ) /* scanning all the nodes  
                                            exept the last*/
  {

    last = ( ( nodes + i + 1 ) -> first ) - arcs;
                             /* arcs outgoing from i must be cited    
                              from position arc_first[i] to the position
                              equal to initial value of arc_first[i+1]-1  */

    for ( arc_num = arc_first[i]; arc_num < last; arc_num ++ )
      { tail = arc_tail[arc_num];

	while ( tail != i )
          /* the arc no  arc_num  is not in place because arc cited here
             must go out from i;
             we'll put it to its place and continue this process
             until an arc in this position would go out from i */

	  { Arc arc_tmp;        /* used in swapping below */
            arc_new_num  = arc_first[tail];
	    arc_current  = arcs + arc_num;
	    arc_new      = arcs + arc_new_num;
	    
	    /* arc_current must be cited in the position arc_new    
	       swapping these arcs:                                 */

	    arc_tmp.head         = arc_new -> head;
	    arc_new -> head      = arc_current -> head;
	    arc_current -> head = arc_tmp.head;

	    arc_tmp.len         = arc_new -> len;
	    arc_new -> len      = arc_current -> len;
	    arc_current -> len = arc_tmp.len;

	    arc_tail[arc_num] = arc_tail[arc_new_num];

	    /* we increase arc_first[tail] but label previous position */

	    arc_tail[arc_new_num] = tail;
	    arc_first[tail] ++ ;

            tail = arc_tail[arc_num];
	  }
      }
    /* all arcs outgoing from  i  are in place */
  }       

/* -----------------------  arcs are ordered  ------------------------- */

/* assigning output values */
*m_ad = m;
*n_ad = node_max - node_min + 1;
*node_min_ad = node_min;
*nodes_ad = nodes + node_min;
*arcs_ad = arcs;

/* free internal memory */
free ( arc_first ); free ( arc_tail );

 fclose(gFile);

/* Uff! all is done */
return (0);

/* ---------------------------------- */
 error:  /* error found reading input */

printf ( "\nPrs%d: line %ld of input - %s\n", 
         err_no, no_lines, err_message[err_no] );

exit (1);

}
/* --------------------   end of parser  -------------------*/

void calculate_bucket(int bucket[], Node* nodes, int n, int Delta)
{
	for (int i = 0; i < n; ++i)
		bucket[i] = (nodes + i)->dist / Delta;
}