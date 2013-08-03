    /* message.c */
    /*************/


#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include "message.h"



void fatal (char format_str[], ... )
{
/*  va_list argptr;

  fprintf(stderr, "*** Error: ");
  va_start(argptr, format_str);
  vfprintf(stderr, format_str, argptr);
  va_end(argptr);
  fprintf(stderr, "\n");
  exit(-1);
*/  
}


void warning (char format_str[], ... )
{
 /* va_list argptr;

  fprintf(stderr, "*** Warning: ");
  va_start(argptr, format_str);
  vfprintf(stderr, format_str, argptr);
  va_end(argptr);
  fprintf(stderr, "\n");
  */
}


void message (char format_str[], ... )
{
/*  va_list argptr;

  fprintf(stderr, "*** Note: ");
  va_start(argptr, format_str);
  vfprintf(stderr, format_str, argptr);
  va_end(argptr);
  fprintf(stderr, "\n");
 */
}


void print (char format_str[], ... )
{
  /*va_list argptr;

  va_start(argptr, format_str);
  vfprintf(stdout, format_str, argptr);
  va_end(argptr);
  fprintf(stdout, "\n");
  */
}



