  /* readin.c */

#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "message.h"

char *strdup2(char *c)
{
  char *d = NULL;

  if (c) {
    if ((d = strdup(c)) == NULL) fatal("Cant replicate string \"%s\"", c);
  }
  return d;
}


char **readline(FILE *f, unsigned *n, char escape)
{
  static char      *line         = NULL,
    **tag          = NULL;
  static unsigned   length       = 128,
    tag_num      = 10;

  char             *malloc_error = "Can\'t allocate storage in \"readline()\"",
    *tmp,
    **tag_tmp;
  unsigned          i;

  if (line == NULL) {
    /* If this is the first call of this function, we allocate storage for
     * several variables.  */
    line = (char*)  malloc(length + 1);
    tag  = (char**) malloc((tag_num + 1) * sizeof(char*));
    if (( ! line) || ( ! tag))
      fatal(malloc_error);
  }
  if ( ! f) {
    /* This is a special case:  If "f" is NULL, all internally allocated
     * memory is released and size delimiters are reset to their default
     * values.  This might come quite handy if the program tends to run
     * low on memory.  */
    if (line)    free(line);   line = NULL;
    if (tag)     free(tag);    tag  = NULL;
    length  = 128;
    tag_num =  10;
    return NULL;
  }
  if (fgets(line, length + 1, f)) {
    while (strchr(line, '\n') == NULL) {
      /* fgets() did not read the complete line from the input file.  So we
       * adjust the size of the line buffer (we make it twice as large) and
       * read another chunk of chars to fill it.  */
      tmp = line;
      length *= 2;
      if ((line = (char*) malloc(length + 1)) == NULL)
        fatal(malloc_error);
      strcpy(line, tmp);
      free(tmp);
      fgets(&line[length / 2], (length / 2) + 1, f);
    }
    /* Next, we look for chars which mark the beginning of comments.  If they
     * are present, the string to be examined stops there.  */
    if (escape)
      if ((tmp = strchr(line, escape)) != NULL)
        *tmp = '\0';
    /* A complete text line has been read from "f".  We will now extract
     * text tags and store pointers to them in "tag[]".  */
    tmp = line;
    *n = 0;
    while (*tmp) {
      while (isspace((int)*tmp))   /* 12.12.2005 ROHDE */
	/* skipping leading blanks */
        tmp++;
      if (*tmp) {
        /* remember the start of a tag */
        tag[*n] = tmp;
        *n += 1;
        if (*n == tag_num) {
          /* We encountered more text tags than we expected.  So we replace
           * the tag array by something twice as large.  */
          tag_tmp = tag;
          tag_num *= 2;
          if ((tag = (char**) malloc(tag_num * sizeof(char*)))  == NULL)
            fatal(malloc_error);
          for (i = 0; i < (tag_num / 2); i++)
            tag[i] = tag_tmp[i];
          free(tag_tmp);
        }
        while (*tmp && ( ! isspace((int)*tmp)))  /* 12.12.2005 ROHDE */
          /* find the end of this tag */
          tmp++;
        if (*tmp) {
          *tmp = '\0';
          tmp++;
        }
      }
    }
    tag[*n] = NULL;
    return tag;
  }
  /* If program execution reaches this part of the function, the end of file
   * "f" is reached.  No text tags are returned.  */
  *n = 0;
  return NULL;
}


void fileLength (char *name, int *persons, int *nfams, int *length, char *fnam)
{
  char     *fam_memo = NULL,
    **substring;
  FILE     *f;
  unsigned  data,
    line     = 1;

  *persons = *nfams = 0;
  if ((f = fopen(name, "r")) == NULL)
    fatal ("Cannot open input file %s", name);

  while ((substring = readline(f, &data, 0)) != NULL) {
    if (data < 3)
      fatal("Error in \"%s\": too few data columns in line %d", name, line);
    if (fam_memo == NULL) {
      fam_memo = strdup2(substring[0]);
      *length = strlen(substring[2]) + 1;
      strcpy(fnam, substring[0]);
      *nfams += 1;
    }
    else if (strcmp(fam_memo, substring[0])) {
      free(fam_memo);
      fam_memo = strdup2(substring[0]);
      *nfams += 1;
    }
    *persons += 1;
    line     += 1;
  }
  fclose(f);
}

void FileLength (char *name, int *persons, int *length)
{
  char     **substring;
  int      leng;
  FILE     *f;
  unsigned  data,
    line     = 1;

  *persons = 0;
  if ((f = fopen(name, "r")) == NULL)
    fatal ("Cannot open input file %s", name);

  while ((substring = readline(f, &data, 0)) != NULL) {
    if (data < 2)
      fatal("Error in \"%s\": too few data columns in line %d", name, line);
      if(line==1)leng = strlen(substring[1]) + 1;
      else{
      if (leng == strlen(substring[1]) + 1) *length =leng;
      else
      fatal ("Error in \"%s\": diff. genotype length in line %d",name, line);
      }
    *persons += 1;
    line     += 1;
    }
  fclose(f);
}

void SNPLength (char *name, int *persons, int *length)
{
  char     **substring;
  int      leng;
  FILE     *f;
  unsigned  data,
    line     = 1;

  *persons = 0;
  if ((f = fopen(name, "r")) == NULL)
    fatal ("Cannot open input file %s", name);

  while ((substring = readline(f, &data, 0)) != NULL) {
    if (data < 2)
      fatal("Error in \"%s\": too few data columns in line %d", name, line);
      if(line==1)leng = data;
      else{
      if (leng == data);
      else
      fatal ("Error in \"%s\": diff. sample length in line %d",name, line);
      }
    *length   = (int)data;
    *persons += 1;
    line     += 1;
    }
  fclose(f);
}
