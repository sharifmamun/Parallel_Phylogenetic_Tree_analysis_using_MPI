#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "phylo.h"

/*This file is used only for handling the arguments in command line*/
static int parse(int ch, struct argument* arg);

struct argument core_args[] = {
    { 'i', SARG, &input, "The distance matrix." },
    { 'o', SARG, &output, "Name for output file"},
    {   0, NARG, NULL, NULL },
};

void options(int argc,char **argv)
{
  int ch;
  extern char *optarg;
  char* option_str;
  char* cptr;
  int str_len = 2;
  const struct argument* arg;

  for (arg = core_args; arg->c != 0; ++arg) {
    ++str_len;
    if (NARG != arg->type) ++str_len;
  }

  option_str = malloc ((str_len+1) * sizeof(char));
  option_str[str_len] = '\000';

  cptr = option_str + 2;
  for (arg = core_args; arg->c != 0; ++arg) {
    *cptr++ = arg->c;
    if (NARG != arg->type) *cptr++ = ':';
  }

// help was taken from :
// http://www.gnu.org/s/libc/manual/html_node/Example-of-Getopt.html

  while ((ch = getopt(argc, argv, option_str))!= -1) {

    if (parse(ch, core_args)) continue;
    printf("Please, see the help file.....\n");	
    exit(-1);
  }
}

static int parse(int ch, struct argument* arg)
{
  extern char *optarg;

  for (; arg->c != 0; ++arg) {
    if (ch == arg->c) {
      switch (arg->type) {
      case NARG:
        {
          *(int*)(arg->val) = 1;
          return 1;
        }
        break;
      case IARG:
        {
          int tmp = atoi(optarg);
          if (tmp <= 0) return 0;
          *(int*)(arg->val) = tmp;
          return 1;
        }
        break;
      case DARG:
        {
          double tmp = strtod(optarg, NULL);
          if (tmp <= 0) return 0;
          *(double*)(arg->val) = tmp;
          return 1;
        }
        break;
      case SARG:
        {
          char* cptr = *(char**)(arg->val);
          size_t len;
          if (cptr)
            free(cptr);
          len = strlen(optarg);
          if (len > SARG_MAXLEN) len = SARG_MAXLEN;
          cptr = malloc (len+1);
          strncpy (cptr, optarg, len+1);
          *(char**)(arg->val) = cptr;
          return 1;
        }
      }
    }
  }
}

