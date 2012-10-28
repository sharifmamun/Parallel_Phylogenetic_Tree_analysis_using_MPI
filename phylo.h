#if !defined(PHYLO_H_)
#define PHYLO_H_

#include <math.h>
#define SARG_MAXLEN 1024
#define NUM 131078

typedef enum { NARG, IARG, DARG, SARG } arg_type_t;

struct argument{
  char c;
  arg_type_t type;
  void* val;
  const char* desc;
};

struct proc_attr{
    int start;
    int end;
    int rank;
    int N;
    double* matrix;
};

typedef struct proc_attr proc_t;

struct node{
    int index;
    int is_leaf; 
    double left_distance;
    double right_distance; 
    struct node* left; 
    struct node* right; 
};

typedef struct node node_t; 

int owns(proc_t, int, int);
void store(proc_t, int, int, double);
double get(proc_t, int, int);

extern char* input;
extern char* output;

/* Process the command line options */
void options(int, char **);

#endif
