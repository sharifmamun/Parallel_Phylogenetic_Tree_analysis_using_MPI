#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mpi.h"
#include "phylo.h"

char* input = 0;
char* output = 0;

//returns only true or false child or root
int owner(proc_t proc, int i, int j) {
    return (((j >= proc.start) && (j < proc.end)) ? 1 : 0);
}

int master(int j, int P, int N) {
    int i;
    int per_proc = (int) (N / P);

    for (i = 0; i < P; i++) { 
	if ((j >= (i*per_proc)) && (j < ((i+1)*per_proc)))
	    return i;
    }
}

int num_conv(int i, int j, int N) {
    return (j*N + i);
}

void print(node_t* node) {
    if (node != NULL) {
	if (node->is_leaf) {
	    printf("%d", node->index);
	}
	else {
	    printf("(");
	    print(node->left);
	    printf(",");
	    print(node->right);
	    printf(")");
	}
    }
}

//save to memory
void save(proc_t proc, int i, int j, double d) {
    proc.matrix[num_conv(i, j - proc.start, proc.N)] = d;
}

//getting values from file
double get(proc_t proc, int i, int j) {
    return (proc.matrix[num_conv(i,j-proc.start, proc.N)]);
}

// returns the number of size to the processors 
int get_size() {

    char number[20];
    unsigned int i = 0;

    FILE* fp = fopen(input, "r");
    if (fp == NULL) {
	printf("Program can't read: %s \n", input);
	return -1;
    }

    while (fscanf(fp, "%s", number) == 1) {
	i++;
    }
    fclose(fp);
    return i;
}

//reads from the hard disk!!
int read_mat(proc_t this_proc) 
{
    double d = 0;
    FILE* fp = NULL;
    int i = 0;
    int j = 0;
    char number[20];
    
    fp = fopen(input, "r");
    if (fp == NULL) {
	printf("cannot open filename: %s for reading.\n", input);
	return -1;
    }

    /** column major layout in the file. */
    for (j = 0; j < this_proc.N; j++) {
	for (i = 0; i < this_proc.N; i++) {
	    if (fscanf(fp, "%s", number) == 1) {
		d = strtod(number, NULL);
		
		if(owner(this_proc, i, j)) {
		    save(this_proc, i, j, d);
		}
	    }
	}
    }
    fclose(fp);
    return 0;
}


int main(int argc, char** argv) {
    int my_rank = 0;
    unsigned N = 0, per_procs = 0, mat_size = 0;
    int P = 1;
    int i,j;
    double minD, tmpD;
    proc_t this_proc;
    double global_min = 0;
    int n_messages = 0;
    double start,end;
    int low_val = -1; 	// initialization for working with data
    
    node_t** tree;

    start = MPI_Wtime();

//initialization for MPI 
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &P);

    /** working with command line arguments **/
    options(argc, argv);	// options locates in phylo-handler.c
    if (my_rank == 0) {

        if ((mat_size = get_size()) < 0) {
            printf("Problem getting size of distance matrix.\n");
            exit(1);
        }
    }

    /** get the total size. **/
    MPI_Bcast(&mat_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    N = (int) sqrt(mat_size);// typecasting is used so that N can be divided by P
    per_procs = (int) (N / P);
    
    if (N % P != 0) {
	printf("N must be divisable by P.\n");
	exit(1);
    }

    this_proc.start = (per_procs*my_rank) + ((my_rank == 0) ? 0 : (N % P));
    this_proc.end = this_proc.start + per_procs;
    this_proc.rank = my_rank;
    this_proc.matrix = (double*) malloc(sizeof(double)*per_procs*N);
    this_proc.N = N;

    if (my_rank == 0) {
	printf("\t No. of matrix size:  %d\n", N);
	printf("\t Total number of processors: %d\n", P);
	printf("\t Each processor calculates:  %d\n", per_procs);
    }

    // 
    short int* valid = (short int*) malloc(sizeof(short int)*N);
    double* R = (double*) malloc(sizeof(double) * N);
    double* minimums = (double*) malloc(sizeof(double) * P);
    double D_min = 0.0;
    int* minij = (int*) malloc(sizeof(int) * 2);
    double* buf = (double*) malloc(sizeof(double) * N);
    int processors = N;

//For constructing tree help was taken from :   http://cslibrary.stanford.edu/110/BinaryTrees.pdf
    if (my_rank == 0) {
	tree = (node_t **) malloc(sizeof(node_t *) * N);

	for (i = 0; i < N; i++) {
	    node_t* leaf = (node_t*) malloc(sizeof(node_t));
	    
	    leaf->left = NULL;
	    leaf->right = NULL;
	    leaf->left_distance = 0;
	    leaf->right_distance = 0;
	    leaf->index = i + 1;
	    leaf->is_leaf = 1;
	    tree[i] = leaf;
	}
    }
  
    for (i = 0; i < N; i++) {
	valid[i] = 1;
	R[i] = 0;
	buf[i] = 0;
    }
    for (i = 0; i < P; i++) {
	minimums[i] = NUM;
    }

    //checking
    if (read_mat(this_proc) < 0) {
	printf("problem reading distance matrix.");
    }

    //calculates the neew matrix [Step 1]
    while (processors > 2) {

	for (j = this_proc.start; j < this_proc.end; j++) {
	    if (valid[j]) {
		R[j] = 0.0;
		for (i = 0; i < N; i++) {
		    if (valid[i] && i != j) {
			R[j] += get(this_proc, i, j);
		    }
		}
		R[j] /= ((double) processors - 2.0);
	    }
     }
	/**  N must be divisable by P	*/
     MPI_Allgather(&R[this_proc.start], (this_proc.end - this_proc.start), MPI_DOUBLE,R, (this_proc.end - this_proc.start), MPI_DOUBLE, MPI_COMM_WORLD);

	/** Find a minimum D_min **/
	minij[0] = minij[1] = -1;
	minD = NUM;
    
	for (j = this_proc.start; j < this_proc.end; j++) {
	    if (valid[j]) {
		for (i = j + 1; i < N; i++) {
		    if (valid[i]) {
			tmpD = get(this_proc, i, j) - (R[i] + R[j]);
			if (tmpD < minD) {
			    minij[0] = i;
			    minij[1] = j;
			    minD = tmpD;
			}
		    }
		}
	    }
	}

	//getting minimum from all processors
	minimums[this_proc.rank] = minD;
	MPI_Allgather(&minimums[this_proc.rank], 1, MPI_DOUBLE, minimums, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	global_min = NUM;
	low_val = -1;
	for (i = 0; i < P; i++) {
	    if (minimums[i] < global_min) {
		low_val = i;
		global_min = minimums[i];
	    }
	}

	// getting the lowest value [Step 2]
	global_min = (this_proc.rank == low_val) ? get(this_proc, minij[0], minij[1]) : 0;
    
	MPI_Bcast(minij, 2, MPI_INT, low_val, MPI_COMM_WORLD);
	MPI_Bcast(&global_min, 1, MPI_DOUBLE, low_val, MPI_COMM_WORLD);

	//invalidate the  i'th value
	valid[minij[0]] = 0;
    
	//calculates the distances to the new node [Step 3]
	for (j = this_proc.start; j < this_proc.end; j++) {
	    if (valid[j]) {

		double tmp = .5*(get(this_proc, minij[0], j) + get(this_proc, minij[1], j) - global_min);
		// saving the specific node
		save(this_proc, minij[1], j, tmp);
		buf[j] = tmp;
	    }
	    else {
		buf[j] = get(this_proc, minij[1], j);
	    }
	}

	int to_proc = master(minij[1], P, N);

	//gathers distances from all the processors
	MPI_Gather(&buf[this_proc.start], (this_proc.end - this_proc.start), MPI_DOUBLE, buf, (this_proc.end - this_proc.start),MPI_DOUBLE, to_proc, MPI_COMM_WORLD);

	if (to_proc == this_proc.rank) {
	    for (i = 0; i < N; i++) {
		if (valid[i]) {
		    save(this_proc, i, minij[1], buf[i]);
		}
	    }
	}
	
	if (this_proc.rank == 0) {
	   		    
            // tree starts here -- calculates distance from all other taxa [Step 4]
	    double dik = .5*(global_min + R[minij[0]] - R[minij[1]]);
	    double djk = global_min - dik;
	    
	    node_t* internal_node = (node_t*) malloc(sizeof(node_t));
	    internal_node->left = tree[minij[0]];
	    internal_node->right = tree[minij[1]];
	    internal_node->left_distance = dik;
	    internal_node->right_distance = djk;
	    internal_node->is_leaf = 0;
	    
	    /** we invalidate i - so we shove it back in j. **/
	    tree[minij[1]] = internal_node;
	}
	processors--;
    }
    
    // Final join
    int li = -1;
    int lj = -1;

    for (i = 0; i < N; i++) {
	if (valid[i]) {
	    if (li < 0) 
		li = i;
	    else 
		lj = i;
	}
    }
    /** now the owner sends it to the root. **/
    low_val = master(lj, P, N);
    if (this_proc.rank == low_val) {
	global_min = get(this_proc, li, lj);
    }
    
    MPI_Bcast(&global_min, 1, MPI_DOUBLE, low_val, MPI_COMM_WORLD);

    if (my_rank == 0) {
	/** get final distances. **/
	double dik = .5*(global_min + R[li] - R[lj]);
	double djk = global_min - dik;

        end = MPI_Wtime();
	printf("\t Total time: %g\n",(end-start)); 

	node_t* root = (node_t*) malloc(sizeof(node_t));
	
	/** finish up the tree. **/
	root->left = tree[li];
	root->right = tree[lj];
	root->left_distance = dik;
	root->right_distance = djk;
	root->is_leaf = 0;
	
	printf("\t Neighbor-Joining Tree:\n");
	print(root);	//calling function "print"
	printf("\n");
    }
    MPI_Finalize();
    
    return 0;
}
