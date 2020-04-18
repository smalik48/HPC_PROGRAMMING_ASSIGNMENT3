#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include <mpi.h>
#include <cstring>

using namespace std;

inline int block_decompose(const int n, const int p, const int rank)
{
    return n / p + ((rank < n % p) ? 1 : 0);
}

inline int block_decompose(const int n, MPI_Comm comm)
{
    int rank, p;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &p);
    return block_decompose(n, p, rank);
}

inline int block_decompose_by_dim(const int n, MPI_Comm comm, int dim)
{
    // get dimensions
    int dims[2];
    int periods[2];
    int coords[2];
    MPI_Cart_get(comm, 2, dims, periods, coords);
    return block_decompose(n, dims[dim], coords[dim]);
}

// Calculates y = A*x for a square n-by-n matrix A, and n-dimensional vectors x
// and y
void matrix_vector_mult(const int n, const double* A, const double* x, double* y)
{
    //TODO
    for (int j = 0; j < n; j++){
      y[j] = 0.0;
      for(int i = 0; i< n; i++)
        y[j] += A[j*n+i] * x[i];
    }
}

// Calculates y = A*x for a n-by-m matrix A, a m-dimensional vector x
// and a n-dimensional vector y
void matrix_vector_mult(const int n, const int m, const double* A, const double* x, double* y)
{
    // TODO
    for (int j = 0; j < n; j++){
      y[j] = 0.0;
      for(int i = 0; i< m; i++)
        y[j] += A[j*m+i] * x[i];
    }
}

void transpose_bcast_vector(const int n, double* col_vector, double* row_vector, MPI_Comm comm)
{
    // TODO
    //process rank in the cartesian coordinates
    int grid_rank;
    //MPI Recv status
    MPI_Status status;
    //process coordinates in cartesia
    int grid_coord[2];
    //size of processors
    int p;
    //dimension
    int dims[2];
    //period
    int periods[2];

    MPI_Comm_rank(comm, &grid_rank);
    MPI_Comm_size(comm, &p);

    //find grid_coord in the cartesian top
    MPI_Cart_get(comm, 2, dims, periods, grid_coord);

    //create Row communicator to send to diagonal elements
    MPI_Comm row_comm;
    int remain_dims[2] = {false, true};
    MPI_Cart_sub(comm, remain_dims, &row_comm);

    //get rank 00  as no send is required, just a memcpy will suffice
    int rank00;
    int coords[2] = {0, 0};
    MPI_Cart_rank(comm, coords, &rank00);
    MPI_Cart_get(comm, 2, dims, periods, grid_coord);

    //{0,0} cpy from colvector row vector.
    if(rank00 == grid_rank){
      memcpy(row_vector, col_vector, (block_decompose(n, dims[0], grid_coord[0])*sizeof(double)));
    }
    //the process is a column process that sends the data to the diagonal element using row comm
    //rank in row communicator is grid_coord[0] for the diagonal element
    else if(grid_coord[1] == 0){
      MPI_Send(col_vector, block_decompose(n, dims[0], grid_coord[0]), MPI_DOUBLE, grid_coord[0], 111, row_comm);
    }
    //the diagonal process thats going to recv from the sender(the column sender will always be zero is row_comm as the first element)
    else if(grid_coord[1] == grid_coord[0]){
      MPI_Recv(row_vector, block_decompose(n, dims[0], grid_coord[0]), MPI_DOUBLE, 0, 111, row_comm, &status);
    }

    //create column communicator for bcast
    MPI_Comm col_comm;
    remain_dims[0] = true;
    remain_dims[1] = false;
    MPI_Cart_sub(comm, remain_dims, &col_comm);

    //bcast the vector from diagonal elements to the column
    MPI_Bcast(row_vector, block_decompose(n, dims[1], grid_coord[1]), MPI_DOUBLE, grid_coord[1], col_comm);

    /*cout<<grid_coord[0]<<" "<<grid_coord[1]<<endl;
    for(int i = 0; i < block_decompose(n, dims[1], grid_coord[1]); i++){
      cout<<row_vector[i]<<" ";
    }
    cout<<endl;*/
}

void distribute_vector(const int n, double* input_vector, double** local_vector, MPI_Comm comm)
{
    // TODO
    //process rank in cartesian coordinates
    int grid_rank;
    //process coordinates
    int grid_coord[2];
    //size of processors
    int p;
    //dimension
    int dims[2];
    //period
    int periods[2];

    MPI_Comm_rank(comm, &grid_rank);
    MPI_Comm_size(comm, &p);

    //find grid_coord in the cartesian top
    MPI_Cart_get(comm, 2, dims, periods, grid_coord);



    //make a communicator in just the columns to communicate the vector
    MPI_Comm col_comm;
    int remain_dims[2] = {true, false};
    MPI_Cart_sub(comm, remain_dims, &col_comm);

    if(grid_coord[1] == 0){
      //calculate block of vector distribute in the first column
      vector<int> sendcnt(dims[0], 0);
      vector<int> displ(dims[0], 0);

      //update sendcnts and displaces cnts
      for(int i = 0; i < dims[0]; i++){
        sendcnt[i] = block_decompose(n, dims[0], i);
        if(i == 0)displ[0] = 0;
        else displ[i] = displ[i-1] + sendcnt[i-1];
      }

      //calculate recv size for first row
      int recvcnt = block_decompose_by_dim(n, comm, 0);

      //allocate size for the local vector
      (*local_vector) = new double[recvcnt];

      //find the root in first col
      int root;
      int root_coords[] = {0};
      MPI_Cart_rank(col_comm, root_coords, &root);

      //scatter value to different processors
      MPI_Scatterv(input_vector, &sendcnt[0], &displ[0], MPI_DOUBLE, *local_vector, recvcnt, MPI_DOUBLE, root, col_comm);

      //print the vector;
      /*for(int i = 0;i<recvcnt;i++){
        cout<<(*local_vector)[i]<<" ";
      }
      cout<<endl;*/
    }

}
// gather the local vector distributed among (i,0) to the processor (0,0)
void gather_vector(const int n, double* local_vector, double* output_vector, MPI_Comm comm)
{
    // TODO
    //process rank in cartesian coordinates
    int grid_rank;
    //process coordinates
    int grid_coord[2];
    //size of processors
    int p;
    //dimension
    int dims[2];
    //period
    int periods[2];

    MPI_Comm_rank(comm, &grid_rank);
    MPI_Comm_size(comm, &p);

    //find grid_coord in the cartesian top
    MPI_Cart_get(comm, 2, dims, periods, grid_coord);

    //make a communicator in just the columns to communicate the vector
    MPI_Comm col_comm;
    int remain_dims[2] = {true, false};
    MPI_Cart_sub(comm, remain_dims, &col_comm);

    if(grid_coord[1] == 0){
      //calculate block of vector distribute in the first column
      vector<int> recvcnt(dims[0], 0);
      vector<int> displ(dims[0], 0);

      //update sendcnts and displaces cnts
      for(int i = 0; i < dims[0]; i++){
        recvcnt[i] = block_decompose(n, dims[0], i);
        if(i == 0)displ[0] = 0;
        else displ[i] = displ[i-1] + recvcnt[i-1];
      }

      //calculate send size for first row
      int sendcnt = block_decompose_by_dim(n, comm, 0);

      //find the root in first col{0,0}
      int root;
      int root_coords[] = {0};
      MPI_Cart_rank(col_comm, root_coords, &root);

      // Gather values from different processors
      MPI_Gatherv(local_vector, sendcnt, MPI_DOUBLE, output_vector, &recvcnt[0], &displ[0], MPI_DOUBLE, root, col_comm);

      if(root==grid_rank){
        for(int i=0;i<n;i++){
          cout<<output_vector[i]<<" ";
        }
        cout<<endl;
      }
    }

}


void distribute_matrix(const int n, double* input_matrix, double** local_matrix, MPI_Comm comm)
{
  //TODO
  //process rank in the cartesian coordinates
  int grid_rank;
  //number of local_rows(local matrix)
  int n_local_rows = 0;
  //number of local_columns(local matrix)
  int n_local_cols = 0;
  //MPI Recv status
  MPI_Status status;
  //rank of recv process in cart common
  int dest_rank;
  //process coordinates in cartesia
  int grid_coord[2];
  //size of processors
  int p;
  //dimension
  int dims[2];
  //period
  int periods[2];
  //start_address of the block to be send
  double* start_address;

  MPI_Comm_rank(comm, &grid_rank);
  MPI_Comm_size(comm, &p);

  //find the number of rows and cols owned by each processor
  n_local_rows = block_decompose_by_dim(n, comm, 0);
  n_local_cols = block_decompose_by_dim(n, comm, 1);

  //local matrix size allocated for each processor
  /* = (double **)malloc( n_local_rows* sizeof(double *));
  for (int i=0; i<n_local_rows; i++)
    local_matrix[i] = (double *)malloc(n_local_cols * sizeof(double));*/

  // Allocate spaces for local matrix
  (*local_matrix) = new double[n_local_cols * n_local_rows];

  //get dimensions and grid_coord
  int rank00;
  int coords[2] = {0, 0};
  MPI_Cart_rank(comm, coords, &rank00);
  MPI_Cart_get(comm, 2, dims, periods, grid_coord);

 
  int row_count = 0;
 

  //iterate over grid row
  for(int i = 0; i < dims[0]; i++){
    grid_coord[0] = i;
    if(i!=0) row_count+=block_decompose(n, dims[0], i-1);
    //go over each row within a block of matrix
    for(int j = 0; j < block_decompose(n, dims[0], i); j++){
      //go over each column in grid row
      int col_count = 0;
      for(int k = 0; k < dims[1]; k++){
        grid_coord[1] = k;
        if(k!=0) col_count += block_decompose(n, dims[1], k-1);
        //calculate the destination rank where the block is supposed to be delivered
        MPI_Cart_rank(comm, grid_coord, &dest_rank);

        //see if its master
        if(rank00 == grid_rank){
          //if(i == 0 && k == 0) start_address = &input_matrix[ 0 + (j*n)];
          //else if(i==0 && k!=0) start_address = &input_matrix[0 + (j*n) + (col_count)];
	 // else if(i!=0 && k==0) start_address = &input_matrix[(row_count*n)+(j*n)];
          start_address = &input_matrix[ (row_count*n)+(j*n)+(col_count)];
          //no need to send the matrix and keep it here in new location
          if(dest_rank == rank00){
            memcpy((*local_matrix + j*n_local_cols), start_address, n_local_cols * sizeof(double));
            /*for(int i = 0; i <n_local_cols;i++){
              //local_matrix[i+(j*n_local_cols)] = start_address[i];
              (*local_matrix + j*n_local_cols)[i] = start_address[i];
            }*/
          }
          //need to send the matrix to other processor
          else{
            MPI_Send(start_address, block_decompose(n, dims[1], k), MPI_DOUBLE, dest_rank, 111, comm);
          }
        }
        //recv from zero
        else if(grid_rank == dest_rank){
          MPI_Recv((*local_matrix + j*n_local_cols), n_local_cols, MPI_DOUBLE, rank00, 111, comm, &status);
        }
      }//end for k
    }//end for j
  }//end for i

   for (int i = 0; i < n_local_rows; i++) {
        for (int j = 0; j < n_local_cols; j++) {
             std::cout << *(*local_matrix + i * n_local_cols + j) << " ";
         }
         std::cout << std::endl;
     }

}

void distributed_matrix_vector_mult(const int n, double* local_A, double* local_x, double* local_y, MPI_Comm comm)
{
    // TODO
    //process rank in the cartesian coordinates
    int grid_rank;
    //process coordinates in cartesia
    int grid_coord[2];
    //number of local_rows(local matrix)
    int n_local_rows = 0;
    //number of local_columns(local matrix)
    int n_local_cols = 0;
    //size of processors
    int p;
    //dimension
    int dims[2];
    //period
    int periods[2];

    MPI_Comm_rank(comm, &grid_rank);
    MPI_Comm_size(comm, &p);

    //find dims in the cartesian top
    MPI_Cart_get(comm, 2, dims, periods, grid_coord);

    //find the number of rows and cols owned by each processor
    n_local_rows = block_decompose_by_dim(n, comm, 0);
    n_local_cols = block_decompose_by_dim(n, comm, 1);

    //local x vector needs to be provided to each processor
    vector<double> local_x_new(n_local_cols);

    transpose_bcast_vector(n, local_x, &local_x_new[0], comm);

    vector<double> local_y_temp(n_local_rows);

    //nxn nx1-> matrix multiply
    if(n_local_cols==n_local_rows){
      matrix_vector_mult(n_local_cols, local_A, &local_x_new[0], &local_y_temp[0]);
    }
    //nxm mx1-> matrix multiply
    else{
      matrix_vector_mult(n_local_rows, n_local_cols, local_A, &local_x_new[0], &local_y_temp[0]);
    }

    //Create a row communicator to reduce the result to first columns
    MPI_Comm row_comm;
    int remain_dims[2] = {false, true};
    MPI_Cart_sub(comm, remain_dims, &row_comm);

    //get the rank of the first cols
    int first_col_rank;
    int first_col_coords[ ] = {0};
    MPI_Cart_rank(row_comm, first_col_coords, &first_col_rank);

    //reduction on first col
    MPI_Reduce(&local_y_temp[0], local_y, n_local_rows, MPI_DOUBLE, MPI_SUM, first_col_rank, row_comm);

    //print the vector
    /*if(grid_coord[1] == 0){
      for(int i=0;i<n_local_rows;i++){
        cout<<local_y[i]<<" ";
      }
      cout<<endl;
    }*/

}

void parallel_update_x(const int n, double* local_R, double* local_x, double* local_b, double* local_D, MPI_Comm comm){
  //process rank in the cartesian coordinates
  int grid_rank;
  //process coordinates in cartesia
  int grid_coord[2];
  //number of local_rows(local matrix)
  int n_local_rows = 0;
  //number of local_columns(local matrix)
  int n_local_cols = 0;
  //size of processors
  int p;
  //dimension
  int dims[2];
  //period
  int periods[2];

  MPI_Comm_rank(comm, &grid_rank);
  MPI_Comm_size(comm, &p);

  //find dims in the cartesian top
  MPI_Cart_get(comm, 2, dims, periods, grid_coord);

  //find the number of rows and cols owned by each processor
  n_local_rows = block_decompose_by_dim(n, comm, 0);
  n_local_cols = block_decompose_by_dim(n, comm, 1);
  //multiply R*x
  vector<double> temp_local_y(n_local_rows);
  distributed_matrix_vector_mult(n, local_R, local_x, &temp_local_y[0], comm);
  //first column will update the x
  if(grid_coord[1] == 0){
    for(int i = 0; i < n_local_rows; i++){
      local_x[i] = (local_b[i] - temp_local_y[i]) / local_D[i];
    }
  }
}

double parallel_L2_norm(const int n, double* local_A, double* local_x, double* local_b, MPI_Comm comm){
  //process rank in the cartesian coordinates
  int grid_rank;
  //process coordinates in cartesia
  int grid_coord[2];
  //number of local_rows(local matrix)
  int n_local_rows = 0;
  //number of local_columns(local matrix)
  int n_local_cols = 0;
  //size of processors
  int p;
  //dimension
  int dims[2];
  //period
  int periods[2];

  MPI_Comm_rank(comm, &grid_rank);
  MPI_Comm_size(comm, &p);

  //find dims in the cartesian top
  MPI_Cart_get(comm, 2, dims, periods, grid_coord);

  //find the number of rows and cols owned by each processor
  n_local_rows = block_decompose_by_dim(n, comm, 0);
  n_local_cols = block_decompose_by_dim(n, comm, 1);
  //multiply A*x
  vector<double> temp_local_y(n_local_rows);
  distributed_matrix_vector_mult(n, local_A, local_x, &temp_local_y[0], comm);

  //create a column communicator
  MPI_Comm col_comm;
  int remain_dims[2] = {true, false};
  MPI_Cart_sub(comm, remain_dims, &col_comm);

  //calculate (Ax-b)^2 in the first columns
  double l2_norm;
  if(grid_coord[1] == 0){
    double accum = 0.0;
    for(int i = 0; i < n_local_rows; i++){
      accum+=(temp_local_y[i]-local_b[i])*(temp_local_y[i]-local_b[i]);
    }
    MPI_Allreduce(&accum, &l2_norm, 1, MPI_DOUBLE, MPI_SUM, col_comm);
    l2_norm = sqrt(l2_norm);
  }
  //Bcast this value to the row so that every one has the value of l2_norm
  //row communicator
  MPI_Comm row_comm;
  remain_dims[0] = false;
  remain_dims[1] = true;
  MPI_Cart_sub(comm, remain_dims, &row_comm);

  //Bcasting root
  int bcast_root;
  int bcast_coords[ ] = {0};
  MPI_Cart_rank(row_comm, bcast_coords, &bcast_root);
  //bcast l2norm
  MPI_Bcast(&l2_norm, 1, MPI_DOUBLE, bcast_root, row_comm);
  return l2_norm;
}


// Solves Ax = b using the iterative jacobi method
void distributed_jacobi(const int n, double* local_A, double* local_b, double* local_x,
                MPI_Comm comm, int max_iter, double l2_termination)
{
    // TODO
    //process rank in the cartesian coordinates
    int grid_rank;
    //process coordinates in cartesia
    int grid_coord[2];
    //number of local_rows(local matrix)
    int n_local_rows = 0;
    //number of local_columns(local matrix)
    int n_local_cols = 0;
    //status for MPI_Recv
    MPI_Status status;
    //size of processors
    int p;
    //dimension
    int dims[2];
    //period
    int periods[2];

    MPI_Comm_rank(comm, &grid_rank);
    MPI_Comm_size(comm, &p);

    //find dims in the cartesian top
    MPI_Cart_get(comm, 2, dims, periods, grid_coord);

    //find the number of rows and cols owned by each processor
    n_local_rows = block_decompose_by_dim(n, comm, 0);
    n_local_cols = block_decompose_by_dim(n, comm, 1);

    //calculate local_D and local_R matrix
    vector<double>local_D(n_local_rows);
    vector<double>local_R(n_local_rows* n_local_cols);

    //only processor with coord[0]==coord[1] will own diagonal Elements
    if(grid_coord[0] == grid_coord[1]){
      /*memcpy(&local_R[0], local_A, n_local_cols*n_local_rows*sizeof(double));
      //now make diagonal elements zero in R and store diagonal elements in local vector
      for(int i = 0; i < n_local_rows;i++){
        local_D[i] = local_R[i*n_local_cols + i];
        local_R[i*n_local_cols + i] = 0;
      }*/
      //now make diagonal elements zero in R and store diagonal elements in local vector
      int j = 0;
      for(int i = 0 ;i <(n_local_cols*n_local_rows);i++){
        if(i % (n_local_cols + 1) == 0) {
          local_D[j] = local_A[i];
          local_R[i] = 0;
          j++;
        }
        else{
          local_R[i] = local_A[i];
        }
      }
    }
    //no diagnal elements so R = A
    else{
      memcpy(&local_R[0], local_A, n_local_cols*n_local_rows*sizeof(double));
    }
    //now collect the Diagonal elements into the first column using row communicators
    MPI_Comm row_comm;
    int remain_dims[2] = {false, true};
    MPI_Cart_sub(comm, remain_dims, &row_comm);

    //the diagonal processor will send the local_D to the first processor in the comm except(0,0) no need to send in that case
    if(grid_coord[0] == grid_coord[1] && grid_coord[0]!=0){
      int recv_rank;
      int recv_coord[] = {0};
      MPI_Cart_rank(row_comm, recv_coord, &recv_rank);
      //send local_D to the rank calculated above
      MPI_Send(&local_D[0], n_local_rows, MPI_DOUBLE, recv_rank, 111, row_comm);
    }
    //The first column need to recv the local_D(except(0,0 as it already has the data))
    if(grid_coord[1]==0 && grid_coord[0]!=0){
      int send_rank;
      int send_coord[] = {grid_coord[0]};
      MPI_Cart_rank(row_comm, send_coord, &send_rank);
      MPI_Recv(&local_D[0], n_local_rows, MPI_DOUBLE, send_rank, 111, row_comm, &status);
    }

    cout<<"hello"<<grid_coord[0]<<" "<<grid_coord[1]<<endl;
    for(int i = 0; i<n_local_rows;i++){
      cout<<local_D[i]<<endl;
    }

    //all the data is laid out, now starting Jacobi
    //initialize x to zero
    if(grid_coord[1] == 0){
      for(int i = 0;i < n_local_rows; i++){
        local_x[i] = 0;
      }
    }
    //iteration start
    int iter = 0;
    //check for condition
    while(parallel_L2_norm(n, local_A, local_x, local_b, comm) > l2_termination && iter < max_iter){
      parallel_update_x(n, &local_R[0], local_x, local_b, &local_D[0], comm);
      iter++;
    }
}

int main(int argc, char *argv[])
{
   // set up MPI
   MPI_Init(&argc, &argv);

   // get communicator size
   MPI_Comm comm = MPI_COMM_WORLD;
   int p;
   MPI_Comm_size(comm, &p);
   // get the dimensions
   int q = (int)sqrt(p);


   // create 2D cartesian grid for the processors (enable reordering)
   MPI_Comm grid_comm;
   int dims[2] = {q, q};
   int periods[2] = {0, 0};
   MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &grid_comm);
   // get the rank of process with coordinates (0,0)


     // simple 4 by 4 input matrix
  /*double A[5*5] = {10., -1., 2., 0.,10.,
                            -1., 11., -1., 3.,11.,
                            2., -1., 10., -1.,12.,
                            0.0, 3., -1., 8.,13.,
                            10.0, 13., -11., 18.,19.};

  double x[5] = {6., 25., -11., 15., 10.};

   int n = 5;*/
   /*double* local_A = NULL;
   double* local_x = NULL;
   distribute_matrix(n, &A[0], &local_A, grid_comm);
   distribute_vector(n, &x[0], &local_x, grid_comm);
   // Get the local vector for matrix multiplication
   //int n_local_cols = block_decompose_by_dim(n, grid_comm, 1);
   //double* local_xx = new double[n_local_cols];
   //transpose_bcast_vector(n, local_x, local_xx, grid_comm);

   // allocate local result space
   double* local_y = new double[block_decompose_by_dim(n, grid_comm, 0)];
   distributed_matrix_vector_mult(n, local_A, local_x, local_y, grid_comm);

   // gather results back to rank 0
   vector<double> y(n,0);
   gather_vector(n, local_y, &y[0], grid_comm);


   //if(local_x != NULL)cout<<local_x[0]<<endl;

   // finalize MPI*/

   // distribute the array onto local processors!
   // simple 4 by 4 input matrix
   double A[4*4] = {10., -1., 2., 0.,
                          -1., 11., -1., 3.,
                          2., -1., 10., -1.,
                          0.0, 3., -1., 8.};
   double b[4] =  {6., 25., -11., 15.};
   double x[4];
   int n=4;
   double* local_A = NULL;
   double* local_b = NULL;
   distribute_matrix(n, &A[0], &local_A, grid_comm);
//   distribute_vector(n, &b[0], &local_b, grid_comm);

   // allocate local result space
  // double* local_x = new double[block_decompose_by_dim(n, grid_comm, 0)];
   //distributed_jacobi(n, local_A, local_b, local_x, grid_comm, 100, 1e-10);

   // gather results back to rank 0
   //gather_vector(n, local_x, x, grid_comm);
   MPI_Finalize();
   return 0;
}
