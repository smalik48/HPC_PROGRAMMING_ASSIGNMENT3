/**
 * @file    mpi_jacobi.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements MPI functions for distributing vectors and matrixes,
 *          parallel distributed matrix-vector multiplication and Jacobi's
 *          method.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */

#include "mpi_jacobi.h"
#include "jacobi.h"
#include "utils.h"

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <vector>


using namespace std;

/*
 * TODO: Implement your solutions here
 */


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

      /*//print the vector;
      for(int i = 0;i<recvcnt;i++){
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

      /*if(root==grid_rank){
        for(int i=0;i<n;i++){
          cout<<output_vector[i]<<" ";
        }
        cout<<endl;
      }*/
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

    // Allocate spaces for local matrix(this way was nice else I will change the pointer)
    (*local_matrix) = new double[n_local_cols * n_local_rows];

    //get dimensions and grid_coord
    int rank00;
    int coords[2] = {0, 0};
    MPI_Cart_rank(comm, coords, &rank00);
    MPI_Cart_get(comm, 2, dims, periods, grid_coord);

    //iterate over grid row
    for(int i = 0; i < dims[0]; i++){
      grid_coord[0] = i;
      //go over each row within a block of matrix
      for(int j = 0; j < block_decompose(n, dims[0], i); j++){
        //go over each column in grid row
        for(int k = 0; k < dims[1]; k++){
          grid_coord[1] = k;
          //calculate the destination rank where the block is supposed to be delivered
          MPI_Cart_rank(comm, grid_coord, &dest_rank);

          //see if its master
          if(rank00 == grid_rank){
            if(i == 0 && k == 0) start_address = &input_matrix[ 0 + (j*n)];
            else if(i==0 && k!=0) start_address = &input_matrix[0 + (j*n) + (k * block_decompose(n, dims[1], k-1))];
            else start_address = &input_matrix[ (i*block_decompose(n, dims[0], i-1)*n)+(j*n)+(k * block_decompose(n, dims[1], k-1))];
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

    //result temp storing before reduction in each processor
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

    /*//print the vector
    if(grid_coord[1] == 0){
      for(int i=0;i<n_local_rows;i++){
        cout<<local_y[i]<<" ";
      }
      cout<<endl;
    }*/
}

// Solves Ax = b using the iterative jacobi method
void distributed_jacobi(const int n, double* local_A, double* local_b, double* local_x,
                MPI_Comm comm, int max_iter, double l2_termination)
{
    // TODO
}


// wraps the distributed matrix vector multiplication
void mpi_matrix_vector_mult(const int n, double* A,
                            double* x, double* y, MPI_Comm comm)
{
    // distribute the array onto local processors!
    double* local_A = NULL;
    double* local_x = NULL;
    distribute_matrix(n, &A[0], &local_A, comm);
    distribute_vector(n, &x[0], &local_x, comm);

    // allocate local result space
    double* local_y = new double[block_decompose_by_dim(n, comm, 0)];
    distributed_matrix_vector_mult(n, local_A, local_x, local_y, comm);

    // gather results back to rank 0
    gather_vector(n, local_y, y, comm);
}

// wraps the distributed jacobi function
void mpi_jacobi(const int n, double* A, double* b, double* x, MPI_Comm comm,
                int max_iter, double l2_termination)
{
    // distribute the array onto local processors!
    double* local_A = NULL;
    double* local_b = NULL;
    distribute_matrix(n, &A[0], &local_A, comm);
    distribute_vector(n, &b[0], &local_b, comm);

    // allocate local result space
    double* local_x = new double[block_decompose_by_dim(n, comm, 0)];
    distributed_jacobi(n, local_A, local_b, local_x, comm, max_iter, l2_termination);

    // gather results back to rank 0
    gather_vector(n, local_x, x, comm);
}
