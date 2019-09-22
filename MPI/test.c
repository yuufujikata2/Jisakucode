#include<mpi.h>
#include<stdio.h>
#include<unistd.h>

void print(int numprocs, int myid);


int main(int argc, char *argv[]){

  int myid, numprocs;
  
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  if (myid==0)
    print(numprocs,myid);
  
  MPI_Finalize();
  
  return 0;
}

void print(numprocs,myid){

  int n, i, l;
  int t=0;

  l=100;

  if (myid == 0 ){
    sleep(3);
    printf("a\n");
    t=1;
  }
//  MPI_Barrier(MPI_COMM_WORLD);
//  MPI_Bcast(&t, 1, MPI_INT, 0, MPI_COMM_WORLD);

/*  while (t == 0){

  }
*/ 
  printf("myid= %d t= %d\n",myid,t);   
  
  if(myid < numprocs -1){
    for (i = l / numprocs * myid ; i < l / numprocs * (myid + 1) ; i++){
      printf ("myid = %d i = %d\n", myid, i);
    }
  }
  else{
    for (i = l / numprocs * myid ; i < l ; i++){
      printf ("myid = %d i = %d\n", myid, i);
    }
  }
  //printf ("Hello World %d\n",myid);

}

