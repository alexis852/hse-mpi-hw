#include <iostream>
#include <mpi.h>

using namespace std;

int main(int argc, char **argv) {
    int my_rank, total_ranks;

    MPI_Init(&argc, &argv);
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &total_ranks);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    cout << "Hello from rank " << my_rank << " out of " << total_ranks << endl;

    int num = 0;
    if (my_rank == 0) {
        num = 42;
        MPI_Send(&num, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
    } else if (my_rank == 1) {
        cout << "before " << num << endl;
        MPI_Recv(&num, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        cout << "after " << num << endl;
    }

    MPI_Finalize();
    return 0;
}
