#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <mpi.h>

using namespace std;

const int l = 1;
const int k = 1;

#define PRINT 1
#define PRINT_VAR(x) cout << #x" = " << x << endl

int main(int argc, char **argv) {
    if (argc != 4) {
        cerr << "Wrong number of args." << endl;
        return 1;
    }
    int n_proc = atoi(argv[1]);
    int N = atoi(argv[2]);
    double T = atof(argv[3]);

    double h = (double)l / (N - 1);
    double tau = 0.5 * h * h / k;
    double m = T / tau;
    double m_intpart;
    assert(modf(m, &m_intpart) == 0.0 && "M is not integer.");
    int M = (int)m_intpart + 1;

    int my_rank, size;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    assert(size == n_proc && "size != n_proc for whatever reason.");

    int N_local = N / n_proc + (my_rank < N % n_proc);
    double u_0, start, finish;
    int i, j;
    vector<vector<double>> u(M, vector<double>(N_local));

    if (my_rank == 0) {
        // PRINT_VAR(h);
        // PRINT_VAR(tau);
        // PRINT_VAR(M);

        start = MPI_Wtime();

        // for (int i = 0; i < n_proc; ++i) {
        //     left[i] = (i == 0) ? 1 : right[i - 1];
        //     right[i] = left[i] + N / n_proc;
        //     if (i < N % n_proc) {
        //         right[i] += 1;
        //     }
        //     cout << i << ": " << left[i] << " " << right[i] << endl;
        // }

        u_0 = 1;
        for (int rank = 1; rank < size; ++rank) {
            MPI_Send(&u_0, 1, MPI_DOUBLE, rank, 0, MPI_COMM_WORLD);
        }
        // MPI_Bcast(&buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
    } else {
        MPI_Recv(&u_0, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
    }

    for (j = 0; j < N_local; ++j) {
        u[0][j] = u_0;
    }

    // cout << "OK " << my_rank << endl;

    double u_left, u_right;
    for (i = 1; i < M; ++i) {
        j = 0;
        u_left = 0;
        if (i != 1 && my_rank != 0) {
            MPI_Recv(&u_left, 1, MPI_DOUBLE, my_rank - 1, 1, MPI_COMM_WORLD, &status);
        } else if (i == 1 && my_rank != 0) {
            u_left = u_0;
        }
        u[i][j] = u[i - 1][j] + 0.5 * (u[i - 1][j + 1] - 2 * u[i - 1][j] + u_left);
        // cout << "OK " << i << " 0 " << my_rank << endl;

        if (i != M - 1 && my_rank != 0) {
            MPI_Send(&u[i][j], 1, MPI_DOUBLE, my_rank - 1, 2, MPI_COMM_WORLD);
        }
        // cout << "OK " << i << " 1 " << my_rank << endl;

        for (j = 1; j < N_local - 1; ++j) {
            u[i][j] = u[i - 1][j] + 0.5 * (u[i - 1][j + 1] - 2 * u[i - 1][j] + u[i - 1][j - 1]);
        }
        // cout << "OK " << i << " 2 " << my_rank << endl;

        u_right = 0;
        if (i != 1 && my_rank != n_proc - 1) {
            MPI_Recv(&u_right, 1, MPI_DOUBLE, my_rank + 1, 2, MPI_COMM_WORLD, &status);
        } else if (i == 1 && my_rank != n_proc - 1) {
            u_right = u_0;
        }
        // cout << "OK " << i << " 3 " << my_rank << endl;

        u[i][j] = u[i - 1][j] + 0.5 * (u_right - 2 * u[i - 1][j] + u[i - 1][j - 1]);
        if (i != M - 1 && my_rank != n_proc - 1) {
            MPI_Send(&u[i][j], 1, MPI_DOUBLE, my_rank + 1, 1, MPI_COMM_WORLD);
        }
        // cout << "OK " << i << " 4 " << my_rank << endl;
    }

    if (my_rank == 0) {
        finish = MPI_Wtime();
        cout << finish - start << endl;
        if (PRINT) {
            int q = 0;
            for (int rank = 0; rank < size; ++rank) {
                for (j = 0; j < N / n_proc + (rank < N % n_proc); ++j) {
                    if (rank != 0) {
                        MPI_Recv(&u[M - 1][j], 1, MPI_DOUBLE, rank, j, MPI_COMM_WORLD, &status);
                    }
                    if (q % 5 == 0) {
                        cout << fixed << setprecision(5)
                             << "u(" << q * h << ", " << T << ") = "
                             << u[M - 1][j] << endl;
                    }
                    ++q;
                }
            }
        }
    } else {
        if (PRINT) {
            for (j = 0; j < N_local; ++j) {
                MPI_Send(&u[M - 1][j], 1, MPI_DOUBLE, 0, j, MPI_COMM_WORLD);
            }
        }
    }

    MPI_Finalize();
    return 0;
}
