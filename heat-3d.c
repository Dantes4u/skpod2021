#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include "mpi.h"

// Сохранение контрольной точки
static void save_control_point(int rank, int n, int in_use, double A[in_use + 2][n][n], double B[in_use + 2][n][n], int last_step, int last_block) {
    char filename[42];
    snprintf(filename, 42, "%d_%d_%d.txt", rank, last_step, last_block);
    FILE *point = fopen(filename, "w");
    // Записываем массив A
    for (int i = 0; i < in_use + 2; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                fprintf(point, "%0.2lf ", A[i][j][k]);
            }
            fprintf(point, "\n");
        }
    }

    // Записываем массив B
    for (int i = 0; i < in_use + 2; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                fprintf(point, "%0.2lf ", B[i][j][k]);
            }
            fprintf(point, "\n");
        }
    }
    fclose(point);
}


// Иницилизация массивов
static void init_array(int rank, int processes, int n, int in_use, double A[in_use + 2][n][n], double B[in_use + 2][n][n]) {
    int start_use = (rank - 1) * n / (processes - 1);
    for (int i = 1; i <= in_use; i++) {
        int i_global = start_use + i - 1;
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                A[i][j][k] = B[i][j][k] = (double) (i_global + j + (n - k)) * 10 / n;
            }
        }
    }
}

int rank_to_kill;
int block_to_kill;
int step_to_kill;
//Убийство процесса
static void kill_process(int rank, int curr_block, int curr_t_step) {
    if ((rank == rank_to_kill) && (curr_block == block_to_kill) && (curr_t_step == step_to_kill)){
		raise(SIGKILL);
    }
}

// Вычисления A и B
static void kernel_heat_3d(int rank, int processes, int steps, int n, int in_use, double A[in_use + 2][n][n], double B[in_use + 2][n][n], int start_step, int start_block) {
    int killed_rank = -1, t = start_step, rc;
    MPI_Request sender_first, sender_second, rec;
    MPI_Status getter_first, getter_second;
    switch (start_block + 1) {
        case 0:
        case 1:
            goto first_block;
        case 2:
            goto second_block;
        case 3:
            goto third_block;
        case 4:
            goto fourth_block;
        case 5:
            goto fifth_block;
        case 6:
            goto sixth_block;
        default:
            goto seventh_block;
    }
    for (; start_step <= steps; start_step++) {

first_block:
        kill_process(rank, 1, start_step);
        if (rank != 1) {
            MPI_Isend(&A[1][0][0], n * n, MPI_DOUBLE, (rank - 1 != killed_rank) ? (rank - 1) : processes, 112, MPI_COMM_WORLD, &sender_first);
            rc = MPI_Recv(&A[0][0][0], n * n, MPI_DOUBLE, (rank - 1 != killed_rank) ? (rank - 1) : processes, 112, MPI_COMM_WORLD, &getter_first);
            if (rc != 0) {
                killed_rank = rank - 1;
                MPI_Isend(&killed_rank, 1, MPI_INT, processes, 400, MPI_COMM_WORLD, &rec);
                goto first_block;
            }
        }
        save_control_point(rank, n, in_use, A, B, start_step, 1);

second_block:
        kill_process(rank, 2, start_step);
        if (rank != processes - 1) {
            MPI_Isend(&A[in_use][0][0], n * n, MPI_DOUBLE, (rank + 1 != killed_rank) ? (rank + 1) : processes, 112, MPI_COMM_WORLD, &sender_second);
            rc = MPI_Recv(&A[in_use + 1][0][0], n * n, MPI_DOUBLE, (rank + 1 != killed_rank) ? (rank + 1) : processes, 112,
                     MPI_COMM_WORLD, &getter_second);
            if (rc != 0) {
                killed_rank = rank + 1;
                MPI_Isend(&killed_rank, 1, MPI_INT, processes, 400, MPI_COMM_WORLD, &rec);
                goto second_block;
            }
        }
        save_control_point(rank, n, in_use, A, B, start_step, 2);

third_block:
        kill_process(rank, 3, start_step);
        for (int i = 1; i <= in_use; i++) {
            if ((i == 1 && rank == 1) || (i == in_use && rank == processes - 1)) {
                continue;
            }
            for (int j = 1; j < n - 1; j++) {
                for (int k = 1; k < n - 1; k++) {
                    double residue;
                    B[i][j][k] = A[i][j][k] * 0.25;
                    residue = A[i + 1][j][k] + A[i - 1][j][k]
                              + A[i][j + 1][k] + A[i][j - 1][k]
                              + A[i][j][k + 1] + A[i][j][k - 1];
                    residue *= 0.125;
                    B[i][j][k] += residue + 1.0;
                }
            }
        }
        save_control_point(rank, n, in_use, A, B, start_step, 3);


fourth_block:
        kill_process(rank, 4, start_step);
        if (rank != 1) {
            MPI_Isend(&B[1][0][0], n * n, MPI_DOUBLE, (rank - 1 != killed_rank) ? (rank - 1) : processes, 112, MPI_COMM_WORLD, &sender_first);
            rc = MPI_Recv(&B[0][0][0], n * n, MPI_DOUBLE, (rank - 1 != killed_rank) ? (rank - 1) : processes, 112, MPI_COMM_WORLD, &getter_first);
            if (rc != 0) {
                killed_rank = rank - 1;
                MPI_Isend(&killed_rank, 1, MPI_INT, processes, 400, MPI_COMM_WORLD, &rec);
                goto fourth_block;
            }
        }
        save_control_point(rank, n, in_use, A, B, start_step, 4);

fifth_block:
        kill_process(rank, 5, start_step);
        if (rank != processes - 1) {
            MPI_Isend(&B[in_use][0][0], n * n, MPI_DOUBLE, (rank + 1 != killed_rank) ? (rank + 1) : processes, 112, MPI_COMM_WORLD, &sender_second);
            rc = MPI_Recv(&B[in_use + 1][0][0], n * n, MPI_DOUBLE, (rank + 1 != killed_rank) ? (rank + 1) : processes, 112, MPI_COMM_WORLD, &getter_second);
            if (rc != 0) {
                killed_rank = rank + 1;
                MPI_Isend(&killed_rank, 1, MPI_INT, processes, 400, MPI_COMM_WORLD, &rec);
                goto fifth_block;
            }
        }
        save_control_point(rank, n, in_use, A, B, start_step, 5);

sixth_block:
        kill_process(rank, 6, start_step);
        for (int i = 1; i <= in_use; i++) {
            if ((i == 1 && rank == 1) || (i == in_use && rank == processes - 1)) {
                continue;
            }
            for (int j = 1; j < n - 1; j++) {
                for (int k = 1; k < n - 1; k++) {
                    double residue;
                    A[i][j][k] = B[i][j][k] * 0.25;
                    residue = B[i + 1][j][k] + B[i - 1][j][k]
                              + B[i][j + 1][k] + B[i][j - 1][k]
                              + B[i][j][k + 1] + B[i][j][k - 1];
                    residue *= 0.125;
                    A[i][j][k] += residue + 2.0;
                }
            }
        }
        save_control_point(rank, n, in_use, A, B, start_step, 6);

seventh_block:
        kill_process(rank, 7, start_step);
    }

}

// Сбор массивов воедино
static void all_mas_to_one(int processes, int n, int in_use, double A_all[n][n][n], double B_all[n][n][n]) {
    MPI_Request rec;
    int killed_rank = -1, rc;
    for (int i = 1; i < processes; i++) {
        int start_use = ((i - 1) * n) / (processes - 1);
problem_A:
        rc = MPI_Recv(&A_all[start_use][0][0], in_use * n * n, MPI_DOUBLE, (i != killed_rank) ? i : processes, 999, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (rc != 0) {
            killed_rank = i;
            MPI_Isend(&killed_rank, 1, MPI_INT, processes, 400, MPI_COMM_WORLD, &rec);
            goto problem_A;
        }
    }
    for (int i = 1; i < processes; i++) {
        int start_use = ((i - 1) * n) / (processes - 1);
problem_B:
        rc = MPI_Recv(&B_all[start_use][0][0], in_use * n * n, MPI_DOUBLE, (i != killed_rank) ? i : processes, 999, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (rc != 0) {
            killed_rank = i;
            MPI_Isend(&killed_rank, 1, MPI_INT, processes, 400, MPI_COMM_WORLD, &rec);
            goto problem_B;
        }
    }
}


// Печать массива в файл
static void print_array(int n, double A[n][n][n], FILE *f) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                fprintf(f, "%0.2lf ", A[i][j][k]);
            }
            fprintf(f, "\n");
        }
    }
}

// Рабочие процессы отправляют результаты процессу-координатору
static void send_results(int rank, int n, int in_use, double A[in_use + 2][n][n], double B[in_use + 2][n][n]) {
    MPI_Send(&A[1][0][0], in_use * n * n, MPI_DOUBLE, 0, 999, MPI_COMM_WORLD);
    MPI_Send(&B[1][0][0], in_use * n * n, MPI_DOUBLE, 0, 999, MPI_COMM_WORLD);
}

// Восстановление из последней контрольной точки
static void load_control_point(int rank, int n, int in_use, double A[in_use + 2][n][n], double B[in_use + 2][n][n], int *last_step, int *last_block) {
    char filename[42];
    int max_t_step = *last_step, max_block = *last_block;
    int flag = 0;
    printf("Загружаем контрольную точку из файла: %d_%d_%d.txt\n", rank, *last_step, *last_block);
    snprintf(filename, 42, "%d_%d_%d.txt", rank, *last_step, *last_block);
    FILE *point = fopen(filename, "r");
    // Читаем массив A
    for (int i = 0; i < in_use + 2; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                fscanf(point, "%lf", &A[i][j][k]);
            }
        }
    }
    // Читаем массив B
    for (int i = 0; i < in_use + 2; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                fscanf(point, "%lf", &B[i][j][k]);
            }
        }
    }
    fclose(point);
    printf("Контрольная точка загружена!\n");
}

int main(int argc, char** argv) {

    int num_processes, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);


    int n = 8;
    int steps = 3;

    if (argc == 3) {
        steps = atoi(argv[1]);
        n = atoi(argv[2]);
    }
    int processes = num_processes - 1;
    rank_to_kill = 1;
    block_to_kill = 4;
    step_to_kill = 2;



    int in_use = n / (processes - 1);

    if (rank == 0) {
    	printf("\n\tSTART\n");
		printf("Число шагов = %d\n", steps);
		printf("Размер измерения 3х-мерной кубической матрицы = %d\n", n);
		printf("Число процессов - num_processes = %d\n", num_processes);
		double (*A_all)[n][n][n];
        double (*B_all)[n][n][n];

        A_all = (double (*)[n][n][n]) malloc((n) * (n) * (n) * sizeof(double));
        B_all = (double (*)[n][n][n]) malloc((n) * (n) * (n) * sizeof(double));

        all_mas_to_one(processes, n, in_use, *A_all, *B_all);

        // Записываем массивы в файл "matrix.txt"
        FILE *f;
        f = fopen("matrix.txt", "w");
        print_array(n, *A_all, f);
        print_array(n, *B_all, f);
        fclose(f);

        free((void *) A_all);
        free((void *) B_all);

        // Отправляем резервному процессу, что никто не умер
        int success = 0;
        MPI_Request rec;
        MPI_Isend(&success, 1, MPI_INT, processes, 400, MPI_COMM_WORLD, &rec);

    } else {

        double (*A)[in_use + 2][n][n];
        double (*B)[in_use + 2][n][n];

        A = (double (*)[in_use + 2][n][n]) malloc((in_use + 2) * (n) * (n) * sizeof(double));
        B = (double (*)[in_use + 2][n][n]) malloc((in_use + 2) * (n) * (n) * sizeof(double));

        int first_step = 1, start_block = 0;

        if (rank == processes) {

            // Ищем умерший процесс
            MPI_Status status_test;
            MPI_Request request_test[processes];
            int test[processes];
            for (int i = 0; i < processes; i++) {
                MPI_Irecv(&test[i], 1, MPI_INT, i, 404, MPI_COMM_WORLD, &request_test[i]);
            }
            int fail = -1;
            MPI_Waitany(processes, request_test, &fail, &status_test);

            printf("Процесс %d умер в блоке %d на шаге %d\n", fail, block_to_kill, step_to_kill);

            if (fail == 0) {
                goto end_block;
            }
            if (fail == 1) {
                fail = fail + 1;
            } else {
                fail = fail - 1;
            }
            MPI_Request request_rank[2];
            MPI_Status status;
            int ranks[2];
            MPI_Irecv(&ranks[0], 1, MPI_INT, 0, 400, MPI_COMM_WORLD, &request_rank[0]);
            MPI_Irecv(&ranks[1], 1, MPI_INT, fail, 400, MPI_COMM_WORLD, &request_rank[1]);
            int idx_success = -1;
            MPI_Waitany(2, request_rank, &idx_success, &status);
            int new_rank = ranks[idx_success];
            printf("Резервный процесс заменил процесс №%d!\n", new_rank);
            rank_to_kill = -1;
            block_to_kill = -1;
            step_to_kill = -1;
            first_step = 2;
            start_block = 4;
            //load_control_point(new_rank, n, in_use, *A, *B, 2, 4);
            load_control_point(new_rank, n, in_use, *A, *B, &first_step, &start_block);
            rank = new_rank;
            switch (start_block + 1) {
                case 0:
                    goto zero_block;
                default:
                    goto execute_block;
            }
        }
zero_block:
        kill_process(rank, 0, 0);
        init_array(rank, processes, n, in_use, *A, *B);
        save_control_point(rank, n, in_use, *A, *B, 1, 0);
execute_block:
        kernel_heat_3d(rank, processes, steps, n, in_use, *A, *B, first_step, start_block);
		send_results(rank, n, in_use, *A, *B);

end_block:
        free((void *) A);
        free((void *) B);
    }
    if (rank == 0) {
		printf("\n\tEND\n");
    }
    MPI_Finalize();
}

