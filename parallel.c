////////////////////////////////////////////////////////////////////////
// Time-Dependent 2D Heat Equation Solution (Runge-Kutta method)
// ECE 5720 - Cornell University
// Author: Alberto Lopez Delgado
//
// Instructions 
// 1. Compile using the command below:
// gcc -fopenmp parallel.c -o parallel -lm
// 2. Run using the command below:
// ./parallel Nx Ny num_threads
// 3. Follow instructions
// Output data will be placed in .txt files
////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <omp.h>
#include "sequential.h"

#define BILLION 1000000000.0

int main(int argc, char **argv)
{
	////////////////////////////////////////////////////////////////////
	// Material selection
	////////////////////////////////////////////////////////////////////
	
	int material;
	float conductivity;	// thermal conductivity (j/m*C*sec)
	float specific_heat;	// specific heat (j/kg*C)
	float density;			// density (kg/m^3)
	
	printf("This program computes the heat distribution in a rectangular plate over time.\n\n");
	printf("Available materials for the plate are:\n");
	printf("1 - Aluminum\n2 - Copper\n3 - Silver\n4 - Custom Material\n\n");
	printf("Please enter the number of your material of choice: ");
	scanf("%d", &material);
	
	// validate input
	while (material<1 || material>4){
		printf("Please enter a valid selection(1-4): ");
		scanf("%d", &material);
	}
		
	switch (material){
		case 1: // Aluminum
			conductivity = 204.3;
			specific_heat = 910;
			density = 2700;
			break;
		case 2: // Copper
			conductivity = 401;
			specific_heat = 390;
			density = 8940;
			break;
		case 3: // Silver
			conductivity = 629;
			specific_heat = 233;
			density = 10490;
			break;
		case 4: // Custom Material
			printf("Please enter the conductivity of your material: ");
			scanf("%f", &conductivity);
			printf("Please enter the specific heat of your material: ");
			scanf("%f", &specific_heat);
			printf("Please enter the density of your material: ");
			scanf("%f", &density);
			break;
		default:
			printf("Fatal Error: Invalid selection.\n");
	}
	
	////////////////////////////////////////////////////////////////////
	// Variable declarations
	////////////////////////////////////////////////////////////////////
	
	float Lx = 1;			// plate width (m)
	float Ly = 1;			// plate length (m)
	int Nx = atoi(argv[1]);	// number of nodes in x direction
	int Ny = atoi(argv[2]);	// number of nodes in y direction
	float T_initial = 0;	// initial temperature in all the nodes
	float T_up = 50;		// temperature at boundaries (Celsius)
	float T_down = 100;
	float T_left = 300;
	float T_right = 150; 
	float T_max = 300;
	float T_min = 0;
	float dt = 0.6;			// time step (1 sec)
	float tol = 0.5;		// tolerance for numerical sim (Celsius)
	float tol_ss = 0.001;	// tolerance for steady state (Celsius)
	int i, j, k;			// iteration counters
	float err_ss_max = 1;	// initial steady state max error
	float err_ss_min = 1;	// initial steady state min error
	float dx = Lx/Nx;		// delta x
	float dy = Ly/Ny;		// delta y
	float alpha = conductivity/(specific_heat*density);
	
	// variables for OpenMP
	int num_threads = atoi(argv[3]);
	
	// timekeeping variables
	struct timespec start, end;
	double ss_time, fd_time, sec, nsec;
	
	// test stability condition
	if (dt > 1/(2*alpha*((1/(dx*dx))+(1/(dy*dy))))){
		printf("Stability condition not met, choose a smaller dt.\n");
		return 1;
	}
	
	////////////////////////////////////////////////////////////////////
	// Computation of Steady State
	////////////////////////////////////////////////////////////////////	
	
	float *Tss[Nx+2], *Tss2[Nx+2], *diff[Nx+2];
	malloc_2d_mtx(Tss, Nx+2, Ny+2);	 // add extra rows for boundary temp
	malloc_2d_mtx(Tss2, Nx+2, Ny+2); // add extra rows for boundary temp	
	malloc_2d_mtx(diff, Nx+2, Ny+2); // add extra rows for boundary temp	
	init_2d_mtx(Tss, Nx+2, Ny+2);
	init_2d_mtx(Tss2, Nx+2, Ny+2);
	init_2d_mtx(diff, Nx+2, Ny+2);
	
	// initialize boundaries of temperature matrix
	init_col_temp_ss(Tss, Ny+2, 0, T_down);
	init_col_temp_ss(Tss2, Ny+2, 0, T_down);
	init_col_temp_ss(Tss, Ny+2, Ny, T_up);
	init_col_temp_ss(Tss2, Ny+2, Ny, T_up);
	init_col_temp_ss(Tss, Ny+2, Ny+1, T_up);	// used for graphing
	init_col_temp_ss(Tss2, Ny+2, Ny+1, T_up);	// used for graphing
	init_row_temp_ss(Tss, Nx+2, 0, T_left);
	init_row_temp_ss(Tss2, Nx+2, 0, T_left);
	init_row_temp_ss(Tss, Nx+2, Nx, T_right);
	init_row_temp_ss(Tss2, Nx+2, Nx, T_right);
	init_row_temp_ss(Tss, Nx+2, Nx+1, T_right);	// used for graphing
	init_row_temp_ss(Tss2, Nx+2, Nx+1, T_right);// used for graphing

	// compute the steady state
	printf("Computing Steady State...\n\n");
	float *temp_diff[Nx+2];
	malloc_2d_mtx(temp_diff, Nx+2, Ny+2);
	
	clock_gettime(CLOCK_MONOTONIC, &start); // start counting time
	
	omp_set_num_threads(num_threads);
	while (err_ss_max>=tol_ss || err_ss_min>=tol_ss){
		#pragma omp parallel
		{
			int nthreads = omp_get_num_threads();
			int ID = omp_get_thread_num();
			int ii, jj;
			
			for (ii=ID+1; ii<Nx; ii+=nthreads){ // ID+1 bc 1st row not touched
				for (jj=1; jj<Ny; jj++){
					Tss2[ii][jj] = 0.25*(Tss[ii+1][jj]+Tss[ii][jj+1]+Tss[ii-1][jj]+Tss[ii][jj-1]);
				}
			}
			
			// get Tss2(ii,jj)-Tss(ii,jj)
			for (ii=ID+1; ii<Nx+2; ii+=nthreads){
				for (jj=1; jj<Ny+2; jj++){
					temp_diff[ii][jj] = Tss2[ii][jj]-Tss[ii][jj];
				}
			}	
			
			err_ss_max = fabs(find_max(temp_diff, Nx+2, Ny+2));
			err_ss_min = fabs(find_min(temp_diff, Nx+2, Ny+2));
						
			for (ii=ID+1; ii<Nx+2; ii+=nthreads){
				for (jj=1; jj<Ny+2; jj++){
					Tss[ii][jj] = Tss2[ii][jj];
				}
			}
		}
	}
	
	clock_gettime(CLOCK_MONOTONIC, &end);
	sec = end.tv_sec-start.tv_sec; 
	nsec = (end.tv_nsec - start.tv_nsec)/BILLION; // converted to s
	ss_time = sec + nsec;
		
	printf("...finished Steady State.\n\n");
	
	////////////////////////////////////////////////////////////////////
	// Computation of Finite Difference Solution 
	////////////////////////////////////////////////////////////////////
	printf("Computing Finite Difference Solution...\n\n");
	float k1, k2;
	int num_steps = 75000;			// number of time steps for array
	float err_R_k_max[num_steps];		// initial Runge-Kutta max error
	float err_R_k_min[num_steps];		// initial Runge-Kutta min error
	k = 0;							// used to count while loop cycles
	err_R_k_max[k] = 100;
	err_R_k_min[k] = 100;
	
	float **T[num_steps];
	malloc_3d_mtx(T, Nx+2, Ny+2, num_steps);
	init_3d_mtx(T, Nx+2, Ny+2, num_steps);
	
	float *Tk[Nx+2];
	malloc_2d_mtx(Tk, Nx+2, Ny+2);
	
	// initialize boundaries of temperature matrix
	init_col_temp_fd(T, Ny+2, 0, num_steps, T_down);
	init_col_temp_fd(T, Ny+2, Ny, num_steps, T_up);
	init_col_temp_fd(T, Ny+2, Ny+1, num_steps, T_up);
	init_row_temp_fd(T, Nx+2, 0, num_steps, T_left);
	init_row_temp_fd(T, Nx+2, Ny, num_steps, T_right);
	init_row_temp_fd(T, Nx+2, Ny+1, num_steps, T_right);
	init_timestep_fd(T, 0, Nx+2, Ny+2, T_initial);
	
	clock_gettime(CLOCK_MONOTONIC, &start); // start time	
	
	while (err_R_k_max[k]>=tol || err_R_k_min[k]>=tol){
		#pragma omp parallel
		{
			int nthreads = omp_get_num_threads();
			int ID = omp_get_thread_num();
			int ii, jj;
			
			for (ii=ID+1; ii<Nx; ii+=nthreads){
				for (jj=1; jj<Ny; jj++){
					k1 = alpha*((T[ii-1][jj][k]-2*T[ii][jj][k]+T[ii+1][jj][k])/(dx*dx)+
								(T[ii][jj-1][k]-2*T[ii][jj][k]+T[ii][jj+1][k])/(dy*dy));
					update_Tk(T, Nx+2, Ny+2, Tk, k, k1, dt);
					k2 = alpha*((Tk[ii-1][jj]-2*Tk[ii][jj]+Tk[ii+1][jj])/(dx*dx)+
								(Tk[ii][jj-1]-2*Tk[ii][jj]+Tk[ii][jj+1])/(dy*dy));
					T[ii][jj][k+1] = T[ii][jj][k] + (dt/2)*(k1+k2);
				}
			}
		}

		k++;
			
		#pragma omp parallel
		{
			int nthreads = omp_get_num_threads();
			int ID = omp_get_thread_num();
			int ii, jj;
			
			// get T(:,:,k)-Tss
			for (ii=ID+1; ii<Nx+2; ii+=nthreads){
				for (jj=1; jj<Ny+2; jj++){
					temp_diff[ii][jj] = T[ii][jj][k]-Tss[ii][jj];
				}
			}	
		}
				
		err_R_k_max[k] = fabs(find_max(temp_diff, Nx+2, Ny+2));
		err_R_k_min[k] = fabs(find_min(temp_diff, Nx+2, Ny+2));
		
		// Test for convergence
		if ((int)(err_R_k_max[k]*10000) == (int)(err_R_k_max[k-1]*10000) && err_R_k_max[k]!=0){
			printf("Max Error is not converging.\n");
		}
			
		if ((int)(err_R_k_min[k]*10000) == (int)(err_R_k_min[k-1]*10000) && err_R_k_min[k]!=0){
			printf("Min Error is not converging.\n");
		}	
	}
	
	clock_gettime(CLOCK_MONOTONIC, &end); // finish time
	sec = end.tv_sec-start.tv_sec; 
	nsec = (end.tv_nsec - start.tv_nsec)/BILLION; // converted to s
	fd_time = sec + nsec;
	
	float SStime = k*dt;	// used for graph
	
	printf("...finite Difference calculation finished\n\n");
	
	////////////////////////////////////////////////////////////////////
	// Output files creation
	////////////////////////////////////////////////////////////////////
	
	printf("Creating output files...\n\n");
	
	// create files
	FILE *tss;
	FILE *t;
	FILE *int_vars;
	FILE *f_vars;
		
	
	// print Tss matrix
	tss = fopen("tss.txt", "w+");

	if (tss == NULL){
		printf("Error: File not created.\n");
		return 1;
	}
	
	for (i=0; i<Nx+2; i++){
		for (j=0; j<Ny+2; j++){
			fprintf(tss, "%f\t", Tss[i][j]);
		}
		fprintf(tss, "\n");
	}
	fclose(tss);
	
	// print T matrix 
	int step;
	
	t = fopen("t.txt", "w+");
	
	for (step=0; step<k; step++){
		for (i=0; i<Nx+2; i++){
			for (j=0; j<Ny+2; j++){
				fprintf(t, "%f\t", T[i][j][step]);
			}
			fprintf(t, "\n");
		}
	}
	
	// print integer variables associated with Tss and T
	int_vars = fopen("int_vars.txt", "w+");
	
	if (int_vars == NULL){
		printf("Error: File not created.\n");
		return 1;
	}
	
	fprintf(int_vars, "%d\n", Nx);
	fprintf(int_vars, "%d\n", Ny);
	fprintf(int_vars, "%d\n", k);
	
	fclose(int_vars);
	
	// print float variables associated with Tss and T
	f_vars = fopen("f_vars.txt", "w+");
	
	if (f_vars == NULL){
		printf("Error: File not created.\n");
		return 1;
	}
	
	fprintf(f_vars, "%f\n", dx);
	fprintf(f_vars, "%f\n", dy);
	fprintf(f_vars, "%f\n", SStime);
	fprintf(f_vars, "%f\n", Lx);
	fprintf(f_vars, "%f\n", Ly);
	fprintf(f_vars, "%f\n", T_max);	// max temp in boundaries
	fprintf(f_vars, "%f\n", T_min);	// min temp in boundaries
	
	printf("...output files successfully created. Program finished.\n\n");
	
	printf("Execution times (seconds):\n");
	printf("Steady State\t\tFinite Difference\n");
	printf("%.9f\t\t%.9f\n", ss_time, fd_time);
	
	return 0;
}
