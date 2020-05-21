void malloc_2d_mtx(float **src, int rows, int cols)
{
	int row;
	
	for (row=0; row<rows; row++){
		src[row] = (float *)malloc(cols * sizeof(float));
	}
}

void init_2d_mtx(float **src, int rows, int cols)
{
	int row, col;
	
	for (row=0; row<rows; row++){
		for (col=0; col<cols; col++){
			src[col][row] = 0;
		}
	}
}

void print_2d_mtx(float **mat, int rows, int cols)
{
	int row, col;
	
	if (rows<=8 && cols <= 8){
		for (col=0; col<cols*16; col++){
			printf("-");
		}
		
		printf("-\n");
		
		for (row=0; row<rows; row++){
			for (col=0; col<cols; col++){
				if (col==0){
					printf("| ");
				}
				
				printf("%f   \t|", mat[row][col]);
			}
			
			printf("\n");
			
			for (col=0; col<cols*16; col++){
				printf("-");
			}
			printf("-\n");
		}
	}
	else{
		printf("\nMatrix too large to print!\n");
	}
}

void init_col_temp_ss(float **src, int col_size, int col_idx, float temp)
{
	int row;
	
	for (row=0; row<col_size; row++){
		src[row][col_idx] = temp;
	}
}

void init_row_temp_ss(float **src, int row_size, int row_idx, float temp)
{
	int col;
	
	for (col=0; col<row_size; col++){
		src[row_idx][col] = temp;
	}
}

float find_max(float **src, int rows, int cols)
{
	int row, col;
	float max = src[0][0];
	
	for (row=0; row<rows; row++){
		for (col=0; col<cols; col++){
			if (src[row][col] > max) max = src[row][col];
		}
	}
	
	return max;
}

float find_min(float **src, int rows, int cols)
{
	int row, col;
	float min = src[0][0];
	
	for (row=0; row<rows; row++){
		for (col=0; col<cols; col++){
			if (src[row][col] < min) min = src[row][col];
		}
	}
	
	return min;
}

void malloc_3d_mtx(float ***src, int rows, int cols, int height)
{
	int row, col, step;
	
	for (row=0; row<rows; row++){
		src[row] = (float **) malloc(sizeof(float*)*cols);
		for (col=0; col<cols; col++){
			src[row][col] = (float *) malloc(sizeof(float)*height);
		}
	}
}

void init_3d_mtx(float ***src, int rows, int cols, int height)
{
	int row, col, step;
	
	for (row=0; row<rows; row++){
		for (col=0; col<cols; col++){
			for (step=0; step<height; step++){
				src[row][col][step] = 0;
			}
		}
	}
}

void init_col_temp_fd(float ***src, int col_size, int col_idx, int num_steps, float temp)
{
	int row, step;
	
	for (row=0; row<col_size; row++){
		for (step=0; step<num_steps; step++){
			src[row][col_idx][step] = temp;
		}
	}
}

void init_row_temp_fd(float ***src, int row_size, int row_idx, int num_steps, float temp)
{
	int row, col, step;
	
	for (col=0; col<row_size; col++){
		for (step=0; step<num_steps; step++){
			src[row_idx][col][step] = temp;
		}
	}
}

void init_timestep_fd(float ***src, int time_idx, int rows, int cols, float temp)
{
	int row, col, step;
	
	for (row=0; row<rows; row++){
		for (col=0; col<cols; col++){
			src[row][col][time_idx] = temp;
		}
	}
}

void print_3d_mtx(float ***src, int rows, int cols, int time)
{
	int row, col;
	
	if (rows<=8 && cols <= 8){		
		for (row=0; row<rows; row++){
			for (col=0; col<cols; col++){
				printf("%.2f    \t", src[row][col][time]);
			}
			printf("\n");
		}
	}
	else{
		printf("\nMatrices too large to print!\n");
	}
}

void update_Tk(float ***src, int rows, int cols, float **dest, int k, float k1, float dt)
{
	int row, col;
	
	for (row=0; row<rows; row++){
		for (col=0; col<cols; col++){
			dest[row][col] = src[row][col][k] + k1*dt;
		}
	}
}
