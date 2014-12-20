#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdlib.h>
#include <pthread.h>


#define RADIUS 1000

int debug = 0;
int serial = 0;
int openMp = 0;
int pthread = 0;

int size;
int thread_count = 8;

double drand();
void print_command_error_message();
int check_arguments(int argc, char** argv);
void debug_message(const char *format, ...);
void debug_matrix(double *factors, double* values);
int matrix_index(int row, int column);
void excute_elimination(const char* method_name, double *factors, double *values, void (*executor)(double*, double*));
int clock_diff_in_milies(clock_t start, clock_t end);

void serial_elimination(double* factors, double* values);
void open_mp_elimination(double* factors, double* values);
void pthread_elimination(double* factors, double* values);

void* pthread_run(void* parameter);


int main(int argc, char **argv)
{

	if(!check_arguments(argc, argv))
	{
		print_command_error_message();
		return 1;
	}


	double *factors = (double*) malloc(sizeof(double) * size*size);
	double *values = (double*) malloc(sizeof(double) * size);

	srand(time(NULL));

	//init factor matrix with random values
	int i,j;
	for(i=0; i<size; i++)
		for(j=0; j<size; j++)
			factors[matrix_index(i,j)] = drand();

	//init values vector with random values
	for(i=0; i<size; i++)
		values[i] = drand();

	debug_message("----------------------- init matrix -----------------------\n");
	debug_matrix(factors, values);
  	debug_message("-----------------------------------------------------------\n");
	

  	if (serial)
		excute_elimination("serial ", factors, values, &serial_elimination);

	if (openMp)
		excute_elimination("open mp", factors, values, &open_mp_elimination);

	if (pthread)
		excute_elimination("pthread", factors, values, &pthread_elimination);
  	
}


double drand()
{
	return (rand()/(RAND_MAX*1.0))*RADIUS;
}

void print_command_error_message()
{
	printf("%s\n", "error usage the command shuld be in form :\n command [options] [size] [thread_count] \n options : \n \t -d \t debug \n \t -s \t serial execution \n \t -o \t OpenMp \n \t -p \t Pthread");
}

int check_arguments(int argc, char** argv)
{
	if (argc != 4)
	{
		printf("invalid parameter size\n");
		return 0;
	}

	if(argv[1][0] != '-')
	{
		printf("wrong option format\n");
		return 0;
	}
	
	int i;
	for(i=1; argv[1][i] != 0; i++)
	{
		if (argv[1][i] == 'd')
			debug = 1;
		else if (argv[1][i] == 'o')
			openMp = 1;
		else if (argv[1][i] == 'p')
			pthread = 1;
		else if (argv[1][i] == 's')
			serial = 1;
		else
		{
			printf("unknown option : %c\n", argv[1][i]);
			return 0;
		}
	}

	size = atoi(argv[2]);
	if (size == 0)
	{
		printf("%s\n", "invalid size");
		return 0;
	}

	thread_count = atoi(argv[3]);
	if (thread_count == 0)
	{
		printf("%s\n", "invalid thread_count");
		return 0;
	}

	return 1;
}

void debug_message(const char *format, ...)
{
	if (debug)
	{
		va_list args;
	    va_start( args, format );
	    vprintf( format, args );
	    va_end( args );
	}
}

void debug_matrix(double *factors, double* values)
{
	if (debug)
	{
		int i,j;
		for(i=0; i<size; i++){
			for(j=0; j<size; j++)
				printf("%.2f ", factors[matrix_index(i,j)]);

			printf("| %.2f \n", values[i]);
	  	}
	}	
}

int matrix_index(int row, int column)
{
	return row * size + column;
}

void excute_elimination(const char* method_name, double *factors, double *values, void (*executor)(double*, double*))
{
	printf("=============================== start %s =================================\n", method_name);

	double *factors_copy = (double*) malloc(sizeof(double) * size*size);
	double *values_copy = (double*) malloc(sizeof(double) * size);

	memcpy(factors_copy, factors, size*size * sizeof(double));
	memcpy(values_copy, values, size * sizeof(double));

	debug_message("------------------------------- start matrix ---------------------------\n");
	debug_matrix(factors_copy, values_copy);
	debug_message("-------------------------------------------------------------------------\n");

	time_t start_time = time(NULL);

	executor(factors_copy, values_copy);

	time_t end_time = time(NULL);

	debug_message("------------------------------- result matrix ---------------------------\n");
	debug_matrix(factors_copy, values_copy);
	debug_message("-------------------------------------------------------------------------\n");

	free(factors_copy);
	free(values_copy);

	printf("total time : %d s\n", end_time - start_time);

	printf("================================= end %s =================================\n", method_name);
}

int clock_diff_in_milies(clock_t start, clock_t end)
{
	return (end - start)/(CLOCKS_PER_SEC/1000);
}

void serial_elimination(double *factors, double *values){
	int i,j,k;
	double ratio;
	for (i = 0; i < size-1; i++)
	{
	   for (j = i + 1; j < size; j++) 
	   {
	   		if (factors[matrix_index(i,i)] == 0)
	   		{
	   			printf("can not countinue cause (%d,%d) is ziro\n", i, i);
	   			return;
	   		}
	        ratio = factors[matrix_index(j,i)]/factors[matrix_index(i,i)];
	        if (debug)
	        	debug_message("ratio [(%d, %d) : (%f, %f)] -> %f\n", i, j, factors[matrix_index(j,i)], factors[matrix_index(i,i)], ratio);

	        for (k = i; k < size; k++) 
	        	factors[matrix_index(j,k)] -= (ratio*factors[matrix_index(i,k)]);
	        values[j] -= (ratio*values[i]);
	   }
  	}	
}

void open_mp_elimination(double* factors, double* values)
{
	int i,j,k;
	double ratio;
	{
		for (i = 0; i < size-1; i++) 
		{
		   #pragma omp parallel for private(j, k, ratio) shared(i, factors, values) num_threads(thread_count)
	       for (j = i + 1; j < size; j++) 
	       {
	       		if (factors[matrix_index(i,i)] == 0)
	       		{
	       			printf("can not countinue cause (%d,%d) is ziro\n", i, i);
	       		}
	       		else
	       		{
		            ratio = factors[matrix_index(j,i)]/factors[matrix_index(i,i)];
		            if (debug)
		            	debug_message("ratio [(%d, %d) : (%f, %f)] -> %f\n", i, j, factors[matrix_index(j,i)], factors[matrix_index(i,i)], ratio);

		            for (k = i; k < size; k++) 
		            	factors[matrix_index(j,k)] -= (ratio*factors[matrix_index(i,k)]);
		            values[j] -= (ratio*values[i]);
	            }
	       }
	  	}
  	}		
}

double* pthread_factors;
double* pthread_values;

int pthread_current;
int pthread_row_count;
int pthread_over_row_count;
int pthread_finished;
pthread_barrier_t pthread_barrier;

void pthread_elimination(double* factors, double* values)
{
	pthread_barrier_init(&pthread_barrier, NULL, thread_count+1);

	pthread_factors = factors;
	pthread_values = values;

	pthread_t *pthread_threads = malloc(sizeof(pthread_t) * thread_count);
	int i;
	for(i=0; i<thread_count; i++)	
		pthread_create(&pthread_threads[i], NULL, &pthread_run, (void*)i);

	pthread_finished = 0;
	for(pthread_current = 0; pthread_current < size-1; pthread_current++)
	{
		pthread_row_count = (size - pthread_current - 1)/thread_count;
		pthread_over_row_count = (size - pthread_current - 1)%thread_count;

		pthread_barrier_wait(&pthread_barrier);

		if(pthread_current == size-2)
			pthread_finished = 1;

		pthread_barrier_wait(&pthread_barrier);
	}

	for(i=0; i<thread_count; i++)	
		pthread_join(pthread_threads[i], NULL);
}



void* pthread_run(void* parameter)
{
	int id = (int)parameter;
	int j_start, j_end, j, k;
	double ratio;
	while(!pthread_finished)
	{
		pthread_barrier_wait(&pthread_barrier);

		if (id < pthread_over_row_count)
		{
			j_start = pthread_current + 1 + id*(pthread_row_count+1);
			j_end = pthread_current + 1 + (id+1)*(pthread_row_count+1);
		}
		else
		{
			j_start = pthread_current + 1 + (pthread_over_row_count * (pthread_row_count + 1)) + (id-pthread_over_row_count)*pthread_row_count;
			j_end = pthread_current + 1 + (pthread_over_row_count * (pthread_row_count + 1)) + (id-pthread_over_row_count+1)*pthread_row_count;
		}
		
		if (j_end > size)
			j_end = size;

	    for (j = j_start; j < j_end; j++) 
	    {
	   		if (pthread_factors[matrix_index(pthread_current,pthread_current)] == 0)
	   		{
	   			printf("can not countinue cause (%d,%d) is ziro\n", pthread_current, pthread_current);
	   			return;
	   		}
	   		else
	   		{
	            ratio = pthread_factors[matrix_index(j,pthread_current)]/pthread_factors[matrix_index(pthread_current,pthread_current)];
	            if (debug)
	            	debug_message("thread_id : %d, ratio [(%d, %d) : (%f, %f)] -> %f\n", id, pthread_current, j,
	            		 pthread_factors[matrix_index(j,pthread_current)], pthread_factors[matrix_index(pthread_current,pthread_current)], ratio);

	            for (k = pthread_current; k < size; k++) 
	            	pthread_factors[matrix_index(j,k)] -= (ratio*pthread_factors[matrix_index(pthread_current,k)]);
	            pthread_values[j] -= (ratio*pthread_values[pthread_current]);
	        }
	    }
		pthread_barrier_wait(&pthread_barrier);
	}
}