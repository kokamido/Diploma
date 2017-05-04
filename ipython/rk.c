#include <stdio.h>
#include <stdlib.h> 
#include <unistd.h> 
#include <string.h> 
#include <time.h>
#include <limits.h>
#include <math.h>

double P = 1.5;
double Q = 1;
double STEP = 0.01;
double NOISE_LVL = 0.00001;
	

double f(double x, double y)
{
	return 1-x*y;
}

double g(double x, double y)
{
	return P*y*(x-(1+Q)/(Q+y));
}

double k1(double x0, double y0)
{
	return STEP*f(x0,y0);
}
   
double l1(double x0, double y0)
{
	return STEP*g(x0,y0);
}

double k2(double x0, double y0)
{
    	return STEP*f(x0+k1(x0,y0)/2,y0+l1(x0,y0)/2);	
}
 
double l2(double x0, double y0)
{
    	return STEP*g(x0+k1(x0,y0)/2,y0+l1(x0,y0)/2);	
}   

double k3(double x0, double y0)
{
    	return STEP*f(x0+k2(x0,y0)/2,y0+l2(x0,y0)/2);	
}

double l3(double x0, double y0)
{
    	return STEP*g(x0+k2(x0,y0)/2,y0+l2(x0,y0)/2);	
}

double k4(double x0, double y0)
{
    	return STEP*f(x0+k3(x0,y0)/2,y0+l3(x0,y0)/2);	
}

double l4(double x0, double y0)
{
    	return STEP*g(x0+k3(x0,y0)/2,y0+l3(x0,y0)/2);	
}

void rk(double* pair)
{
	double x0 = pair[0];
    	double y0 = pair[1];
    	double x = pair[0] +1.0/6*(k1(x0,y0)+2*k2(x0,y0)+2*k3(x0,y0)+k4(x0,y0));
   	double y = pair[1] +1.0/6*(l1(x0,y0)+2*l2(x0,y0)+2*l3(x0,y0)+l4(x0,y0));
    	pair[0] = x;
    	pair[1] = y;
}
    
double** curve(double x0, double y0, int n)
{
	int i = 0;
  	double** result = (double**)malloc(2*sizeof(double*));
    	result[0] = (double*)malloc(n*sizeof(double));
    	result[1] = (double*)malloc(n*sizeof(double));
    	double* current_point = (double*)malloc(2*sizeof(double));
    	for(i = 0; i < n; i++)
   	{
 		result[0][i] = x0;
     		result[1][i] = y0;
	        	current_point[0] = x0;
	        	current_point[1] = y0;
	        	rk(current_point);
	        	x0 = current_point[0];
	        	y0 = current_point[1];
    	}
    	free(current_point);
  	return result;
}

//simple r-k curve
int draw_curve(double h, double new_P, double new_Q, double x, double y, int steps)
{
	STEP = h;
	P = new_P;
	Q = new_Q;
	double** crv = (double**)malloc(2*sizeof(double*));
	crv = curve(x,y,steps);
	FILE* output = fopen("C_output.txt","w");
	int i = 0;
	for(i = 0; i <=steps; i++)
	{
		fprintf(output,"%.16lf %.16lf\n", crv[0][i],crv[1][i]);
	}
	free(crv);
	fclose(output);
	return 0;
}

double norm_rand_value_to_abs(double val, double abs_val)
{
	double new_val = val;
	while(fabs(new_val) > abs_val)
	{
		new_val = new_val/10;
	}
	return new_val;
}


void gauss_random_values(double* result)
{
	double x = norm_rand_value_to_abs(((double)abs(rand()))/INT_MAX, 1.0);
	double y = norm_rand_value_to_abs(((double)abs(rand()))/INT_MAX, 1.0);
	double x1 = sqrt(-2*log(x))*cos(2*M_PI*y);
	double y1 = sqrt(-2*log(x))*sin(2*M_PI*y);
	result[0]=x1;
	result[1]=y1;
}

double** noisy_rk(double x0, double y0, int n)
{
	int i = 0;
  	double** result = (double**)malloc(2*sizeof(double*));
    	result[0] = (double*)malloc(n*sizeof(double));
    	result[1] = (double*)malloc(n*sizeof(double));
    	double* current_point = (double*)malloc(2*sizeof(double));
    	double* noise = (double*)malloc(2*sizeof(double));
    	for(i = 0; i < n; i++)
   	{
 		result[0][i] = x0;
     		result[1][i] = y0;
	        	current_point[0] = x0;
	        	current_point[1] = y0;
	        	rk(current_point);
	        	gauss_random_values(noise);
	        	x0 = current_point[0]+norm_rand_value_to_abs(noise[0],NOISE_LVL)*sqrt(STEP);
	        	y0 = current_point[1]+norm_rand_value_to_abs(noise[1],NOISE_LVL)*sqrt(STEP);
    	}
    	free(current_point);
    	free(noise);
  	return result;
}


int draw_noisy_curve(double noise, double h, double new_P, double new_Q, double x, double y, int steps)
{
	NOISE_LVL = noise;
	STEP = h;
	P = new_P;
	Q = new_Q;
	double** crv = (double**)malloc(2*sizeof(double*));
	crv = noisy_rk(x,y,steps);
	FILE* output = fopen("C_output.txt","w");
	int i = 0;
	for(i = 0; i <=steps; i++)
	{
		fprintf(output,"%.16lf %.16lf\n", crv[0][i],crv[1][i]);
	}
	free(crv[0]);
	free(crv[1]);
	free(crv);
	fclose(output);
	return 0;
}

int call_function(int func_id, FILE* input)
{
	if(func_id == 1)
	{
		double new_P;
		double new_Q;
		double x;
		double y;
		double h; // step for R-K
		int steps;
		fscanf(input, "%lf %lf %lf %lf %lf %d",&h, &new_P, &new_Q, &x, &y, &steps);
		return draw_curve(h, new_P, new_Q,x,y,steps);
	}
	if(func_id == 2)
	{
		double noise;
		double new_P;
		double new_Q;
		double x;
		double y;
		double h; // step for R-K
		int steps;
		fscanf(input, "%lf %lf %lf %lf %lf %lf %d",&noise, &h, &new_P, &new_Q, &x, &y, &steps);
		return draw_noisy_curve(noise, h, new_P, new_Q,x,y,steps);
	}
}

    
 //1 - curve
//2 - noisy curve
int main()
{
	time_t t;
	srand((unsigned)time(&t));
	FILE* input = fopen("C_input.txt","r");
	int func_id;
	fscanf(input, "%d\n", &func_id);	
	int exit_code =  call_function(func_id, input);
	fclose(input);
	return exit_code;
}
	    

