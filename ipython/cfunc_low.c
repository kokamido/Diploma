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
double** out_crv;	
double singular_point_x = 1.0;
double singular_point_y = 1.0;

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
double** draw_curve(double h, double new_P, double new_Q, double x, double y, int steps)
{
	STEP = h;
	P = new_P;
	Q = new_Q;
	out_crv = curve(x,y,steps);
	return out_crv;
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
	        	x0 = current_point[0]+(noise[0]*NOISE_LVL)*sqrt(STEP);
	        	y0 = current_point[1]+(noise[1]*NOISE_LVL)*sqrt(STEP);
    	}
    	free(current_point);
    	free(noise);
  	return result;
}

double min(double a, double b)
{
	return (a<b)? a : b;
}

double max(double a, double b)
{
	return (a>b)? a : b;
}

double** find_cycle(double p, double q, double step, double eps, int iterations)
{
	P = p;
	Q = q;
	STEP = step;
	int is_first_iter = 1;
	int cycle_start_pos;
	int cycle_finish_pos;
	int i = 0;
	double x0;
	double y0;
	double x1;
	double y1;
	double cross = 1.1;
	
	double** crv;
	out_crv = (double**)malloc(2*sizeof(double*));
	while(1)
	{
		cycle_start_pos = 0;
		cycle_finish_pos = 0;
		if(is_first_iter)
		{
			crv = curve(1.0,1.1,iterations);
			is_first_iter = 0;
		}
		else
		{
			double buf_x = crv[0][iterations-1];
			double buf_y = crv[1][iterations-1];
			free(crv);
			crv = curve(buf_x,buf_y,iterations);
		}
		for(i = 0; i < iterations-1; i++)
		{
			x0 = crv[0][i];
           			y0 = crv[1][i];
             		x1 = crv[0][i+1];
             		y1 = crv[1][i+1];
             		if((y0>=singular_point_y && y1>=singular_point_y) && ((x0-singular_point_x) * (x1-singular_point_x) < 0))
             		{
             			if (fabs(cross-(min(y0,y1)+fabs(y0-y1)/(fabs(x0-x1)/
                                                     fabs((y0<=y1 ? x0 : x1)-singular_point_x))))<eps)
             			{
             				cycle_finish_pos+=1;
             				out_crv[0] = (double*)malloc((cycle_finish_pos - cycle_start_pos+2)*sizeof(double));
             				out_crv[1] = (double*)malloc((cycle_finish_pos - cycle_start_pos+2)*sizeof(double));
             				out_crv[0][0]=cycle_finish_pos - cycle_start_pos+1;
             				out_crv[1][0]=cycle_finish_pos - cycle_start_pos+1;
             				int j;
             				for(j = cycle_start_pos+1; j<=cycle_finish_pos; j++)
             				{
             					out_crv[0][j - cycle_start_pos] = crv[0][j];
             					out_crv[1][j - cycle_start_pos] = crv[1][j];
             				}
             				free(crv[0]);
             				free(crv[1]);
             				free(crv);
             				return out_crv;
             			}
             			else
             			{
             				cross = min(y0,y1)+fabs(y0-y1)/(fabs(x0-x1)/fabs((y0<=y1 ? x0 : x1)-singular_point_x));
             				cycle_start_pos = i;
             				cycle_finish_pos = i;	
             			}
             		}
             		else
             		{
             			cycle_finish_pos+=1;
             		}
		}  
	}
}

double** draw_noisy_curve(double noise, double h, double new_P, double new_Q, double x, double y, int steps)
{
	NOISE_LVL = noise;
	STEP = h;
	P = new_P;
	Q = new_Q;
	out_crv = noisy_rk(x,y,steps);
	return out_crv;
}

void free_crv_memory()
{
	free(out_crv[0]);
	free(out_crv[1]);
	free(out_crv);
}