#include <stdio.h>
#include <stdlib.h> 
#include <unistd.h> 
#include <string.h> 
#include <time.h>
#include <limits.h>
#include <math.h>

long double P = 1.5;
long double Q = 1;
long double STEP = 0.01;
long double NOISE_LVL = 0.00001;
long double** out_crv;	
long double singular_point_x = 1.0;
long double singular_point_y = 1.0;

long double f(long double x, long double y)
{
	return 1-x*y;
}

long double g(long double x, long double y)
{
	return P*y*(x-(1+Q)/(Q+y));
}

long double k1(long double x0, long double y0)
{
	return STEP*f(x0,y0);
}
   
long double l1(long double x0, long double y0)
{
	return STEP*g(x0,y0);
}

long double k2(long double x0, long double y0)
{
    	return STEP*f(x0+k1(x0,y0)/2,y0+l1(x0,y0)/2);	
}
 
long double l2(long double x0, long double y0)
{
    	return STEP*g(x0+k1(x0,y0)/2,y0+l1(x0,y0)/2);	
}   

long double k3(long double x0, long double y0)
{
    	return STEP*f(x0+k2(x0,y0)/2,y0+l2(x0,y0)/2);	
}

long double l3(long double x0, long double y0)
{
    	return STEP*g(x0+k2(x0,y0)/2,y0+l2(x0,y0)/2);	
}

long double k4(long double x0, long double y0)
{
    	return STEP*f(x0+k3(x0,y0)/2,y0+l3(x0,y0)/2);	
}

long double l4(long double x0, long double y0)
{
    	return STEP*g(x0+k3(x0,y0)/2,y0+l3(x0,y0)/2);	
}

void rk(long double* pair)
{
	long double x0 = pair[0];
    	long double y0 = pair[1];
    	long double x = pair[0] +1.0/6*(k1(x0,y0)+2*k2(x0,y0)+2*k3(x0,y0)+k4(x0,y0));
   	long double y = pair[1] +1.0/6*(l1(x0,y0)+2*l2(x0,y0)+2*l3(x0,y0)+l4(x0,y0));
    	pair[0] = x;
    	pair[1] = y;
}
    
long double** curve(long double x0, long double y0, int n)
{
	int i = 0;
  	long double** result = (long double**)malloc(2*sizeof(long double*));
    	result[0] = (long double*)malloc(n*sizeof(long double));
    	result[1] = (long double*)malloc(n*sizeof(long double));
    	long double* current_point = (long double*)malloc(2*sizeof(long double));
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
long double** draw_curve(long double h, long double new_P, long double new_Q, long double x, long double y, int steps)
{
	STEP = h;
	P = new_P;
	Q = new_Q;
	out_crv = curve(x,y,steps);
	return out_crv;
}

long double norm_rand_value_to_abs(long double val, long double abs_val)
{
	long double new_val = val;
	while(fabs(new_val) > abs_val)
	{
		new_val = new_val/10;
	}
	return new_val;
}


void gauss_random_values(long double* result)
{
	long double x = norm_rand_value_to_abs(((long double)abs(rand()))/INT_MAX, 1.0);
	long double y = norm_rand_value_to_abs(((long double)abs(rand()))/INT_MAX, 1.0);
	long double x1 = sqrt(-2*log(x))*cos(2*M_PI*y);
	long double y1 = sqrt(-2*log(x))*sin(2*M_PI*y);
	result[0]=x1;
	result[1]=y1;
}

long double** noisy_rk(long double x0, long double y0, int n)
{
	int i = 0;
  	long double** result = (long double**)malloc(2*sizeof(long double*));
    	result[0] = (long double*)malloc(n*sizeof(long double));
    	result[1] = (long double*)malloc(n*sizeof(long double));
    	long double* current_point = (long double*)malloc(2*sizeof(long double));
    	long double* noise = (long double*)malloc(2*sizeof(long double));
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

long double min(long double a, long double b)
{
	return (a<b)? a : b;
}

long double max(long double a, long double b)
{
	return (a>b)? a : b;
}

long double** find_cycle(long double p, long double q, long double step, long double eps, int iterations)
{
	P = p;
	Q = q;
	STEP = step;
	int is_first_iter = 1;
	int cycle_start_pos;
	int cycle_finish_pos;
	int i = 0;
	long double x0;
	long double y0;
	long double x1;
	long double y1;
	long double cross = 1.1;
	
	long double** crv;
	out_crv = (long double**)malloc(2*sizeof(long double*));
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
			long double buf_x = crv[0][iterations-1];
			long double buf_y = crv[1][iterations-1];
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
             				out_crv[0] = (long double*)malloc((cycle_finish_pos - cycle_start_pos+2)*sizeof(long double));
             				out_crv[1] = (long double*)malloc((cycle_finish_pos - cycle_start_pos+2)*sizeof(long double));
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

long double** draw_noisy_curve(long double noise, long double h, long double new_P, long double new_Q, long double x, long double y, int steps)
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