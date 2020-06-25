#if !defined(_MSC_VER)
#include <config.h>
#else
#include <config_msvc.h>
#endif

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <minmax.h>
#include <models/equations.h>

double sq(double x) { return x * x; }

//Type 252
//Contains 3 layers on hillslope: ponded, top layer, soil
//Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
//The numbering is:	0   1   2     3    4   5   6   7
//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent
//The numbering is:        0      1        2     3   4     5        6   7  8 9  10
void TopLayerHillslope(double t, const double *const y_i, unsigned int dim, const double *const y_p, unsigned short num_parents, unsigned int max_dim, const double *const global_params, const double *const params, const double *const forcing_values, const QVSData *const qvs, int state, void *user, double *ans, double h)
{
	unsigned short i;

	double lambda_1 = global_params[1];
	double k_3 = global_params[4]; //[1/min]
	double h_b = global_params[6]; //[m]
	double S_L = global_params[7]; //[m]
	double A = global_params[8];
	double B = global_params[9];
	double exponent = global_params[10];
	double e_pot = forcing_values[1] * (1e-3 / (30.0 * 24.0 * 60.0)); //[mm/month] -> [m/min]
																	  //double e_pot = global_params[11];	//[m/min]
																	  //double e_pot = global_params[11] * (1e-3*60.0);	//[m/min]
																	  //double e_pot = 0.0;

	double L = params[1];	   //[m]
	double A_h = params[2];	   //[m^2]
							   //double h_r = params[3];	//[m]
	double invtau = params[3]; //[1/min]
	double k_2 = params[4];	   //[1/min]
	double k_i = params[5];	   //[1/min]
	double c_1 = params[6];
	double c_2 = params[7];

	double q = y_i[0];	 //[m^3/s]
	double s_p = y_i[1]; //[m]
	double s_t = y_i[2]; //[m]
	double s_s = y_i[3]; //[m]

	//Evaporation
	double e_p, e_t, e_s;
	double Corr = s_p + s_t / S_L + s_s / (h_b - S_L);
	if (e_pot > 0.0 && Corr > 1e-12)
	{
		e_p = s_p * 1e3 * e_pot / Corr;
		e_t = s_t / S_L * e_pot / Corr;
		e_s = s_s / (h_b - S_L) * e_pot / Corr;
	}
	else
	{
		e_p = 0.0;
		e_t = 0.0;
		e_s = 0.0;
	}

	double pow_term = (1.0 - s_t / S_L > 0.0) ? pow(1.0 - s_t / S_L, exponent) : 0.0;
	double k_t = (A + B * pow_term) * k_2;

	//Fluxes
	double q_pl = k_2 * s_p;
	double q_pt = k_t * s_p;
	double q_ts = k_i * s_t;
	double q_sl = k_3 * s_s;

	//Discharge
	ans[0] = -q + (q_pl + q_sl) * c_2;
	for (i = 0; i < num_parents; i++)
		ans[0] += y_p[i * dim];
	ans[0] = invtau * pow(q, lambda_1) * ans[0];

	//Hillslope
	ans[1] = forcing_values[0] * c_1 - q_pl - q_pt - e_p;
	ans[2] = q_pt - q_ts - e_t;
	ans[3] = q_ts - q_sl - e_s;
}

//Type 254
//Contains 3 layers on hillslope: ponded, top layer, soil. Also has 3 extra states: total precip, total runoff, base flow
//Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
//The numbering is:	0   1   2     3    4   5   6   7
//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent,v_B
//The numbering is:        0      1        2     3   4     5        6   7  8 9  10       11
void TopLayerHillslope_extras(double t, const double *const y_i, unsigned int dim, const double *const y_p, unsigned short num_parents, unsigned int max_dim, const double *const global_params, const double *const params, const double *const forcing_values, const QVSData *const qvs, int state, void *user, double *ans, double h)
{
	unsigned short i;

	double lambda_1 = global_params[1];
	double k_3 = global_params[4]; //[1/min]
	double h_b = global_params[6]; //[m]
	double S_L = global_params[7]; //[m]
	double A = global_params[8];
	double B = global_params[9];
	double exponent = global_params[10];
	double v_B = global_params[11];
	double e_pot = forcing_values[1] * (1e-3 / (30.0 * 24.0 * 60.0)); //[mm/month] -> [m/min]

	double L = params[1];	   //[m]
	double A_h = params[2];	   //[m^2]
							   //double h_r = params[3];	//[m]
	double invtau = params[3]; //[1/min]
	double k_2 = params[4];	   //[1/min]
	double k_i = params[5];	   //[1/min]
	double c_1 = params[6];
	double c_2 = params[7];

	double q = y_i[0];	 //[m^3/s]
	double s_p = y_i[1]; //[m]
	double s_t = y_i[2]; //[m]
	double s_s = y_i[3]; //[m]
						 //double s_precip = y_i[4];	//[m]
						 //double V_r = y_i[5];	//[m^3]
	double q_b = y_i[6]; //[m^3/s]

	//Evaporation
	double e_p, e_t, e_s;
	double Corr = s_p + s_t / S_L + s_s / (h_b - S_L);
	if (e_pot > 0.0 && Corr > 1e-12)
	{
		e_p = s_p * 1e3 * e_pot / Corr;
		e_t = s_t / S_L * e_pot / Corr;
		e_s = s_s / (h_b - S_L) * e_pot / Corr;
	}
	else
	{
		e_p = 0.0;
		e_t = 0.0;
		e_s = 0.0;
	}

	double pow_term = (1.0 - s_t / S_L > 0.0) ? pow(1.0 - s_t / S_L, exponent) : 0.0;
	double k_t = (A + B * pow_term) * k_2;

	//Fluxes
	double q_pl = k_2 * s_p;
	double q_pt = k_t * s_p;
	double q_ts = k_i * s_t;
	double q_sl = k_3 * s_s; //[m/min]

	//Discharge
	ans[0] = -q + (q_pl + q_sl) * c_2;
	for (i = 0; i < num_parents; i++)
		ans[0] += y_p[i * dim];
	ans[0] = invtau * pow(q, lambda_1) * ans[0];

	//Hillslope
	ans[1] = forcing_values[0] * c_1 - q_pl - q_pt - e_p;
	ans[2] = q_pt - q_ts - e_t;
	ans[3] = q_ts - q_sl - e_s;

	//Additional states
	ans[4] = forcing_values[0] * c_1;
	ans[5] = q_pl;
	ans[6] = q_sl * A_h - q_b * 60.0;
	for (i = 0; i < num_parents; i++)
		ans[6] += y_p[i * dim + 6] * 60.0;
	//ans[6] += k_3*y_p[i].ve[3]*A_h;
	ans[6] *= v_B / L;
}

//Type 253 / 255
//Contains 3 layers on hillslope: ponded, top layer, soil
//Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
//The numbering is:	0   1   2     3    4   5   6   7
//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent
//The numbering is:        0      1        2     3   4     5        6   7  8 9  10
void TopLayerHillslope_Reservoirs(double t, const double *const y_i, unsigned int dim, const double *const y_p, unsigned short num_parents, unsigned int max_dim, const double *const global_params, const double *const params, const double *const forcing_values, const QVSData *const qvs, int state, void *user, double *ans, double h)
{
	ans[0] = forcing_values[2];
	ans[1] = 0.0;
	ans[2] = 0.0;
	ans[3] = 0.0;
}

//Type 8
//Contains 3 layers on hillslope: ponded, top layer, soil. Also has 3 extra states: total precip, total runoff, base flow
//Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
//The numbering is:	0   1   2     3    4   5   6   7
//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent,v_B
//The numbering is:        0      1        2     3   4     5        6   7  8 9  10       11
void Eight_TopLayers(double t, const double *const y_i, unsigned int dim, const double *const y_p, unsigned short num_parents, unsigned int max_dim, const double *const global_params, const double *const params, const double *const forcing_values, const QVSData *const qvs, int state, void *user, double *ans, double h)
{
	unsigned short i;
	int num_T_layers = 8;
	double lambda_1 = global_params[1];
	double k_3 = global_params[4]; //[1/min]
	double h_b = global_params[6]; //[m]
	double S_L = global_params[7]; //[m]
	double A = global_params[8];
	double B = global_params[9];
	double exponent = global_params[10];
	double v_B = global_params[11];
	double e_pot = forcing_values[1] * (1e-3 / (30.0 * 24.0 * 60.0)); //[mm/month] -> [m/min]

	double L = params[1];	   //[m]
	double A_h = params[2];	   //[m^2]
							   //double h_r = params[3];	//[m]
	double invtau = params[3]; //[1/min]
	double k_2 = params[4];	   //[1/min]
	double k_i = params[5];	   //[1/min]
	double c_1 = params[6];
	double c_2 = params[7];

	double q = y_i[0];	  //[m^3/s]
	double s_p = y_i[1];  //[m]
	double s_t1 = y_i[2]; //[m]
	double s_t2 = y_i[3]; //[m]
	double s_t3 = y_i[4]; //[m]
	double s_t4 = y_i[5]; //[m]
	double s_t5 = y_i[6]; //[m]
	double s_t6 = y_i[7]; //[m]
	double s_t7 = y_i[8]; //[m]
	double s_t8 = y_i[9]; //[m]
	double s_s = y_i[10]; //[m]
	double q_t[7];
	//double s_precip = y_i[4];	//[m]
	//double V_r = y_i[5];	//[m^3]
	double q_b = y_i[11]; //[m^3/s]

	//Evaporation
	double e_p, e_t, e_s;
	double Corr = s_p + s_t1 / S_L + s_s / (h_b - S_L);
	if (e_pot > 0.0 && Corr > 1e-12)
	{
		e_p = s_p * 1e3 * e_pot / Corr;
		e_t = s_t1 / S_L * e_pot / Corr;
		e_s = s_s / (h_b - S_L) * e_pot / Corr;
	}
	else
	{
		e_p = 0.0;
		e_t = 0.0;
		e_s = 0.0;
	}

	double pow_term = (1.0 - s_t1 / S_L > 0.0) ? pow(1.0 - s_t1 / S_L, exponent) : 0.0;
	double k_t = (A + B * pow_term) * k_2;
	double et_t = e_t / num_T_layers;
	//Fluxes
	double q_pl = k_2 * s_p;
	double q_pt1 = k_t * s_p;
	//double infilt = k_t * s_p;
	//q_t[0] = k_i  * s_t1;
	//q_t[1] = k_t  * s_t2;
	//q_t[2] = k_t  * s_t3;
	//q_t[3] = k_t  * s_t4;
	//q_t[4] = k_t  * s_t5;
	//q_t[5] = k_t  * s_t6;
	//q_t[6] = k_t  * s_t7;

	if (forcing_values[0] > 0) // Rainfall period
	{
		q_t[0] = s_t1 * 1e-6;
		q_t[1] = s_t2 * 1e-6;
		q_t[2] = s_t3 * 1e-6;
		q_t[3] = s_t4 * 1e-6;
		q_t[4] = s_t5 * 1e-6;
		q_t[5] = s_t6 * 1e-6;
		q_t[6] = s_t7 * 1e-6;
	}
	else if (forcing_values[0] == 0) // Dry-out  period
	{

		// test_2
		//q_t[0] = s_t1 * (1 - s_t2 / S_L)*0.02;
		//q_t[1] = s_t2 * (1 - s_t3 / S_L)*0.02;
		//q_t[2] = s_t3 * (1 - s_t4 / S_L)*0.02;
		//q_t[3] = s_t4 * (1 - s_t5 / S_L)*0.02;
		//q_t[4] = s_t5 * (1 - s_t6 / S_L)*0.02;
		//q_t[5] = s_t6 * (1 - s_t7 / S_L)*0.02;
		//q_t[6] = s_t7 * (1 - s_t8 / S_L)*0.02;

		//test_3
		//q_t[0] = s_t1 *	0.02 / (24 * 60);
		//q_t[1] = s_t2 * 0.01 / (24 * 60);
		//q_t[2] = s_t3 * 0.005 / (24 * 60);
		//q_t[3] = s_t4 * 0.0001 / (24 * 60);
		//q_t[4] = s_t5 * 0.0008 / (24 * 60);
		//q_t[5] = s_t6 * 0.00005 / (24 * 60);
		//q_t[6] = s_t7 * 0.00002 / (24 * 60);

		// test_4
		//q_t[0] = (s_t1 - s_t2) / (S_L) *  1.73e-4;
		//q_t[1] = (s_t2 - s_t3) / (S_L) *  1.73e-5;
		//q_t[2] = (s_t3 - s_t4) / (S_L) *  1.73e-6;
		//q_t[3] = (s_t4 - s_t5) / (S_L) *  8.73e-6;
		//q_t[4] = (s_t5 - s_t6) / (S_L) *  5.93e-7;
		//q_t[5] = (s_t6 - s_t7) / (S_L) *  9.30e-7;
		//q_t[6] = (s_t7 - s_t8) / (S_L) *  1.93e-8;

		//test_1
		//q_t[0] = (s_t1 - s_t2) / (S_L) *  1.73e-4;
		//q_t[1] = (s_t2 - s_t3) / (S_L) *  1.73e-4;
		//q_t[2] = (s_t3 - s_t4) / (S_L) *  1.73e-4;
		//q_t[3] = (s_t4 - s_t5) / (S_L) *  8.73e-4;
		//q_t[4] = (s_t5 - s_t6) / (S_L) *  5.93e-4;
		//q_t[5] = (s_t6 - s_t7) / (S_L) *  9.30e-4;
		//q_t[6] = (s_t7 - s_t8) / (S_L) *  1.93e-4;

		// test_5
		q_t[0] = s_t1 > s_t2 ? (s_t1 - s_t2) / (S_L)*1.73e-4 : 0;
		q_t[1] = s_t2 > s_t3 ? (s_t2 - s_t3) / (S_L)*1.73e-5 : 0;
		q_t[2] = s_t3 > s_t4 ? (s_t3 - s_t4) / (S_L)*1.73e-6 : 0;
		q_t[3] = s_t4 > s_t5 ? (s_t4 - s_t5) / (S_L)*8.73e-6 : 0;
		q_t[4] = s_t5 > s_t6 ? (s_t5 - s_t6) / (S_L)*5.93e-7 : 0;
		q_t[5] = s_t6 > s_t7 ? (s_t6 - s_t7) / (S_L)*9.30e-7 : 0;
		q_t[6] = s_t7 > s_t8 ? (s_t7 - s_t8) / (S_L)*1.93e-8 : 0;
	}
	/*1 / (3 * 60);
	1 / (5 * 60);
	1 / (10 * 60);
	1 / (20 * 60);
	1 / (24 * 60);
	1 / (48 * 60);
	1 / (96 * 60);*/

	//double q_ts = s_t8 * (1 - (s_s / (h_b - S_L*num_T_layers)))*k_3;
	double q_ts = s_t8 * 1.93e-6;
	//double q_sl = k_3 * s_s;
	double q_sl = k_3 * s_s; //[m/min]
	////Uniform ET from all layers

	//Discharge
	ans[0] = -q + (q_pl + q_sl) * c_2;
	for (i = 0; i < num_parents; i++)
		ans[0] += y_p[i * dim];
	ans[0] = invtau * pow(q, lambda_1) * ans[0];

	//Hillslope
	ans[1] = forcing_values[0] * c_1 - q_pl - q_pt1 - e_p; //S_p
	ans[2] = q_pt1 - q_t[0] - e_t;						   //S_t1
	ans[3] = q_t[0] - q_t[1];							   //S_t2
	ans[4] = q_t[1] - q_t[2];							   //S_t3
	ans[5] = q_t[2] - q_t[3];							   //S_t4
	ans[6] = q_t[3] - q_t[4];							   //S_t5
	ans[7] = q_t[4] - q_t[5];							   //S_t6
	ans[8] = q_t[5] - q_t[6];							   //S_t7
	ans[9] = q_t[6] - q_ts;								   //S_t8
	ans[10] = q_ts - q_sl;								   //S_s

	//Additional states
	ans[11] = forcing_values[0] * c_1;
	ans[12] = q_pl;
	ans[13] = q_sl * A_h - q_b * 60.0;
	for (i = 0; i < num_parents; i++)
		ans[13] += y_p[i * dim + 13] * 60.0;
	//ans[6] += k_3*y_p[i].ve[3]*A_h;
	ans[13] *= v_B / L;
}

//Type 1
//Contains 3 layers on hillslope: ponded, top layer, soil. Also has 3 extra states: total precip, total runoff, base flow
//Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
//The numbering is:	0   1   2     3    4   5   6   7
//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent,v_B
//The numbering is:        0      1        2     3   4     5        6   7  8 9  10       11
void Soil_moisture_velocity(double t, const double *const y_i, unsigned int dim, const double *const y_p, unsigned short num_parents, unsigned int max_dim, const double *const global_params, const double *const params, const double *const forcing_values, const QVSData *const qvs, int state, void *user, double *ans, double h)
{
	unsigned short i;
	int num_T_layers = 8;
	double lambda_1 = global_params[1];
	double k_3 = global_params[4]; //[1/min]
	double h_b = global_params[6]; //[m]
	double S_L = global_params[7]; //[m]
	double A = global_params[8];
	double B = global_params[9];
	double exponent = global_params[10];
	double v_B = global_params[11];
	double e_pot = forcing_values[1] * (1e-3 / (30.0 * 24.0 * 60.0)); //[mm/month] -> [m/min]

	double L = params[1];	   //[m]
	double A_h = params[2];	   //[m^2]
							   //double h_r = params[3];	//[m]
	double invtau = params[3]; //[1/min]
	double k_2 = params[4];	   //[1/min]
	double k_i = params[5];	   //[1/min]
	double c_1 = params[6];
	double c_2 = params[7];

	double q = y_i[0];	  //[m^3/s]
	double s_p = y_i[1];  //[m]
	double s_t1 = y_i[2]; //[m]
	double s_t2 = y_i[3]; //[m]
	double s_t3 = y_i[4]; //[m]
	double s_t4 = y_i[5]; //[m]
	double s_t5 = y_i[6]; //[m]
	double s_t6 = y_i[7]; //[m]
	double s_t7 = y_i[8]; //[m]
	double s_t8 = y_i[9]; //[m]
	double s_s = y_i[10]; //[m]
	double q_t[8];
	//double s_precip = y_i[4];	//[m]
	//double V_r = y_i[5];	//[m^3]
	double q_b = y_i[11]; //[m^3/s]

	//Evaporation
	double e_p, e_t, e_s;
	double Corr = s_p + s_t1 / S_L + s_s / (h_b - S_L);
	if (e_pot > 0.0 && Corr > 1e-12)
	{
		e_p = s_p * 1e3 * e_pot / Corr;
		e_t = s_t1 / S_L * e_pot / Corr;
		e_s = s_s / (h_b - S_L) * e_pot / Corr;
	}
	else
	{
		e_p = 0.0;
		e_t = 0.0;
		e_s = 0.0;
	}

	double pow_term = (1.0 - s_t1 / S_L > 0.0) ? pow(1.0 - s_t1 / S_L, exponent) : 0.0;
	double k_t = (A + B * pow_term) * k_2;
	double et_t = e_t / num_T_layers;
	//Fluxes
	double q_pl = k_2 * s_p;
	double q_pt1 = k_t * s_p;
	//double infilt = k_t * s_p;
	//q_t[0] = k_i  * s_t1;
	//q_t[1] = k_t  * s_t2;
	//q_t[2] = k_t  * s_t3;
	//q_t[3] = k_t  * s_t4;
	//q_t[4] = k_t  * s_t5;
	//q_t[5] = k_t  * s_t6;
	//q_t[6] = k_t  * s_t7;

	q_t[0] = s_t1 * (1 - s_t2 / S_L) * 6.47E-5;
	q_t[1] = s_t2 * (1 - s_t3 / S_L) * 6.00E-5;
	q_t[2] = s_t3 * (1 - s_t4 / S_L) * 5.47E-5;
	q_t[3] = s_t4 * (1 - s_t5 / S_L) * 5.00E-5;
	q_t[4] = s_t5 * (1 - s_t6 / S_L) * 4.47E-5;
	q_t[5] = s_t6 * (1 - s_t7 / S_L) * 4.00E-5;
	q_t[6] = s_t7 * (1 - s_t8 / S_L) * 3.47E-5;

	double q_ts = s_t8 * (1 - (s_s / (h_b - S_L * num_T_layers))) * k_3;
	double q_sl = k_3 * s_s; //[m/min]
							 //Uniform ET from all layers

	//Discharge
	ans[0] = -q + (q_pl + q_sl) * c_2;
	for (i = 0; i < num_parents; i++)
		ans[0] += y_p[i * dim];
	ans[0] = invtau * pow(q, lambda_1) * ans[0];

	//Hillslope
	ans[1] = forcing_values[0] * c_1 - q_pl - q_pt1 - e_p; //S_p
	ans[2] = q_pt1 - q_t[0] - et_t;						   //S_t1
	ans[3] = q_t[0] - q_t[1] - et_t;					   //S_t2
	ans[4] = q_t[1] - q_t[2] - et_t;					   //S_t3
	ans[5] = q_t[2] - q_t[3] - et_t;					   //S_t4
	ans[6] = q_t[3] - q_t[4] - et_t;					   //S_t5
	ans[7] = q_t[4] - q_t[5] - et_t;					   //S_t6
	ans[8] = q_t[5] - q_t[6] - et_t;					   //S_t7
	ans[9] = q_t[6] - q_ts - et_t;						   //S_t8
	ans[10] = q_ts - e_s - q_sl;						   //S_s

	//Additional states
	ans[11] = forcing_values[0] * c_1;
	ans[12] = q_pl;
	ans[13] = q_sl * A_h - q_b * 60.0;
	for (i = 0; i < num_parents; i++)
		ans[13] += y_p[i * dim + 13] * 60.0;
	//ans[6] += k_3*y_p[i].ve[3]*A_h;
	ans[13] *= v_B / L;
}

double sum(double arr[])
{

	int i;
	int n = (int)sizeof(arr);
	double sum = 0;
	for (i = 0; i < n; i++)
	{
		sum = sum + arr[i];
	}
	return sum;
}

// Van Genuchten model
double K_UNSAT(double s_t, double K_SAT, double m_VG, double S_L, int num_T_layers)
{
	double theta = (s_t / S_L) > 1 ? 1 : (s_t / S_L);
	return K_SAT * pow(theta, 0.5) * pow(1 - pow(1 - pow(theta, 1 / m_VG), m_VG), 2);
}

// Model_uid = 1004
void Four_TopLayers(double t, const double *const y_i, unsigned int dim, const double *const y_p, unsigned short num_parents, unsigned int max_dim, const double *const global_params, const double *const params, const double *const forcing_values, const QVSData *const qvs, int state, void *user, double *ans, double h)
{
	double q = y_i[0];	  //[m^3/s]
	double s_p = y_i[1];  //[m]
	double s_t1 = y_i[2]; //[m]
	double s_t2 = y_i[3]; //[m]
	double s_t3 = y_i[4]; //[m]
	double s_t4 = y_i[5]; //[m]
	double s_s = y_i[6];  //[m]
	double q_b = y_i[9];  //[m^3/s]
	double q_t[3];
	unsigned short i;

	int num_T_layers = 4;
	double lambda_1 = global_params[1];
	double k_3 = global_params[4]; //[1/min]
	double h_b = global_params[6]; //[m]
	double T_L = global_params[7]; //[m] total layer depth physically air, soil and water
	double A = global_params[8];
	double B = global_params[9];
	double exponent = global_params[10];
	double v_B = global_params[11];

	double e_pot = forcing_values[1] * (1e-3 / (30.0 * 24.0 * 60.0)); //[mm/month] -> [m/min]

	double L = params[1];	//[m]
	double A_h = params[2]; //[m^2]
	double K_SAT = params[3];
	double m_VG = params[4];
	double SM_sat = params[5];
	double SM_r = params[6];
	double invtau = params[7]; //[1/min]
	double k_2 = params[8];	   //[1/min]
	double k_i = params[9];	   //[1/min]
	double c_1 = params[10];
	double c_2 = params[11];

	double e_p, e_t, e_s;

	//double D_col = (SM_sat - SM_r) * h_b; // depth of available storage for water in total column
	double D_col = 0.5 * h_b; // depth of available storage for water in total column

	//double S_L = (SM_sat - SM_r) * T_L; // depth of available water storage in Top Layer ( which can be more than one)
	double S_L = 0.5 * T_L; // depth of available water storage in Top Layer ( which can be more than one)

	double st_tot = s_t1 + s_t2 + s_t3 + s_t4;
	double theta_tot = st_tot / S_L;
	//double =

	double Corr = s_p + st_tot / S_L + s_s / (D_col - S_L);
	if (e_pot > 0.0 && Corr > 1e-12)
	{
		e_p = s_p * e_pot / Corr;
		e_t = st_tot / S_L * e_pot / Corr;
		e_s = s_s / (D_col - S_L) * e_pot / Corr;
	}
	else
	{
		e_p = 0.0;
		e_t = 0.0;
		e_s = 0.0;
	}

	double et_t = e_t / num_T_layers;
	//q_t[0] = (s_t1 >0) ? s_t1 / (S_L / num_T_layers) * K_UNSAT(s_t1, K_SAT, m_VG, S_L, num_T_layers) * 0.01 / 60.0 : 0;
	//q_t[1] = (s_t2 >0) ? s_t2 / (S_L / num_T_layers) * K_UNSAT(s_t2, K_SAT, m_VG, S_L, num_T_layers) * 0.01 / 60.0 * 0.1 : 0;
	//q_t[2] = (s_t3 >0) ? s_t3 / (S_L / num_T_layers) * K_UNSAT(s_t3, K_SAT, m_VG, S_L, num_T_layers) * 0.01 / 60.0 * 0.01 : 0;

	//q_t[0] = k_i* s_t1;
	//q_t[1] = k_i* s_t2;
	//q_t[2] = k_i* s_t3;
	/*q_t[0] = (s_t1> s_t2) ? (s_t1 - s_t2) / (S_L) *  1.73e-4 : s_t1 * 1e-6;
	q_t[1] = (s_t2> s_t3) ? (s_t2 - s_t3) / (S_L) *  1.73e-4 : s_t2 * 1e-6;
	q_t[2] = (s_t3> s_t4) ? (s_t3 - s_t4) / (S_L) *  8.73e-4 : s_t3 * 1e-6;
*/
	//q_t[0] = K_UNSAT(s_t1, K_SAT, m_VG, S_L, num_T_layers) * 0.01 / 60.0 * 0.01;
	//q_t[1] = K_UNSAT(s_t2, K_SAT, m_VG, S_L, num_T_layers) * 0.01 / 60.0 * 0.001;
	//q_t[2] = K_UNSAT(s_t3, K_SAT, m_VG, S_L, num_T_layers) * 0.01 / 60.0 * 0.0001;

	//double q_t1 = (s_t1 - s_t2) / (0.025) *  1.73e-4;
	//double q_t2 = (s_t2 - s_t3) / (0.025) *  1.73e-4;
	//double q_t3 = (s_t3 - s_t4) / (0.025) *  8.73e-4;
	double pow_term = (1.0 - s_t1 / S_L > 0.0) ? pow(1.0 - s_t1 / S_L, exponent) : 0.0;
	double k_t = (A + B * pow_term) * k_2;

	double q_pl = k_2 * s_p;
	double q_pt1 = k_t * s_p;
	//double q_t1 = s_t1 * (1 - s_t2 / 0.05)*6.47E-5;
	//double q_t2 = s_t2 * (1 - s_t3 / 0.05)*6.00E-6;
	//double q_t3 = s_t3 * (1 - s_t4 / 0.05)*3.47E-6;

	if (forcing_values[0] > 0) // Rainfall period
	{
		q_t[0] = s_t1 * 1e-6;
		q_t[1] = s_t2 * 1e-6;
		q_t[2] = s_t3 * 1e-6;
	}
	else if (forcing_values[0] == 0) // Dry-out  period
	{
		q_t[0] = (s_t1 - s_t2) / (0.05) * K_UNSAT(s_t1, K_SAT, m_VG, 0.05, num_T_layers) * 0.01 / 60.0;
		q_t[1] = (s_t2 - s_t3) / (0.05) * K_UNSAT(s_t1, K_SAT, m_VG, 0.05, num_T_layers) * 0.01 / 60.0;
		q_t[2] = (s_t3 - s_t4) / (0.05) * K_UNSAT(s_t1, K_SAT, m_VG, 0.05, num_T_layers) * 0.01 / 60.0;
	}
	double q_ts = s_t4 * (1 - (s_s / (D_col - S_L))) * k_3; // k_i * s_t4;
	double q_sl = k_3 * s_s;

	//Discharge
	ans[0] = -q + (q_pl + q_sl) * c_2;
	for (i = 0; i < num_parents; i++)
		ans[0] += y_p[i * dim];
	ans[0] = invtau * pow(q, lambda_1) * ans[0];

	//Hillslope
	ans[1] = forcing_values[0] * c_1 - q_pl - q_pt1 - e_p; //S_p

	// Soil Moisture
	ans[2] = q_pt1 - q_t[0];  //S_t1
	ans[3] = q_t[0] - q_t[1]; //S_t2
	ans[4] = q_t[1] - q_t[2]; //S_t3
	ans[5] = q_t[2] - q_ts;	  //S_t4
	ans[6] = q_ts - q_sl;	  //S_s

	//Additional states
	ans[7] = forcing_values[0] * c_1;
	ans[8] = q_pl;
	ans[9] = q_sl * A_h - q_b * 60.0;
	for (i = 0; i < num_parents; i++)
		ans[9] += y_p[i * dim + 9] * 60.0;
	ans[9] *= v_B / L;
}

float flux(float theta_1, float theta_2, float K_sat, float psi_sat, float lambda, float theta_s, float theta_r, float S_L)
{
	//theta_1 = theta_1 <= theta_r ? theta_r + 0.01 : theta_1;
	//theta_2 = theta_2 <= theta_r ? theta_r + 0.01 : theta_2;
	//theta_1 = theta_1 > theta_s ? theta_s : theta_1;
	//theta_2 = theta_2 > theta_s ? theta_s : theta_2;
	float theta_hat = (theta_2 + theta_1) / 2.0;
	float phi_1 = (theta_1 - theta_r) / (theta_s - theta_r);
	float phi_2 = (theta_2 - theta_r) / (theta_s - theta_r);

	// K
	float phi_k = (theta_hat - theta_r) / (theta_s - theta_r);
	float power_k = (2.0 + 3.0 * lambda) / lambda;
	float K_phi_hat = K_sat * pow(phi_k, power_k);

	// psi
	float power_psi = -1.0 / lambda;
	float psi_1 = psi_sat * pow(phi_1, power_psi);
	float psi_2 = psi_sat * pow(phi_2, power_psi);

	return K_phi_hat * (1.0 - (psi_2 - psi_1) / S_L);
}

// MODEL UID = 1005
//Type 8 based on the paper Herrada Martin 2014 JHydrol http://dx.doi.org/10.1016/j.jhydrol.2014.04.026
void GAMPT_extras(double t, const double *const y_i, unsigned int dim, const double *const y_p, unsigned short num_parents, unsigned int max_dim, const double *const global_params, const double *const params, const double *const forcing_values, const QVSData *const qvs, int state, void *user, double *ans, double h)
{
	unsigned short i;

	// Model parameters (spatially uniform)
	double lambda_1 = global_params[1];
	double k_3 = global_params[4]; //[1/min]
	double h_b = global_params[6]; //[m]
	double S_L = global_params[7]; //[m]
	double A = global_params[8];
	double B = global_params[9];
	double exponent = global_params[10];
	double v_B = global_params[11];
	double psi_sat = global_params[12];

	// Model parameters (spatially variable)
	double L = params[1];	//[m]
	double A_h = params[2]; //[m^2]
	double K_sat = params[3];
	double invtau = params[4]; //[1/min]
	double k_2 = params[5];	   //[1/min]
	double k_i = params[6];	   //[1/min]
	double c_1 = params[7];
	double c_2 = params[8];
	double theta_s = params[9];
	double theta_r = params[10];
	double bc_lambda = params[11];

	double L_Top = S_L * 10;

	// Initial conditions (or from last iteration)
	double q = y_i[0];	 //[m^3/s]
	double s_p = y_i[1]; //[m]
	double s_1 = y_i[2]; //[m]
	double s_2 = y_i[3]; //[m]
	double s_3 = y_i[4];
	double s_4 = y_i[5];
	double s_5 = y_i[6];
	double s_6 = y_i[7];
	double s_7 = y_i[8];
	double s_8 = y_i[9];
	double s_9 = y_i[10];
	double s_10 = y_i[11];
	double s_s = y_i[12];

	double q_b = y_i[15]; //[m^3/s]

	// Potential Evapotranspiration
	double e_pot = forcing_values[1] * (1e-3 / (30.0 * 24.0 * 60.0)); //[mm/month] -> [m/min]
	//Evaporation
	double e_p, e_t, e_s;
	double avg_s_t = (s_1 + s_2 + s_3 + s_4 + s_5 + s_6 + s_7 + s_8 + s_9 + s_10) / 10;
	double Corr = s_p + avg_s_t / L_Top + s_s / (h_b - L_Top);
	if (e_pot > 0.0 && Corr > 1e-12)
	{
		e_p = s_p * e_pot / Corr;
		e_t = s_1 * e_pot / Corr;
		e_s = s_s / (h_b - L_Top) * e_pot / Corr;
	}
	else
	{
		e_p = 0.0;
		e_t = 0.0;
		e_s = 0.0;
	}
	//Fluxes

	//double pow_term = (1.0 - s_1 / S_L > 0.0) ? pow(1.0 - s_1 / S_L, exponent) : 0.0;
	//double k_t = (A + B * pow_term) * k_2;
	double s_1_hat = (s_1 - theta_r) / (theta_s - theta_r);
	double pow_term = (1.0 - s_1_hat > 0.0) ? pow(1.0 - s_1_hat, exponent) : 0.0;
	double k_t = (A + B * pow_term) * k_2;
	double et_t = e_t / 10;
	double q_pl = k_2 * s_p;
	double q_pt = k_t * s_p;
	double q_ts = s_10 * S_L * k_3;
	double q_sl = k_3 * s_s;				 //[m/min]
	double q_rain = forcing_values[0] * c_1; // m/min
	//double q_excess;
	double q_2 = flux(s_1, s_2, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);	 // [m/min]
	double q_3 = flux(s_2, s_3, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);	 // [m/min]
	double q_4 = flux(s_3, s_4, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);	 // [m/min]
	double q_5 = flux(s_4, s_5, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);	 // [m/min]
	double q_6 = flux(s_5, s_6, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);	 // [m/min]
	double q_7 = flux(s_6, s_7, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);	 // [m/min]
	double q_8 = flux(s_7, s_8, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);	 // [m/min]
	double q_9 = flux(s_8, s_9, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);	 // [m/min]
	double q_10 = flux(s_9, s_10, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L); // [m/min]
	//s_10 = (4 * s_9 - s_8) / 3;
	//double q_ts = k_i * s_10 * S_L; // K_i * s_10 [m/min]

	//double dsp, ds1;
	//if (s_1 >= theta_s && q_rain >= q_2) {
	//	q_excess = q_rain - q_2;
	//		s_1 = theta_s;
	//		ds1 = 0;
	//}
	//else {
	//	q_excess = 0.0;
	//	ds1 = (q_rain - q_2) / S_L;
	//}

	// dx/dt ODE definition for solver (RK45)
	//Discharge
	ans[0] = -q + (q_pl + q_sl) * c_2;
	for (i = 0; i < num_parents; i++)
		ans[0] += y_p[i * dim];
	ans[0] = invtau * pow(q, lambda_1) * ans[0];

	//Hillslope
	ans[1] = q_rain - q_pt - q_pl;				 // -e_p;				// S_p
	ans[2] = q_pt / S_L - q_2 / S_L - e_t / S_L; // S_1
	ans[3] = (q_2 - q_3) / S_L - et_t / S_L;	 // S_2
	ans[4] = (q_3 - q_4) / S_L - et_t / S_L;	 // S_3
	ans[5] = (q_4 - q_5) / S_L - et_t / S_L;	 // S_4
	ans[6] = (q_5 - q_6) / S_L - et_t / S_L;	 // S_5
	ans[7] = (q_6 - q_7) / S_L - et_t / S_L;	 // S_6
	ans[8] = (q_7 - q_8) / S_L - et_t / S_L;	 // S_7
	ans[9] = (q_8 - q_9) / S_L - et_t / S_L;	 // S_8
	ans[10] = (q_9 - q_10) / S_L - et_t / S_L;	 // S_9
	ans[11] = (q_10 - q_ts) / S_L - et_t / S_L;	 // S_10
	ans[12] = q_ts - q_sl;						 //   -e_s;					// S_s
	//Additional states
	ans[13] = forcing_values[0] * c_1;
	ans[14] = q_pl; // q_pl;
	ans[15] = q_sl * A_h - q_b * 60.0;
	for (i = 0; i < num_parents; i++)
		ans[15] += y_p[i * dim + 15] * 60.0;
	//ans[6] += k_3*y_p[i].ve[3]*A_h;
	ans[15] *= v_B / L;
	//unsigned short i;
	//
	//// Model parameters (spatially uniform)
	//double lambda_1 = global_params[1];
	//double k_3 = global_params[4];	//[1/min]
	//double h_b = global_params[6];	//[m]
	//double S_L = global_params[7];	//[m]
	//double A = global_params[8];
	//double B = global_params[9];
	//double exponent = global_params[10];
	//double v_B = global_params[11];
	//double psi_sat = global_params[12];
	//
	//// Model parameters (spatially variable)
	//double L = params[1];	//[m]
	//double A_h = params[2];	//[m^2]
	//double K_sat = params[3];
	//double invtau = params[4];	//[1/min]
	//double k_2 = params[5];	//[1/min]
	//double k_i = params[6];	//[1/min]
	//double c_1 = params[7];
	//double c_2 = params[8];
	//double theta_s = params[9];
	//double theta_r = params[10];
	//double bc_lambda = params[11];
	//
	//double L_Top = S_L * 10;
	//
	//// Initial conditions (or from last iteration)
	//double q = y_i[0];		//[m^3/s]
	//double s_p = y_i[1];	//[m]
	//double s_1 = y_i[2];	//[m]
	//double s_2 = y_i[3];	//[m]
	//double s_3 = y_i[4];
	//double s_4 = y_i[5];
	//double s_5 = y_i[6];
	//double s_6 = y_i[7];
	//double s_7 = y_i[8];
	//double s_8 = y_i[9];
	//double s_9 = y_i[10];
	//double s_10 = y_i[11];
	//double s_s = y_i[12];
	//
	//double q_b = y_i[15];	//[m^3/s]
	//
	//						// Potential Evapotranspiration
	//double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]
	//																//Evaporation
	//double e_p, e_t, e_s;
	//double avg_s_t = (s_1 + s_2 + s_3 + s_4 + s_5 + s_6 + s_7 + s_8 + s_9 + s_10) / 10;
	//double Corr = s_p + avg_s_t + s_s / (h_b - L_Top);
	//if (e_pot > 0.0 && Corr > 1e-12)
	//{
	//	e_p = s_p * e_pot / Corr;
	//	e_t = avg_s_t * e_pot / Corr;
	//	e_s = s_s / (h_b - L_Top) * e_pot / Corr;
	//}
	//else
	//{
	//	e_p = 0.0;
	//	e_t = 0.0;
	//	e_s = 0.0;
	//}
	////Fluxes
	//
	//double pow_term = (1.0 - s_1 / S_L > 0.0) ? pow(1.0 - s_1 / S_L, exponent) : 0.0;
	//double k_t = (A + B * pow_term) * k_2;
	//
	//double q_pl = k_2 * s_p;
	//double q_pt = k_t * s_p;
	//double q_rain = forcing_values[0] * c_1; // m/min
	//double q_excess;
	//double q_2 = flux(s_1, s_2, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);	 // [m/min]
	//double q_3 = flux(s_2, s_3, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);	 // [m/min]
	//double q_4 = flux(s_3, s_4, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);	 // [m/min]
	//double q_5 = flux(s_4, s_5, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);	 // [m/min]
	//double q_6 = flux(s_5, s_6, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);	 // [m/min]
	//double q_7 = flux(s_6, s_7, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);	 // [m/min]
	//double q_8 = flux(s_7, s_8, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);	 // [m/min]
	//double q_9 = flux(s_8, s_9, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);	 // [m/min]
	//double q_10 = flux(s_9, s_10, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L); // [m/min]
	//																				 //s_10 = (4 * s_9 - s_8) / 3;
	//																				 //double q_ts = k_i * s_10 * S_L; // K_i * s_10 [m/min]
	//double q_ts = s_10 * S_L * (1 - (s_s / (h_b - S_L * 10))) * k_3;
	//double q_sl = k_3 * s_s;	//[m/min]
	//
	//double dsp, ds1;
	//if (s_1 >= theta_s && q_pt >= q_2) {
	//	q_excess = q_pt - q_2;
	//	s_1 = theta_s;
	//	ds1 = 0;
	//}
	//else {
	//	q_excess = 0;
	//	ds1 = (q_pt - q_2) / S_L;	// S_1
	//}
	//// dx/dt ODE definition for solver (RK45)
	////Discharge
	//ans[0] = -q + (q_pl + q_sl) * c_2;
	//for (i = 0; i<num_parents; i++)
	//	ans[0] += y_p[i * dim];
	//ans[0] = invtau * pow(q, lambda_1) * ans[0];
	//
	////Hillslope
	//ans[1] = q_rain + q_excess - q_pt - e_p; // S_p
	//ans[2] = ds1;							// S_1
	//ans[3] = (q_2 - q_3) / S_L;				// S_2
	//ans[4] = (q_3 - q_4) / S_L;				// S_3
	//ans[5] = (q_4 - q_5) / S_L;				// S_4
	//ans[6] = (q_5 - q_6) / S_L;				// S_5
	//ans[7] = (q_6 - q_7) / S_L;				// S_6
	//ans[8] = (q_7 - q_8) / S_L;				// S_7
	//ans[9] = (q_8 - q_9) / S_L;				// S_8
	//ans[10] = (q_9 - q_10) / S_L; 			// S_9
	//										// analytical (4 * s_9 - s_8) /3		// S_10
	//ans[11] = (q_10 - q_ts) / S_L;							// S_10
	//ans[12] = q_ts - q_sl;				// S_s
	//										//Additional states
	//ans[13] = forcing_values[0] * c_1;
	//ans[14] = q_pl;// q_pl;
	//ans[15] = q_sl * A_h - q_b*60.0;
	//for (i = 0; i<num_parents; i++)
	//	ans[15] += y_p[i * dim + 15] * 60.0;
	////ans[6] += k_3*y_p[i].ve[3]*A_h;
	//ans[15] *= v_B / L;
}

float flux_ts(float theta_1, float theta_2, float K_sat, float psi_sat, float lambda, float theta_s, float theta_r, float S_L)
{
	float theta_hat = (0.025 * theta_1 + 0.25 * theta_2) / 0.275;
	float phi_1 = (theta_1 - theta_r) / (theta_s - theta_r);
	float phi_2 = (theta_2 - theta_r) / (theta_s - theta_r);

	// K
	float phi_k = (theta_hat - theta_r) / (theta_s - theta_r);
	float power_k = (2 + 3 * lambda) / lambda;
	float K_phi_hat = K_sat * pow(phi_k, power_k);

	// psi
	float power_psi = -1 / lambda;
	float psi_1 = psi_sat * pow(phi_1, power_psi);
	float psi_2 = psi_sat * pow(phi_2, power_psi);

	return K_phi_hat * (1 - (psi_2 - psi_1) / 0.275);
}
// MODEL UID = 1006
//Type 8 based on the paper Herrada Martin 2014 JHydrol http://dx.doi.org/10.1016/j.jhydrol.2014.04.026
void navid(double t, const double *const y_i, unsigned int dim, const double *const y_p, unsigned short num_parents, unsigned int max_dim, const double *const global_params, const double *const params, const double *const forcing_values, const QVSData *const qvs, int state, void *user, double *ans, double h)
{
	unsigned short i;

	// Model parameters (spatially uniform)
	double lambda_1 = global_params[1];
	double k_3 = global_params[4]; //[1/min]
	double h_b = global_params[6]; //[m]
	double S_L = global_params[7]; //[m]
	double A = global_params[8];
	double B = global_params[9];
	double exponent = global_params[10];
	double v_B = global_params[11];
	double psi_sat = global_params[12];

	// Model parameters (spatially variable)
	double L = params[1];	//[m]
	double A_h = params[2]; //[m^2]
	double K_sat = params[3];
	double invtau = params[4]; //[1/min]
	double k_2 = params[5];	   //[1/min]
	double k_i = params[6];	   //[1/min]
	double c_1 = params[7];
	double c_2 = params[8];
	double theta_s = params[9];
	double theta_r = params[10];
	double bc_lambda = params[11];

	double L_Top = S_L * 10;

	// Initial conditions (or from last iteration)
	double q = y_i[0];	 //[m^3/s]
	double s_p = y_i[1]; //[m]
	double s_1 = y_i[2]; //[m]
	double s_2 = y_i[3]; //[m]
	double s_3 = y_i[4];
	double s_4 = y_i[5];
	double s_5 = y_i[6];
	double s_6 = y_i[7];
	double s_7 = y_i[8];
	double s_8 = y_i[9];
	double s_9 = y_i[10];
	double s_10 = y_i[11];
	double s_s = y_i[12];
	//s_s = (4 * s_10 - s_9)/3;
	double q_b = y_i[15]; //[m^3/s]
	printf("test %f %f", K_sat, params[3]);
	// Potential Evapotranspiration
	double e_pot = forcing_values[1] * (1e-3 / (30.0 * 24.0 * 60.0)); //[mm/month] -> [m/min]
																	  //Evaporation
	double e_p, e_t, e_s;
	/* double avg_s_t = (s_1 + s_2 + s_3 + s_4 + s_5 + s_6 + s_7 + s_8 + s_9 + s_10) / 10;
	double Corr = s_p + avg_s_t + s_s / (h_b - L_Top); */

	// if (e_pot > 0.0 && Corr > 1e-12)
	// {
	// e_p = s_p * e_pot / Corr;
	// e_t = avg_s_t * e_pot / Corr;
	// e_s = s_s / (h_b - L_Top) * e_pot / Corr;
	// }
	// else
	// {
	e_s = 0.0;
	e_p = 0.0;
	e_t = 0.0;
	// }
	double pow_term = (1.0 - s_1 / theta_s > 0.0) ? pow(1.0 - s_1 / theta_s, exponent) : 0.0;
	double k_t = (A + B * pow_term) * k_2;

	//Fluxes
	double q_pl = k_2 * s_p;
	double q_pt = k_t * s_p;
	double q_rain = forcing_values[0] * c_1; // m/min
	double q_2 = flux(s_1, s_2, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_3 = flux(s_2, s_3, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_4 = flux(s_3, s_4, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_5 = flux(s_4, s_5, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_6 = flux(s_5, s_6, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_7 = flux(s_6, s_7, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_8 = flux(s_7, s_8, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_9 = flux(s_8, s_9, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_10 = flux(s_9, s_10, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	// s_10 = (4 * s_9 - s_8) / 3;
	//double q_ts = k_i * s_10 * S_L;
	double q_sl = 0; // k_3 * s_s * (h_b - L_Top);	//[m/min]

	// This is second option !!!!!!!!!!
	double q_ts = flux_ts(s_10, s_s, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);

	// dx/dt ODE definition for solver (RK45)
	//Discharge
	ans[0] = -q + (q_pl + q_sl) * c_2;
	for (i = 0; i < num_parents; i++)
		ans[0] += y_p[i * dim];
	ans[0] = invtau * pow(q, lambda_1) * ans[0];

	//Hillslope
	ans[1] = q_rain - q_pl - q_pt; // S_P
	ans[2] = (q_pt - q_2) / S_L;   // S_1
	ans[3] = (q_2 - q_3) / S_L;	   // S_2
	ans[4] = (q_3 - q_4) / S_L;	   // S_3
	ans[5] = (q_4 - q_5) / S_L;	   // S_4
	ans[6] = (q_5 - q_6) / S_L;	   // S_5
	ans[7] = (q_6 - q_7) / S_L;	   // S_6
	ans[8] = (q_7 - q_8) / S_L;	   // S_7
	ans[9] = (q_8 - q_9) / S_L;	   // S_8
	ans[10] = (q_9 - q_10) / S_L;  // S_9
								   // S_10 : analytical (4 * s_9 - s_8) /3
	ans[11] = (4 * ans[10] - ans[9]) / 3;
	ans[12] = 0; // (q_ts * 3 - q_sl) / (h_b - L_Top);
	//S_s
	//Additional states
	ans[13] = forcing_values[0] * c_1;
	ans[14] = q_pl; // q_pl;
	ans[15] = q_sl * A_h - q_b * 60.0;
	for (i = 0; i < num_parents; i++)
		ans[15] += y_p[i * dim + 15] * 60.0;
	ans[15] *= v_B / L;
}

// mode UID 1007
void navid0(double t, const double *const y_i, unsigned int dim, const double *const y_p, unsigned short num_parents, unsigned int max_dim, const double *const global_params, const double *const params, const double *const forcing_values, const QVSData *const qvs, int state, void *user, double *ans, double h)
{
	unsigned short i;

	// Model parameters (spatially uniform)
	double lambda_1 = global_params[1];
	double k_3 = global_params[4]; //[1/min]
	double h_b = global_params[6]; //[m]
	double S_L = global_params[7]; //[m]
	double A = global_params[8];
	double B = global_params[9];
	double exponent = global_params[10];
	double v_B = global_params[11];
	double psi_sat = global_params[12];

	// Model parameters (spatially variable)
	double L = params[1];	//[m]
	double A_h = params[2]; //[m^2]
	double K_sat = params[3];
	double invtau = params[4]; //[1/min]
	double k_2 = params[5];	   //[1/min]
	double k_i = params[6];	   //[1/min]
	double c_1 = params[7];
	double c_2 = params[8];
	double theta_s = params[9];
	double theta_r = params[10];
	double bc_lambda = params[11];

	double L_Top = S_L * 10;

	// Initial conditions (or from last iteration)
	double q = y_i[0];	 //[m^3/s]
	double s_p = y_i[1]; //[m]
	double s_1 = y_i[2]; //[m]
	double s_2 = y_i[3]; //[m]
	double s_3 = y_i[4];
	double s_4 = y_i[5];
	double s_5 = y_i[6];
	double s_6 = y_i[7];
	double s_7 = y_i[8];
	double s_8 = y_i[9];
	double s_9 = y_i[10];
	double s_10 = y_i[11];
	double s_s = y_i[12];
	//s_s = (4 * s_10 - s_9)/3;
	double q_b = y_i[15]; //[m^3/s]

	// Potential Evapotranspiration
	double e_pot = forcing_values[1] * (1e-3 / (30.0 * 24.0 * 60.0)); //[mm/month] -> [m/min]
																	  //Evaporation
	double e_p, e_t, e_s;
	//double avg_s_t = (s_1 + s_2 + s_3 + s_4 + s_5 + s_6 + s_7 + s_8 + s_9 + s_10) / 10;
	double Corr = s_p + (s_1 - theta_r) / (theta_s - theta_r) + s_s / (h_b - L_Top);

	if (e_pot > 0.0 && Corr > 1e-12)
	{
		e_p = s_p * e_pot / Corr;
		e_t = s_1 > theta_r + 0.03 ? (s_1 - theta_r) / (theta_s - theta_r) * (e_pot) / Corr : 0.0;
		e_s = s_s / (h_b - L_Top) * e_pot / Corr;
	}
	else
	{
		e_s = 0.0;
		e_p = 0.0;
		e_t = 0.0;
	}
	double pow_term = (1.0 - (s_1 / theta_s) > 0.03) ? pow(1.0 - (s_1 - theta_r) / (theta_s - theta_r), exponent) : 0.0;
	double k_t = (A + B * pow_term) * k_2;

	//Fluxes
	double q_pl = k_2 * s_p;
	double q_pt = k_t * s_p;
	double q_rain = forcing_values[0] * c_1; // m/min
	double q_2 = flux(s_1, s_2, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_3 = flux(s_2, s_3, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_4 = flux(s_3, s_4, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_5 = flux(s_4, s_5, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_6 = flux(s_5, s_6, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_7 = flux(s_6, s_7, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_8 = flux(s_7, s_8, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_9 = flux(s_8, s_9, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_10 = flux(s_9, s_10, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	// s_10 = (4 * s_9 - s_8) / 3;
	//double q_ts = k_i * s_10 * S_L;
	double q_sl = 0; // k_3 * s_s * (h_b - L_Top);	//[m/min]

	// This is second option !!!!!!!!!!
	//double q_ts = flux_ts(s_10, s_s, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);

	// dx/dt ODE definition for solver (RK45)
	//Discharge
	ans[0] = -q + (q_pl + q_sl) * c_2;
	for (i = 0; i < num_parents; i++)
		ans[0] += y_p[i * dim];
	ans[0] = invtau * pow(q, lambda_1) * ans[0];

	//Hillslope
	ans[1] = q_rain - q_pl - q_pt - e_p; // S_P
	ans[2] = (q_pt - q_2) / S_L;		 // -e_t / S_L;		// S_1
	ans[3] = (q_2 - q_3) / S_L;			 // S_2
	ans[4] = (q_3 - q_4) / S_L;			 // S_3
	ans[5] = (q_4 - q_5) / S_L;			 // S_4
	ans[6] = (q_5 - q_6) / S_L;			 // S_5
	ans[7] = (q_6 - q_7) / S_L;			 // S_6
	ans[8] = (q_7 - q_8) / S_L;			 // S_7
	ans[9] = (q_8 - q_9) / S_L;			 // S_8
	ans[10] = (q_9 - q_10) / S_L;		 // S_9
										 // S_10 : analytical (4 * s_9 - s_8) /3
	ans[11] = (4 * ans[10] - ans[9]) / 3;
	ans[12] = 0; // (q_ts * 3 - q_sl) / (h_b - L_Top);
				 //S_s
				 //Additional states
	ans[13] = forcing_values[0] * c_1;
	ans[14] = q_pl; // q_pl;
	ans[15] = q_sl * A_h - q_b * 60.0;
	for (i = 0; i < num_parents; i++)
		ans[15] += y_p[i * dim + 15] * 60.0;
	ans[15] *= v_B / L;
}

// mode UID 1008
void navid1(double t, const double *const y_i, unsigned int dim, const double *const y_p, unsigned short num_parents, unsigned int max_dim, const double *const global_params, const double *const params, const double *const forcing_values, const QVSData *const qvs, int state, void *user, double *ans, double h)
{
	unsigned short i;

	// Model parameters (spatially uniform)
	double lambda_1 = global_params[1];
	double k_3 = global_params[4]; //[1/min]
	double h_b = global_params[6]; //[m]
	double S_L = global_params[7]; //[m]
	double A = global_params[8];
	double B = global_params[9];
	double exponent = global_params[10];
	double v_B = global_params[11];
	double psi_sat = global_params[12];

	// Model parameters (spatially variable)
	double L = params[1];	//[m]
	double A_h = params[2]; //[m^2]
	double K_sat = params[3];
	double invtau = params[4]; //[1/min]
	double k_2 = params[5];	   //[1/min]
	double k_i = params[6];	   //[1/min]
	double c_1 = params[7];
	double c_2 = params[8];
	double theta_s = params[9];
	double theta_r = params[10];
	double bc_lambda = params[11];

	double L_Top = S_L * 10;

	// Initial conditions (or from last iteration)
	double q = y_i[0];	 //[m^3/s]
	double s_p = y_i[1]; //[m]
	double s_1 = y_i[2]; //[m]
	double s_2 = y_i[3]; //[m]
	double s_3 = y_i[4];
	double s_4 = y_i[5];
	double s_5 = y_i[6];
	double s_6 = y_i[7];
	double s_7 = y_i[8];
	double s_8 = y_i[9];
	double s_9 = y_i[10];
	double s_10 = y_i[11];
	double s_s = y_i[12];
	//s_s = (4 * s_10 - s_9)/3;
	double q_b = y_i[15]; //[m^3/s]

	// Potential Evapotranspiration
	double e_pot = forcing_values[1] * (1e-3 / (30.0 * 24.0 * 60.0)); //[mm/month] -> [m/min]
																	  //Evaporation

	double pow_term = (1.0 - (s_1 - theta_r) / (theta_s - theta_r) > 0.01) ? pow(1.0 - (s_1 - theta_r) / (theta_s - theta_r), exponent) : 0.0;
	double k_t = (A + B * pow_term) * k_2;

	//Fluxes
	double q_pl = k_2 * s_p;
	double q_pt = k_t * s_p;
	double q_rain = forcing_values[0] * c_1; // m/min
	double q_2 = flux(s_1, s_2, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_3 = flux(s_2, s_3, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_4 = flux(s_3, s_4, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_5 = flux(s_4, s_5, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_6 = flux(s_5, s_6, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_7 = flux(s_6, s_7, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_8 = flux(s_7, s_8, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_9 = flux(s_8, s_9, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_10 = flux(s_9, s_10, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	// s_10 = (4 * s_9 - s_8) / 3;
	//double q_ts = k_i * s_10 * S_L;
	double q_sl = 0.0; // k_3 * s_s * (h_b - L_Top);	//[m/min]
	double e_p, e_t;   // 1, e_t2, e_t3, e_t4, e_t5, e_t6, e_t7, e_t8, e_t9, e_t10, e_s;
		//double avg_s_t = (s_1 + s_2 + s_3 + s_4 + s_5 + s_6 + s_7 + s_8 + s_9 + s_10) / 10;
	double Corr = s_p + (s_1 - theta_r) / (theta_s - theta_r); // +s_s / (h_b - L_Top);

	if (e_pot > 0.0 && Corr > 1e-2)
	{
		e_p = s_p * e_pot / Corr;
		e_t = s_1 > theta_r + q_2 / S_L ? (s_1 - theta_r) * (e_pot / Corr) : 0.0;
		//e_t5 = s_5  > theta_r + 0.03 ? (s_1 - theta_r) / (theta_s - theta_r) * (e_pot / 10) / Corr : 0.0;
		//e_t6 = s_6  > theta_r + 0.03 ? (s_1 - theta_r) / (theta_s - theta_r) * (e_pot / 10) / Corr : 0.0;
		//e_t7 = s_7  > theta_r + 0.03 ? (s_1 - theta_r) / (theta_s - theta_r) * (e_pot / 10) / Corr : 0.0;
		//e_t8 = s_8  > theta_r + 0.03 ? (s_1 - theta_r) / (theta_s - theta_r) * (e_pot / 10) / Corr : 0.0;
		//e_t9 = s_9  > theta_r + 0.03 ? (s_1 - theta_r) / (theta_s - theta_r) * (e_pot / 10) / Corr : 0.0;
		//e_t10 = s_10  > theta_r + 0.03 ? (s_1 - theta_r) / (theta_s - theta_r) * (e_pot / 10) / Corr : 0.0;
		//e_s = s_s / (h_b - L_Top) * e_pot / Corr;
	}
	else
	{
		//e_s = 0.0;
		e_p = 0.0;
		e_t = 0.0;
	}

	// This is second option !!!!!!!!!!
	//double q_ts = flux_ts(s_10, s_s, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);

	// dx/dt ODE definition for solver (RK45)
	//Discharge
	ans[0] = -q + (q_pl + q_sl) * c_2;
	for (i = 0; i < num_parents; i++)
		ans[0] += y_p[i * dim];
	ans[0] = invtau * pow(q, lambda_1) * ans[0];

	//Hillslope
	ans[1] = q_rain - q_pl - q_pt - e_p;	 // S_P
	ans[2] = (q_pt - q_2) / S_L - e_t / S_L; // S_1
	ans[3] = (q_2 - q_3) / S_L;				 // S_2
	ans[4] = (q_3 - q_4) / S_L;				 // S_3
	ans[5] = (q_4 - q_5) / S_L;				 // S_4
	ans[6] = (q_5 - q_6) / S_L;				 // S_5
	ans[7] = (q_6 - q_7) / S_L;				 // S_6
	ans[8] = (q_7 - q_8) / S_L;				 // S_7
	ans[9] = (q_8 - q_9) / S_L;				 // S_8
	ans[10] = (q_9 - q_10) / S_L;			 // S_9
											 // S_10 : analytical (4 * s_9 - s_8) /3
	ans[11] = (4 * ans[10] - ans[9]) / 3;
	ans[12] = 0; // (q_ts * 3 - q_sl) / (h_b - L_Top);
				 //S_s
				 //Additional states
	ans[13] = forcing_values[0] * c_1;
	ans[14] = q_pl; // q_pl;
	ans[15] = q_sl * A_h - q_b * 60.0;
	for (i = 0; i < num_parents; i++)
		ans[15] += y_p[i * dim + 15] * 60.0;
	ans[15] *= v_B / L;
}

bool is_nan(double x) { return x != x; };
// mode UID 1009
void navid2(double t, const double *const y_i, unsigned int dim, const double *const y_p, unsigned short num_parents, unsigned int max_dim, const double *const global_params, const double *const params, const double *const forcing_values, const QVSData *const qvs, int state, void *user, double *ans, double h)
{
	unsigned short i;

	// Model parameters (spatially uniform)
	double lambda_1 = global_params[1];
	double k_3 = global_params[4]; //[1/min]
	double h_b = global_params[6]; //[m]
	double S_L = global_params[7]; //[m]
	double A = global_params[8];
	double B = global_params[9];
	double exponent = global_params[10];
	double v_B = global_params[11];
	double psi_sat = global_params[12];

	// Model parameters (spatially variable)
	double L = params[1];	//[m]
	double A_h = params[2]; //[m^2]
	double K_sat = params[3];
	double invtau = params[4]; //[1/min]
	double k_2 = params[5];	   //[1/min]
	double k_i = params[6];	   //[1/min]
	double c_1 = params[7];
	double c_2 = params[8];
	double theta_s = params[9];
	double theta_r = params[10];
	double bc_lambda = params[11];

	double L_Top = S_L * 10.0;

	// Initial conditions (or from last iteration)
	double q = y_i[0];	 //[m^3/s]
	double s_p = y_i[1]; //[m]
	double s_1 = y_i[2]; //[m]
	double s_2 = y_i[3]; //[m]
	double s_3 = y_i[4];
	double s_4 = y_i[5];
	double s_5 = y_i[6];
	double s_6 = y_i[7];
	double s_7 = y_i[8];
	double s_8 = y_i[9];
	double s_9 = y_i[10];
	double s_10 = y_i[11];
	double s_s = y_i[12];
	//s_s = (4 * s_10 - s_9)/3;
	double q_b = y_i[15];											  //[m^3/s]
	double e_pot = forcing_values[1] * (1e-3 / (30.0 * 24.0 * 60.0)); //[mm/month] -> [m/min]
	double q_rain = forcing_values[0] * c_1;						  // m/min
	double q_2 = flux(s_1, s_2, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_3 = flux(s_2, s_3, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_4 = flux(s_3, s_4, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_5 = flux(s_4, s_5, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_6 = flux(s_5, s_6, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_7 = flux(s_6, s_7, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_8 = flux(s_7, s_8, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_9 = flux(s_8, s_9, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_10 = flux(s_9, s_10, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_ts = flux_ts(s_10, s_s, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	// s_10 = (4 * s_9 - s_8) / 3;
	//double q_ts = k_i * s_10 * S_L;
	double q_sl = s_s > theta_r ? k_3 * (s_s - theta_r) * (h_b - L_Top) : 0.0; //[m/min]
	double q_pl = k_2 * s_p;
	double q_sp, ds_1, e_p, e_t, ds_extra; // 1, e_t2, e_t3, e_t4, e_t5, e_t6, e_t7, e_t8, e_t9, e_t10, e_s;
		//double avg_s_t = (s_1 + s_2 + s_3 + s_4 + s_5 + s_6 + s_7 + s_8 + s_9 + s_10) / 10;

	float q_inf = q_2; // flux_100(s_t[0], s_t[1], K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	if (q_rain >= q_inf)
	{
		float q_inf_new = flux(theta_s, s_1, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
		q_sp = q_rain - q_inf_new - s_p * k_2;
		ds_1 = (q_inf_new) / S_L;
	}
	else
	{
		if (s_p > 0)
		{
			q_sp = -s_p * k_2;
			ds_1 = (q_rain - q_inf) / S_L;
		}
		else
		{
			ds_1 = (q_rain - q_inf) / S_L;
			q_sp = 0.0;
		}
	}

	double Corr = s_p + (s_1 - theta_r) / (theta_s - theta_r); // +s_s / (h_b - L_Top);

	if (e_pot > 0.0 && Corr > 1e-2)
	{
		e_p = s_p * e_pot / Corr;
		//e_t = s_1  > theta_r + q_2 / S_L ? (-q_2 / (theta_s - theta_r) + (s_1 - theta_r)/(theta_s-theta_r))  * (e_pot ) : 0.0;
		e_t = s_1 > theta_r + q_2 / S_L ? (s_1 - q_2 / S_L - theta_r) / (theta_s - theta_r) * e_pot : 0.0;
	}
	else
	{
		//e_s = 0.0;
		e_p = 0.0;
		e_t = 0.0;
	}
	ans[0] = -q + (q_pl + q_sl) * c_2;
	for (i = 0; i < num_parents; i++)
		ans[0] += y_p[i * dim];
	ans[0] = invtau * pow(q, lambda_1) * ans[0];

	//Hillslope
	ans[1] = q_sp; // q_rain - q_pl - q_pt - e_p;		// S_P
	ans[2] = ds_1; // S_1
	//ans[2] = ds_1 - e_t / ((theta_s - theta_r) *S_L);		// S_1
	ans[3] = (q_2 - q_3) / S_L;				  // S_2
	ans[4] = (q_3 - q_4) / S_L;				  // S_3
	ans[5] = (q_4 - q_5) / S_L;				  // S_4
	ans[6] = (q_5 - q_6) / S_L;				  // S_5
	ans[7] = (q_6 - q_7) / S_L;				  // S_6
	ans[8] = (q_7 - q_8) / S_L;				  // S_7
	ans[9] = (q_8 - q_9) / S_L;				  // S_8
	ans[10] = (q_9 - q_10) / S_L;			  // S_9
											  // S_10 : analytical (4 * s_9 - s_8) /3
	ans[11] = (4.0 * ans[10] - ans[9]) / 3.0; // good SM
	//ans[11] = (4 * ans[10] - ans[9]) / 3.0 -q_ts / (h_b - L_Top); // not good
	ans[12] = 0.0; // (q_ts - q_sl) / (h_b - L_Top);// (q_ts * 3 - q_sl) / (h_b - L_Top);
				   //S_s
				   //Additional states
	ans[13] = forcing_values[0] * c_1;
	ans[14] = q_pl; // q_pl;
	ans[15] = q_sl * A_h - q_b * 60.0;
	for (i = 0; i < num_parents; i++)
		ans[15] += y_p[i * dim + 15] * 60.0;
	ans[15] *= v_B / L;
}
// uid 1010
void navid_gampt(double t, const double *const y_i, unsigned int dim, const double *const y_p, unsigned short num_parents, unsigned int max_dim, const double *const global_params, const double *const params, const double *const forcing_values, const QVSData *const qvs, int state, void *user, double *ans, double h)
{
	unsigned short i;

	// Model parameters (spatially uniform)
	float lambda_1 = global_params[1];
	float k_3 = global_params[4]; //[1/min]
	float h_b = global_params[6]; //[m]
	float S_L = global_params[7]; //[m]
	float A = global_params[8];
	float B = global_params[9];
	float exponent = global_params[10];
	float v_B = global_params[11];
	float psi_sat = global_params[12];

	// Model parameters (spatially variable)
	float L = params[1];   //[m]
	float A_h = params[2]; //[m^2]
	float K_sat = params[3];
	float invtau = params[4]; //[1/min]
	float k_2 = params[5];	  //[1/min]
	float k_i = params[6];	  //[1/min]
	float c_1 = params[7];
	float c_2 = params[8];
	float theta_s = (int)(params[9] * 1000) / 1000.0; // floorf(*1000) / 1000;
	float theta_r = params[10];
	float bc_lambda = params[11];

	float L_Top = S_L * 10.0;

	// Initial conditions (or from last iteration)
	float q = y_i[0];						   //[m^3/s]
	float s_p = y_i[1];						   //[m]
	float s_1 = (int)(y_i[2] * 1000) / 1000.0; //[m]
	float s_2 = y_i[3];						   //[m]
	float s_3 = y_i[4];
	float s_4 = y_i[5];
	float s_5 = y_i[6];
	float s_6 = y_i[7];
	float s_7 = y_i[8];
	float s_8 = y_i[9];
	float s_9 = y_i[10];
	float s_10 = y_i[11];
	float s_s = y_i[12];
	//s_s = (4 * s_10 - s_9)/3;
	float q_b = y_i[15];											 //[m^3/s]
	float e_pot = forcing_values[1] * (1e-3 / (30.0 * 24.0 * 60.0)); //[mm/month] -> [m/min]
	float q_rain = forcing_values[0] * c_1;							 // m/min
	float q_2 = flux(s_1, s_2, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	float q_3 = flux(s_2, s_3, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	float q_4 = flux(s_3, s_4, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	float q_5 = flux(s_4, s_5, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	float q_6 = flux(s_5, s_6, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	float q_7 = flux(s_6, s_7, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	float q_8 = flux(s_7, s_8, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	float q_9 = flux(s_8, s_9, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	float q_10 = flux(s_9, s_10, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	float q_ts = flux_ts(s_10, s_s, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	// s_10 = (4 * s_9 - s_8) / 3;
	//double q_ts = k_i * s_10 * S_L;
	float q_sl = s_s > theta_r ? k_3 * (s_s - theta_r) * (h_b - L_Top) : 0.0; //[m/min]
	float q_pl, q_sp, ds_1, e_p, e_t;										  // 1, e_t2, e_t3, e_t4, e_t5, e_t6, e_t7, e_t8, e_t9, e_t10, e_s;
																			  //double avg_s_t = (s_1 + s_2 + s_3 + s_4 + s_5 + s_6 + s_7 + s_8 + s_9 + s_10) / 10;
	q_pl = s_p * k_2;
	if (s_1 >= theta_s && q_rain >= q_2)
	{
		//ds_extra = (s_1 - theta_s) * S_L;
		//s_1 = floorf(theta_s * 100) / 100;
		s_1 = theta_s;

		q_sp = q_rain - q_2; // +ds_extra;
		ds_1 = 0.0;
		//q_2 = flux(s_1, s_2, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	}
	else
	{
		//q_2 = flux(s_1, s_2, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
		ds_1 = (q_rain - q_2) / S_L;
		q_sp = 0.0;
	}

	float Corr = s_p + (s_1 - theta_r) / (theta_s - theta_r); // +s_s / (h_b - L_Top);

	if (e_pot > 0.0 && Corr > 1e-2)
	{
		e_p = s_p * e_pot / Corr;
		//e_t = s_1  > theta_r + q_2 / S_L ? (-q_2 / (theta_s - theta_r) + (s_1 - theta_r)/(theta_s-theta_r))  * (e_pot ) : 0.0;
		e_t = s_1 > theta_r + q_2 / S_L ? (s_1 - q_2 / S_L - theta_r) * e_pot : 0.0;
	}
	else
	{
		//e_s = 0.0;
		e_p = 0.0;
		e_t = 0.0;
	}
	ans[0] = -q + (q_pl + q_sl) * c_2;
	for (i = 0; i < num_parents; i++)
		ans[0] += y_p[i * dim];
	ans[0] = invtau * pow(q, lambda_1) * ans[0];

	//Hillslope
	ans[1] = q_sp - q_pl;		  // q_rain - q_pl - q_pt - e_p;		// S_P
	ans[2] = ds_1 - e_t / S_L;	  // S_1
								  //ans[2] = ds_1 - e_t / ((theta_s - theta_r) *S_L);		// S_1
	ans[3] = (q_2 - q_3) / S_L;	  // S_2
	ans[4] = (q_3 - q_4) / S_L;	  // S_3
	ans[5] = (q_4 - q_5) / S_L;	  // S_4
	ans[6] = (q_5 - q_6) / S_L;	  // S_5
	ans[7] = (q_6 - q_7) / S_L;	  // S_6
	ans[8] = (q_7 - q_8) / S_L;	  // S_7
	ans[9] = (q_8 - q_9) / S_L;	  // S_8
	ans[10] = (q_9 - q_10) / S_L; // S_9
								  // S_10 : analytical (4 * s_9 - s_8) /3
	ans[11] = (q_10 - q_ts) / S_L;
	(4.0 * ans[10] - ans[9]) / 3.0; // good SM
									//ans[11] = (4 * ans[10] - ans[9]) / 3.0 -q_ts / (h_b - L_Top); // not good
	ans[12] = 0.0;					// (q_ts - q_sl) / (h_b - L_Top);// (q_ts * 3 - q_sl) / (h_b - L_Top);
									//S_s
									//Additional states
	ans[13] = forcing_values[0] * c_1;
	ans[14] = q_pl; // q_pl;
	ans[15] = q_sl * A_h - q_b * 60.0;
	for (i = 0; i < num_parents; i++)
		ans[15] += y_p[i * dim + 15] * 60.0;
	ans[15] *= v_B / L;
}

float flux_new(float theta_1, float theta_2, float K_sat, float psi_sat, float lambda, float theta_s, float theta_r, float S_L_1, float S_L_2)
{
	//theta_1 = theta_1 <= theta_r ? theta_r + 0.01 : theta_1;
	//theta_2 = theta_2 <= theta_r ? theta_r + 0.01 : theta_2;
	//theta_1 = theta_1 > theta_s ? theta_s : theta_1;
	//theta_2 = theta_2 > theta_s ? theta_s : theta_2;
	float delta_z = (S_L_1 + S_L_2) / 2.0;
	float theta_hat = (theta_2 * S_L_2 + theta_1 * S_L_1) / (S_L_1 + S_L_2);
	float phi_1 = (theta_1 - theta_r) / (theta_s - theta_r);
	float phi_2 = (theta_2 - theta_r) / (theta_s - theta_r);

	// K
	float phi_k = (theta_hat - theta_r) / (theta_s - theta_r);
	float power_k = (2.0 + 3.0 * lambda) / lambda;
	float K_phi_hat = K_sat * pow(phi_k, power_k);

	// psi
	float power_psi = -1.0 / lambda;
	float psi_1 = psi_sat * pow(phi_1, power_psi);
	float psi_2 = psi_sat * pow(phi_2, power_psi);

	return K_phi_hat * (1.0 - (psi_2 - psi_1) / delta_z);
}
// model UID 1012
void navid_gampt1012(double t, const double *const y_i, unsigned int dim, const double *const y_p, unsigned short num_parents, unsigned int max_dim, const double *const global_params, const double *const params, const double *const forcing_values, const QVSData *const qvs, int state, void *user, double *ans, double h)
{
	unsigned short i;

	// Model parameters (spatially uniform)
	float lambda_1 = global_params[1];
	float k_3 = global_params[4]; //[1/min]
	float h_b = global_params[6]; //[m]
	//float S_L = global_params[7];	//[m]
	float A = global_params[8];
	float B = global_params[9];
	float exponent = global_params[10];
	float v_B = global_params[11];
	float psi_sat = global_params[12];

	// Model parameters (spatially variable)
	float L = params[1];   //[m]
	float A_h = params[2]; //[m^2]
	float K_sat = params[3];
	float invtau = params[4]; //[1/min]
	float k_2 = params[5];	  //[1/min]
	float k_i = params[6];	  //[1/min]
	float c_1 = params[7];
	float c_2 = params[8];
	float theta_s = (int)(params[9] * 1000) / 1000.0; // floorf(*1000) / 1000;
	float theta_r = params[10];
	float bc_lambda = params[11];
	float s_t[10];
	float q_t[9];
	float S_L[10] = {0.05, 0.05, 0.05, 0.05, 0.1, 0.7, 1.0, 1.0, 1.0, 1.0};
	float L_Top = 5;

	// Initial conditions (or from last iteration)
	float q = y_i[0];	//[m^3/s]
	float s_p = y_i[1]; //[m]

	for (i = 0; i < 10; i++)
		s_t[i] = y_i[i + 2];
	float s_s = y_i[12];
	float q_b = y_i[15];											 //[m^3/s]
	float e_pot = forcing_values[1] * (1e-3 / (30.0 * 24.0 * 60.0)); //[mm/month] -> [m/min]
	float q_rain = forcing_values[0] * c_1;							 // m/min

	for (i = 0; i < 9; i++)
		q_t[i] = flux_new(s_t[i], s_t[i + 1], K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L[i], S_L[i + 1]);
	float q_sl = 0;				// s_s > theta_r ? k_3 * (s_s - theta_r)  * (h_b - L_Top) : 0.0;	//[m/min]
	float q_pl, ds_1, e_p, e_t; //

	q_pl = k_2 * s_p;
	float Corr = s_p + (s_t[0] - theta_r) / (theta_s - theta_r); // +s_s / (h_b - L_Top
	if (e_pot > 0.0 && Corr > 1e-2)
	{
		e_p = s_p * e_pot / Corr;
		e_t = s_t[0] > theta_r + q_t[0] / S_L[0] ? (s_t[0] - q_t[0] / S_L[0] - theta_r) * e_pot : 0.0;
	}
	else
	{
		e_p = 0.0;
		e_t = 0.0;
	}
	ans[0] = -q + (q_pl + q_sl) * c_2;
	for (i = 0; i < num_parents; i++)
		ans[0] += y_p[i * dim];
	ans[0] = invtau * pow(q, lambda_1) * ans[0];

	//Hillslope
	ans[1] = q_rain - q_pl - q_t[0] - e_p;	 // S_P
	ans[2] = q_t[0] - q_t[1] - e_t / S_L[0]; // S_t[0]
	for (i = 0; i < 8; i++)
		ans[i + 3] = (q_t[i + 1] - q_t[i + 2]) / S_L[i + 1]; // S_t[1:8]
	ans[11] = (4.0 * ans[10] - ans[9]) / 3.0;				 // S_t[9] good SM
	ans[12] = 0.0;											 // S_s
															 //Additional states
	ans[13] = forcing_values[0] * c_1;
	ans[14] = q_pl; // q_pl;
	ans[15] = q_sl * A_h - q_b * 60.0;
	for (i = 0; i < num_parents; i++)
		ans[15] += y_p[i * dim + 15] * 60.0;
	ans[15] *= v_B / L;
}

void herrada_gampt1013(double t, const double *const y_i, unsigned int dim, const double *const y_p, unsigned short num_parents, unsigned int max_dim, const double *const global_params, const double *const params, const double *const forcing_values, const QVSData *const qvs, int state, void *user, double *ans, double h)
{
	unsigned short i;

	// Model parameters (spatially uniform)
	float lambda_1 = global_params[1];
	float k_3 = global_params[4]; //[1/min]
	float h_b = global_params[6]; //[m]
	//float S_L = global_params[7];	//[m]
	float A = global_params[8];
	float B = global_params[9];
	float exponent = global_params[10];
	float v_B = global_params[11];
	float psi_sat = global_params[12];

	// Model parameters (spatially variable)
	float L = params[1];   //[m]
	float A_h = params[2]; //[m^2]
	float K_sat = params[3];
	float invtau = params[4]; //[1/min]
	float k_2 = params[5];	  //[1/min]
	float k_i = params[6];	  //[1/min]
	float c_1 = params[7];
	float c_2 = params[8];
	float theta_s = (int)(params[9] * 1000) / 1000.0; // floorf(*1000) / 1000;
	float theta_r = params[10];
	float bc_lambda = params[11];
	float s_t[10];
	float q_t[9];
	float S_L[10] = {0.05, 0.05, 0.05, 0.05, 0.1, 0.7, 1.0, 1.0, 1.0, 1.0};
	//float S_L[10] = { 0.05, 0.05,0.05, 0.05,0.05, 0.05,0.05, 0.05,0.05, 0.05 };// 0.05, 0.05, 0.1, 0.7, 1.0, 1.0, 1.0, 1.0
	float L_Top = 5;

	// Initial conditions (or from last iteration)
	float q = y_i[0];	//[m^3/s]
	float s_p = y_i[1]; //[m]

	for (i = 0; i < 10; i++)
		s_t[i] = y_i[i + 2];
	float s_s = y_i[12];
	float q_b = y_i[15]; //[m^3/s]
	for (i = 0; i < 9; i++)
		q_t[i] = flux_new(s_t[i], s_t[i + 1], K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L[i], S_L[i + 1]);
	float w_front = y_i[16];
	float e_pot = forcing_values[1] * (1e-3 / (30.0 * 24.0 * 60.0));		  //[mm/month] -> [m/min]
	float q_rain = forcing_values[0] * c_1;									  // m/min
	float q_sl = s_s > theta_r ? k_3 * (s_s - theta_r) * (h_b - L_Top) : 0.0; //[m/min]
	float q_pl, ds_1, e_p, e_t;
	float q_infilt = K_sat * (s_p + w_front) / w_front;
	float q_w_front = K_sat * (s_p + w_front) / (w_front * theta_s * (1 - s_t[0]));
	q_pl = k_2 * s_p;
	// Evapotranspiration correction for top layer and ponding states
	float Corr = s_p + (s_t[0] - theta_r) / (theta_s - theta_r); // +s_s / (h_b - L_Top);
	if (e_pot > 0.0 && Corr > 1e-2)
	{
		e_p = s_p * e_pot / Corr;
		e_t = s_t[0] > theta_r + q_t[0] / S_L[0] ? (s_t[0] - q_t[0] / S_L[0] - theta_r) * e_pot : 0.0;
	}
	else
	{
		e_p = 0.0;
		e_t = 0.0;
	}

	ans[0] = -q + (q_pl + q_sl) * c_2;
	for (i = 0; i < num_parents; i++)
		ans[0] += y_p[i * dim];
	ans[0] = invtau * pow(q, lambda_1) * ans[0];

	//Hillslope
	ans[1] = q_rain - q_infilt - q_pl;			 // S_P
	ans[2] = (q_infilt - q_t[0] - e_t) / S_L[0]; // S_t[0]
	for (i = 0; i < 8; i++)
		ans[i + 3] = (q_t[i] - q_t[i + 1]) / S_L[i + 1]; // S_t[1:8]
	ans[11] = (4.0 * ans[10] - ans[9]) / 3.0;			 // good SM SM_10
	ans[12] = 0.0;										 // S_s

	//Additional states
	ans[13] = forcing_values[0] * c_1;
	ans[14] = q_pl; // q_pl;
	ans[15] = q_sl * A_h - q_b * 60.0;
	for (i = 0; i < num_parents; i++)
		ans[15] += y_p[i * dim + 15] * 60.0;
	ans[15] *= v_B / L;
	ans[16] = q_w_front; // wetting front(L) m
}

float flux_100(float theta_1, float theta_2, float K_sat, float psi_sat, float lambda, float theta_s, float theta_r, float S_L)
{
	//theta_1 = theta_1 <= theta_r ? theta_r + 0.01 : theta_1;
	//theta_2 = theta_2 <= theta_r ? theta_r + 0.01 : theta_2;
	//theta_1 = theta_1 > theta_s ? theta_s : theta_1;
	//theta_2 = theta_2 > theta_s ? theta_s : theta_2;
	float theta_hat = (theta_2 + theta_1) / 2.0;
	float phi_1 = (theta_1 - theta_r) / (theta_s - theta_r);
	float phi_2 = (theta_2 - theta_r) / (theta_s - theta_r);

	// K
	float phi_k = (theta_hat - theta_r) / (theta_s - theta_r);
	float power_k = (2.0 + 3.0 * lambda) / lambda;
	float K_phi_hat = K_sat * pow(phi_k, power_k);

	// psi
	float power_psi = -1.0 / lambda;
	float psi_1 = psi_sat * pow(phi_1, power_psi);
	float psi_2 = psi_sat * pow(phi_2, power_psi);

	return K_phi_hat * (1.0 - (psi_2 - psi_1) / S_L);
}

float flux_100_new(float theta_1, float theta_2, float K_sat, float psi_sat, float lambda, float theta_s, float theta_r, float S_L)
{
	float theta_hat = (theta_2 + theta_s) / 2.0;
	float phi_1 = (theta_s - theta_r) / (theta_s - theta_r);
	float phi_2 = (theta_2 - theta_r) / (theta_s - theta_r);

	// K
	float phi_k = (theta_hat - theta_r) / (theta_s - theta_r);
	float power_k = (2.0 + 3.0 * lambda) / lambda;
	float K_phi_hat = K_sat * pow(phi_k, power_k);

	// psi
	float power_psi = -1.0 / lambda;
	float psi_1 = psi_sat * pow(phi_1, power_psi);
	float psi_2 = psi_sat * pow(phi_2, power_psi);

	return K_phi_hat * (1.0 - (psi_2 - psi_1) / S_L);
}

//modeluid = 10000
void navid_100layer(double t, const double *const y_i, unsigned int dim, const double *const y_p, unsigned short num_parents, unsigned int max_dim, const double *const global_params, const double *const params, const double *const forcing_values, const QVSData *const qvs, int state, void *user, double *ans, double h)
{
	unsigned short i;

	// Model parameters (spatially uniform)
	float lambda_1 = global_params[1];
	float k_3 = global_params[4]; //[1/min]
	float h_b = global_params[6]; //[m]
								  //float S_L = global_params[7];	//[m]
	float A = global_params[8];
	float B = global_params[9];
	float exponent = global_params[10];
	float v_B = global_params[11];
	float psi_sat = global_params[12];

	// Model parameters (spatially variable)
	float L = params[1];   //[m]
	float A_h = params[2]; //[m^2]
	float K_sat = params[3];
	float invtau = params[4]; //[1/min]
	float k_2 = params[5];	  //[1/min]
	float k_i = params[6];	  //[1/min]
	float c_1 = params[7];
	float c_2 = params[8];
	float theta_s = (int)(params[9] * 1000) / 1000.0; // floorf(*1000) / 1000;
	float theta_r = params[10];
	float bc_lambda = params[11];
	float s_t[100];
	float q_t[99];
	float S_L = 0.05;
	float L_Top = 5;
	float dsp, ds0;
	// Initial conditions (or from last iteration)
	float q = y_i[0];	//[m^3/s]
	float s_p = y_i[1]; //[m]

	for (i = 0; i < 100; i++)
		s_t[i] = y_i[i + 2];
	float s_s = y_i[102];
	float q_b = y_i[104];											 //[m^3/s]
	float e_pot = forcing_values[1] * (1e-3 / (30.0 * 24.0 * 60.0)); //[mm/month] -> [m/min]
	float q_rain = forcing_values[0] * c_1;							 // m/min
	float q_sl = 0;													 // s_s > theta_r ? k_3 * (s_s - theta_r)  * (h_b - L_Top) : 0.0;	//[m/min]
	float q_pl, e_p, e_t;											 //

	q_pl = k_2 * s_p;
	for (i = 0; i < 99; i++)
		q_t[i] = flux_100(s_t[i], s_t[i + 1], K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);

	float Corr = s_p + (s_t[0] - theta_r) / (theta_s - theta_r); // +s_s / (h_b - L_Top
	if (e_pot > 0.0 && Corr > 1e-2)
	{
		e_p = s_p * e_pot / Corr;
		e_t = s_t[0] > theta_r + q_t[0] / S_L ? (s_t[0] - q_t[0] / S_L - theta_r) * e_pot : 0.0;
	}
	else
	{
		e_p = 0.0;
		e_t = 0.0;
	}
	float q_inf = q_t[0]; // flux_100(s_t[0], s_t[1], K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	float q_inf_new = flux_100_new(theta_s, s_t[1], K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	if (q_rain >= q_inf)
	{
		dsp = q_rain - q_inf - s_p * k_2;
		ds0 = (q_inf) / S_L;
	}
	else
	{
		if (s_p > 0)
		{
			dsp = -s_p * k_2;
			ds0 = (q_rain - q_inf) / S_L;
		}
		else
		{
			ds0 = (q_rain - q_inf) / S_L;
			dsp = 0.0;
		}
	}

	ans[0] = -q + (q_pl + q_sl) * c_2;
	for (i = 0; i < num_parents; i++)
		ans[0] += y_p[i * dim];
	ans[0] = invtau * pow(q, lambda_1) * ans[0];

	//Hillslope
	ans[1] = dsp; // S_P
	ans[2] = ds0; // S_t[0]
	for (i = 1; i < 99; i++)
		ans[i + 2] = (q_t[i - 1] - q_t[i]) / S_L; // S_t[1:98]
	ans[101] = (4.0 * ans[100] - ans[99]) / 3.0;  // S_t[99] good SM
	ans[102] = 0.0;								  // S_s
												  //Additional states
	ans[103] = forcing_values[0] * c_1;
	ans[104] = q_pl; // q_pl;
	ans[105] = q_sl * A_h - q_b * 60.0;
	for (i = 0; i < num_parents; i++)
		ans[105] += y_p[i * dim + 105] * 60.0;
	ans[105] *= v_B / L;
}

// model uid=10001
void navid_10layer(double t, const double *const y_i, unsigned int dim, const double *const y_p, unsigned short num_parents, unsigned int max_dim, const double *const global_params, const double *const params, const double *const forcing_values, const QVSData *const qvs, int state, void *user, double *ans, double h)
{
	unsigned short i;
	double s_t[10];
	double q_t[9];
	double dsp, ds0;

	// Parameters
	double lambda_1 = global_params[1];
	double k_3 = global_params[4]; //[1/min]
	double h_b = global_params[6]; //[m]
	//double S_L = global_params[7];	//[m]
	double A = global_params[8];
	double B = global_params[9];
	double exponent = global_params[10];
	double v_B = global_params[11];
	double e_pot = forcing_values[1] * (1e-3 / (30.0 * 24.0 * 60.0)); //[mm/month] -> [m/min]

	double L = params[1];	   //[m]
	double A_h = params[2];	   //[m^2]
							   //double h_r = params[3];	//[m]
	double invtau = params[3]; //[1/min]
	double k_2 = params[4];	   //[1/min]
	double k_i = params[5];	   //[1/min]
	double c_1 = params[6];
	double c_2 = params[7];

	double K_sat = params[8];
	//double psi_sat = global_params[12];

	double theta_s = (int)(params[9] * 1000) / 1000.0;
	double theta_r = params[10];
	double bc_lambda = params[11];
	double psi_sat = params[12];

	double S_L = 0.5;
	double L_Top = 5;

	// Initial conditions (or from last iteration)
	double q = y_i[0];	  //[m^3/s]
	double s_p = y_i[1];  //[m]
	double s_s = y_i[12]; //[m]
	double q_b = y_i[15]; //[m^3/s]

	for (i = 0; i < 10; i++)
		s_t[i] = y_i[i + 2];

	//
	double e_p, e_t;				  //
	double Corr = s_p + S_L * s_t[0]; // +s_s / (h_b - L_Top
	if (e_pot > 0.0 && Corr > 1e-2)
	{
		e_p = s_p * e_pot / Corr;
		e_t = s_t[0] > theta_r + q_t[0] / S_L ? (s_t[0] - theta_r) * e_pot / Corr : 0.0;
	}
	else
	{
		e_p = 0.0;
		e_t = 0.0;
	}

	double q_rain = forcing_values[0] * c_1; // m/min
	double q_sl = 0;						 // s_s > theta_r ? k_3 * (s_s - theta_r)  * (h_b - L_Top) : 0.0;	//[m/min]

	double q_pl = k_2 * s_p;
	for (i = 0; i < 9; i++)
		q_t[i] = flux_100(s_t[i], s_t[i + 1], K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);

	double q_inf = q_t[0]; // flux_100(s_t[0], s_t[1], K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	double q_inf_new = flux_100_new(theta_s, s_t[0], K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	if (q_rain >= q_inf)
	{
		if (q_rain >= q_inf_new)
		{
			dsp = q_rain - q_inf_new;
			ds0 = (q_inf_new - q_inf) / S_L;
		}
		else
		{
			dsp = 0;
			ds0 = (q_rain - q_inf) / S_L;
		}
	}
	else
	{
		if (s_p > 0)
		{
			dsp = q_rain - q_inf_new;
			ds0 = (q_inf_new - q_inf) / S_L;
		}
		else
		{
			ds0 = (q_rain - q_inf) / S_L;
			dsp = 0.0;
		}
	}

	//if (s_t[0] >= theta_s && q_rain >= q_t[0]) {
	//	dsp = q_rain - q_t[0] - s_p*k_2;
	//		s_t[0] = theta_s;
	//		ds0 = 0;
	//}
	//else {
	//	dsp = -s_p*k_2;
	//	ds0 = (q_rain - q_t[0]) / S_L;
	//}

	ans[0] = -q + (q_pl + q_sl) * c_2;
	for (i = 0; i < num_parents; i++)
		ans[0] += y_p[i * dim];
	ans[0] = invtau * pow(q, lambda_1) * ans[0];

	//Hillslope
	ans[1] = dsp - q_pl; // S_P
	ans[2] = ds0;		 // S_t[0]
	for (i = 1; i < 9; i++)
		ans[i + 2] = (q_t[i - 1] - q_t[i]) / S_L; // S_t[1:98]
	ans[11] = (4.0 * ans[10] - ans[9]) / 3.0;	  // S_t[99] good SM
	ans[12] = 0.0;								  // S_s
												  //Additional states
	ans[13] = forcing_values[0] * c_1;
	ans[14] = q_pl; // q_pl;
	ans[15] = q_sl * A_h - q_b * 60.0;
	for (i = 0; i < num_parents; i++)
		ans[15] += y_p[i * dim + 15] * 60.0;
	ans[15] *= v_B / L;
}

float flux_inf(float theta_1, float theta_2, float K_sat, float psi_sat, float lambda, float theta_s, float theta_r, float S_L1, float S_L2, int K_scheme)
{
	if (K_scheme == 0)
	{
		float depth_total = (S_L2 + S_L1);
		// Geometric average on theta
		float theta_hat = (S_L2 * theta_2 + S_L1 * theta_1) / depth_total;
		float phi_1 = (theta_1 - theta_r) / (theta_s - theta_r);
		float phi_2 = (theta_2 - theta_r) / (theta_s - theta_r);

		// K
		float phi_k = (theta_hat - theta_r) / (theta_s - theta_r);
		float power_k = (2.0 + 3.0 * lambda) / lambda;
		float K_phi_hat = K_sat * pow(phi_k, power_k);

		// psi
		float power_psi = -1.0 / lambda;
		float psi_1 = psi_sat * pow(phi_1, power_psi);
		float psi_2 = psi_sat * pow(phi_2, power_psi);

		return K_phi_hat * (1.0 - (psi_2 - psi_1) / S_L2);
	}
	else
	{
		float depth_total = (S_L2 + S_L1);
		//theta_1 = theta_1 > theta_s ? theta_s : theta_1;
		//theta_2 = theta_2 > theta_s ? theta_s : theta_2;
		float phi_1 = (theta_1 - theta_r) / (theta_s - theta_r);
		float phi_2 = (theta_2 - theta_r) / (theta_s - theta_r);
		// K
		float power_k = (2.0 + 3.0 * lambda) / lambda;
		float K_layer1 = K_sat * pow(phi_1, power_k);
		float K_layer2 = K_sat * pow(phi_2, power_k);
		// Geometric average for K
		float K_hat = (S_L1 * K_layer1 + S_L2 * K_layer2) / depth_total;

		// psi
		float power_psi = -1.0 / lambda;
		float psi_1 = psi_sat * pow(phi_1, power_psi);
		float psi_2 = psi_sat * pow(phi_2, power_psi);
		return K_hat * (1.0 - (psi_2 - psi_1) / S_L2);
	}
}

float flux_inf_pond(float theta_1, float theta_2, float K_sat, float psi_sat, float lambda, float theta_s, float theta_r, float S_L1, float S_L2, int K_scheme)
{

	if (K_scheme == 0)
	{
		float theta_hat = (theta_2 + theta_s) / 2.0;
		float phi_1 = (theta_s - theta_r) / (theta_s - theta_r);
		float phi_2 = (theta_2 - theta_r) / (theta_s - theta_r);

		// K
		float phi_k = (theta_hat - theta_r) / (theta_s - theta_r);
		float power_k = (2.0 + 3.0 * lambda) / lambda;
		float K_phi_hat = K_sat * pow(phi_k, power_k);

		// psi
		float power_psi = -1.0 / lambda;
		float psi_1 = psi_sat * pow(phi_1, power_psi);
		float psi_2 = psi_sat * pow(phi_2, power_psi);
		return K_phi_hat * (1.0 - (psi_2 - psi_1) / S_L2);
	}
	else
	{
		float phi_1 = (theta_s - theta_r) / (theta_s - theta_r);
		float phi_2 = (theta_2 - theta_r) / (theta_s - theta_r);
		// K
		float power_k = (2.0 + 3.0 * lambda) / lambda;
		float K_layer1 = K_sat * pow(phi_1, power_k);
		float K_layer2 = K_sat * pow(phi_2, power_k);
		// Simple averaging for K
		float K_hat = (K_layer1 + K_layer2) / 2.0;
		//TODO:
		//Add geometric Averaging for K

		// psi
		float power_psi = -1.0 / lambda;
		float psi_1 = psi_sat * pow(phi_1, power_psi);
		float psi_2 = psi_sat * pow(phi_2, power_psi);
		return K_hat * (1.0 - (psi_2 - psi_1) / S_L2);
	}
}

//MODEL UID 1011

float flux_gampt(float theta_1, float theta_2, float K_sat, float psi_sat, float lambda, float theta_s, float theta_r, float S_L, double s_p)
{
	float dz = (s_p + S_L) / 2.0;
	float phi_1 = (theta_1 - theta_r) / (theta_s - theta_r);
	float phi_2 = (theta_2 - theta_r) / (theta_s - theta_r);
	// K
	float power_k = (2.0 + 3.0 * lambda) / lambda;
	float K_layer1 = K_sat * pow(phi_1, power_k);
	float K_layer2 = K_sat * pow(phi_2, power_k);
	// Simple averaging for K
	float K_hat = (K_sat * s_p + K_layer2 * S_L) / (s_p + S_L);
	//TODO:
	//Add geometric Averaging for K

	// psi
	float power_psi = -1.0 / lambda;
	float psi_1 = psi_sat * pow(phi_1, power_psi);
	float psi_2 = psi_sat * pow(phi_2, power_psi);
	return K_hat * (-psi_2 - s_p) / dz;
}
void herrada_gampt(double t, const double *const y_i, unsigned int dim, const double *const y_p, unsigned short num_parents, unsigned int max_dim, const double *const global_params, const double *const params, const double *const forcing_values, const QVSData *const qvs, int state, void *user, double *ans, double h)
{
	unsigned short i;
	double S_L[10] = {0.01, 0.05, 0.05, 0.10, 0.30, 0.50, 1.0, 1.0, 1.0, 1.0}; // layer depths [m]
																			   // Model parameters (spatially uniform)
	double s_t[10];
	double q_t[9];
	double dsp, ds0;
	double e_p, e_t; //
	int K_scheme = 1;
	float lambda_1 = global_params[1];
	float k_3 = global_params[4]; //[1/min]
	float h_b = global_params[6]; //[m]
	//float S_L = global_params[7];	//[m]
	float A = global_params[8];
	float B = global_params[9];
	float exponent = global_params[10];
	float v_B = global_params[11];
	float psi_sat = global_params[12];

	// Model parameters (spatially variable)
	float L = params[1];   //[m]
	float A_h = params[2]; //[m^2]
	float K_sat = params[3];
	float invtau = params[4]; //[1/min]
	float k_2 = params[5];	  //[1/min]
	float k_i = params[6];	  //[1/min]
	float c_1 = params[7];
	float c_2 = params[8];
	float theta_s = (int)(params[9] * 1000) / 1000.0; // floorf(*1000) / 1000;
	float theta_r = params[10];
	float bc_lambda = params[11];

	//float L_Top = S_L * 10.0;

	// Initial conditions (or from last iteration)
	float q = y_i[0];	//[m^3/s]
	float s_p = y_i[1]; //[m]
	for (i = 0; i < 10; i++)
	{
		s_t[i] = y_i[i + 2];
		//}
		if (s_t[i] > theta_s)
		{
			s_t[i] = theta_s - 0.01;
			//printf("%s", "hello");
		}
		else if (s_t[i] < theta_r)
		{
			s_t[i] = theta_r + 0.01;
		}
	}
	float s_s = y_i[12];
	float q_b = y_i[15];
	float w_front = y_i[16];
	//[m^3/s]

	float q_rain = forcing_values[0] * c_1; // m/min
											//if (q_rain == 0.0) {
	w_front = w_front > 0.2 ? 0.01 : w_front;
	//}
	float e_pot = forcing_values[1] * (1e-3 / (30.0 * 24.0 * 60.0)); //[mm/month] -> [m/min]
	for (i = 0; i < 9; i++)
	{
		q_t[i] = flux_inf(s_t[i], s_t[i + 1], K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L[i], S_L[i + 1], K_scheme);
	}

	//float q_ts = flux_ts(s_10, s_s, K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L);
	//float q_sl = s_s> theta_r ? k_3 * (s_s - theta_r)  * (h_b - L_Top) : 0.0;	//[m/min]
	float q_pl; // , ds_1, e_p, e_t;
	//float S = (s_1 + s_2) / (2 * theta_s) > 1 ? 0.99 : (s_1 + s_2) / (2 * theta_s);
	//float q_infilt = K_sat * (s_p + w_front) / w_front;
	//float q_infilt = K_sat * s_p + s_1 * (s_p + w_front) / w_front;
	float q_infilt = flux_gampt(theta_s, s_t[0], K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L[0], s_p);
	float q_thin = flux_gampt(theta_s, s_t[0], K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L[0], 0.001);
	//float q_w_front = K_sat * (s_p + w_front) / (w_front*theta_s*(1 - S));
	q_pl = s_p > 0.001 ? k_2 * s_p : 0;
	// Evapotranspiration correction for top layer and ponding states
	//float Corr = s_p + (s_1 - theta_r) / (theta_s - theta_r);// +s_s / (h_b - L_Top);
	//if (e_pot > 0.0 && Corr > 1e-2)
	//{
	//	e_p = s_p * e_pot / Corr;
	//	e_t = s_1 > theta_r + q_2 / S_L ? (s_1 - q_2 / S_L - theta_r)   * e_pot : 0.0;
	//}
	//else
	//{
	//	e_p = 0.0;
	//	e_t = 0.0;
	//}

	double Corr = (s_t[0] - theta_r) / (theta_s - theta_r); // +s_s / (h_b - L_Top
	if (e_pot > 0.0 && Corr > 1e-2)
	{
		e_p = s_p * e_pot / Corr;
		e_t = s_t[0] > theta_r + q_t[0] / S_L[0] ? (s_t[0] - theta_r) / (theta_s - theta_r) * e_pot : 0.0;
	}
	else
	{
		e_p = 0.0;
		e_t = 0.0;
	}
	double extra_water;
	double q_inf = s_p > 0.0 ? q_infilt : q_t[0];
	if (q_rain > q_inf || s_t[0] > theta_s)
	{
		if (s_p == 0.0)
		{

			dsp = q_rain - q_inf - q_pl;
			ds0 = (-e_t) / S_L[0];
			extra_water = ds0 > theta_s - s_t[0] ? ds0 - (theta_s - s_t[0]) : 0.0;
			dsp = dsp + extra_water * S_L[0];
			ds0 = ds0 - extra_water;
		}
		else
		{
			dsp = q_rain - q_inf - q_pl;
			ds0 = 0.0;
		}
	}
	else
	{
		dsp = -q_pl;
		ds0 = (q_rain - q_t[0] - e_t) / S_L[0];
		extra_water = ds0 > theta_s - s_t[0] ? ds0 - (theta_s - s_t[0]) : 0.0;
		dsp = dsp + extra_water * S_L[0];
		ds0 = ds0 - extra_water;
	}

	double q_sl = k_3 * s_t[5] * S_L[5];
	ans[0] = -q + (q_pl + q_sl) * c_2;
	for (i = 0; i < num_parents; i++)
		ans[0] += y_p[i * dim];
	ans[0] = invtau * pow(q, lambda_1) * ans[0];

	//Hillslope
	ans[1] = dsp; // S_P
	ans[2] = ds0; // S_1
	for (i = 1; i < 9; i++)
	{
		ans[i + 2] = (q_t[i - 1] - q_t[i]) / S_L[i]; // S_t[1:98]
	}
	ans[11] = (4.0 * ans[10] - ans[9]) / 3.0; // good SM SM_10
	ans[12] = 0.0;							  // S_s

	//Additional states
	ans[13] = forcing_values[0] * c_1;
	ans[14] = q_pl; // q_pl;
	ans[15] = q_sl * A_h - q_b * 60.0;
	for (i = 0; i < num_parents; i++)
		ans[15] += y_p[i * dim + 15] * 60.0;
	ans[15] *= v_B / L;
	ans[16] = 0; // wetting front(L) m
}

// model_uid = 10002
void navid_10layer_new(double t, const double *const y_i, unsigned int dim, const double *const y_p, unsigned short num_parents, unsigned int max_dim, const double *const global_params, const double *const params, const double *const forcing_values, const QVSData *const qvs, int state, void *user, double *ans, double h)
{

	unsigned short i;
	double s_t[10];
	double q_t[9];
	double dsp, ds0;
	double e_p, e_t; //

	double lambda_1 = global_params[1];
	double k_3 = global_params[4]; //[1/min]
	double h_b = global_params[6]; //[m]

	double A = global_params[8];
	double B = global_params[9];
	double exponent = global_params[10];
	double v_B = global_params[11];
	int K_scheme = 1;

	double e_pot = forcing_values[1] * (1e-3 / (30.0 * 24.0 * 60.0)); //[mm/month] -> [m/min]

	double L = params[1];	   //[m]
	double A_h = params[2];	   //[m^2]
							   //double h_r = params[3];	//[m]
	double invtau = params[3]; //[1/min]
	double k_2 = params[4];	   //[1/min]
	double k_i = params[5];	   //[1/min]
	double c_1 = params[6];
	double c_2 = params[7];

	double K_sat = params[8];
	//double psi_sat = global_params[12];
	//printf("%i", link->h);
	//printf("%s", link_i->ID);
	double theta_s = (int)(params[9] * 1000) / 1000.0;
	double theta_r = params[10];
	double bc_lambda = params[11];
	double psi_sat = params[12];

	double S_L[10] = {0.01, 0.05, 0.05, 0.10, 0.30, 0.50, 1.0, 1.0, 1.0, 1.0}; // layer depths [m]
	double L_Top = 5.0;

	// Initial conditions (or from last iteration)
	double q = y_i[0];	  //[m^3/s]
	double s_p = y_i[1];  //[m]
	double s_s = y_i[12]; //[m]
	double q_b = y_i[15]; //[m^3/s]
	double extra_water = 0.0;
	////////////
	for (i = 0; i < 10; i++)
		s_t[i] = y_i[i + 2];

	//printf("%f\n", t);
	double q_rain = forcing_values[0] * c_1; // m/min
	double q_sl = k_3 * s_t[5] * S_L[5];	 // s_s > theta_r ? k_3 * (s_s - theta_r)  * (h_b - L_Top) : 0.0;	//[m/min]

	double q_pl = k_2 * s_p;

	for (i = 0; i < 10; i++)
	{
		s_t[i] = y_i[i + 2];
		//}
		if (s_t[i] > theta_s)
		{
			s_t[i] = theta_s;
			//printf("%s", "hello");
		}
		else if (s_t[i] < theta_r)
		{
			s_t[i] = theta_r + 0.01;
		}
	}
	for (i = 0; i < 9; i++)
		q_t[i] = flux_inf(s_t[i], s_t[i + 1], K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L[i], S_L[i + 1], K_scheme);

	double Corr = Corr = (s_t[0] - theta_r) / (theta_s - theta_r); // +s_s / (h_b - L_Top
	if (e_pot > 0.0 && Corr > 1e-2)
	{
		e_p = s_p * e_pot / Corr;
		e_t = s_t[0] > theta_r + q_t[0] / S_L[0] ? (s_t[0] - theta_r) / (theta_s - theta_r) * e_pot : 0.0;
	}
	else
	{
		e_p = 0.0;
		e_t = 0.0;
	}

	//double extra_water;
	double q_inf = q_t[0];
	double q_pond_inf = flux_inf_pond(theta_s, s_t[0], K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L[0], S_L[1], K_scheme);
	if (q_rain > q_pond_inf)
	{

		dsp = q_rain - q_pond_inf - q_pl;
		ds0 = (q_pond_inf - q_inf) / S_L[0];
		extra_water = ds0 + s_t[0] > theta_s ? ds0 - (theta_s - s_t[0]) : 0.0;
		dsp = q_rain - q_pond_inf - q_pl + extra_water * S_L[0];
		ds0 = (q_pond_inf - q_inf) / S_L[0] - extra_water;
	}
	else
	{
		if (q_rain > q_inf && s_t[0] >= theta_s)
		{
			s_t[0] = theta_s;
			q_inf = flux_inf(s_t[0], s_t[1], K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L[0], S_L[1], K_scheme);
			dsp = q_rain - q_inf - q_pl;
			ds0 = 0.0;
		}
		else
		{

			dsp = -q_pl;
			ds0 = (q_rain - q_inf) / S_L[0];
			// extra_water = ds0 > theta_s - s_t[0] ? ds0 - (theta_s - s_t[0]) : 0.0;
			// dsp = - q_pl + extra_water * S_L[0];
			// ds0 = (q_rain - q_inf - e_t) / S_L[0] - extra_water;
		}
	}

	ans[0] = -q + (q_pl + q_sl) * c_2;
	for (i = 0; i < num_parents; i++)
		ans[0] += y_p[i * dim];
	ans[0] = invtau * pow(q, lambda_1) * ans[0];

	//Hillslope
	ans[1] = dsp; // S_P
	ans[2] = ds0; // S_t[0]
	for (i = 1; i < 9; i++)
	{
		ans[i + 2] = (q_t[i - 1] - q_t[i]) / S_L[i]; // S_t[1:98]
	}

	// ans[7] = (q_t[4] - q_t[5]) / S_L[5] - k_3 * s_t[5];

	ans[11] = (4.0 * ans[10] - ans[9]) / 3.0; // S_t[99] good SM
	ans[12] = 0.0;							  // S_s
											  //Additional states
	ans[13] = forcing_values[0] * c_1;
	ans[14] = q_pl; // q_pl;
	ans[15] = q_sl * A_h - q_b * 60.0;
	for (i = 0; i < num_parents; i++)
		ans[15] += y_p[i * dim + 15] * 60.0;
	ans[15] *= v_B / L;
}

void changeArray(double *y_i, double s_p)
{
	y_i[1] = s_p;
}

// model_uid = 10005
void navid_extra_water(double t, const double *const y_i, unsigned int dim, const double *const y_p, unsigned short num_parents, unsigned int max_dim, const double *const global_params, const double *const params, const double *const forcing_values, const QVSData *const qvs, int state, void *user, double *ans, double h)
{
	unsigned short i;
	double s_t[10];
	double q_t[9];
	double dsp, ds0;
	double e_p, e_t[5]; //

	double lambda_1 = global_params[1];
	double k_3 = global_params[4]; //[1/min]
	double h_b = global_params[6]; //[m]

	double A = global_params[8];
	double B = global_params[9];
	double exponent = global_params[10];
	double v_B = global_params[11];
	int K_scheme = 1;

	double e_pot = forcing_values[1] * (1e-3 / (30.0 * 24.0 * 60.0)); //[mm/month] -> [m/min]

	double L = params[1];	   //[m]
	double A_h = params[2];	   //[m^2]
							   //double h_r = params[3];	//[m]
	double invtau = params[3]; //[1/min]
	double k_2 = params[4];	   //[1/min]
	double k_i = params[5];	   //[1/min]
	double c_1 = params[6];
	double c_2 = params[7];

	double K_sat = params[8];
	//double psi_sat = global_params[12];
	//printf("%i", link->h);
	//printf("%s", link_i->ID);
	double theta_s = (int)(params[9] * 1000) / 1000.0;
	double theta_r = params[10];
	double bc_lambda = params[11];
	double psi_sat = params[12];

	double S_L[10] = {0.01, 0.05, 0.05, 0.10, 0.30, 0.50, 1.0, 1.0, 1.0, 1.0}; // layer depths [m]
	double L_Top = 5.0;

	// Initial conditions (or from last iteration)
	double q = y_i[0];	  //[m^3/s]
	double s_p = y_i[1];  //[m]
	double s_s = y_i[12]; //[m]
	double q_b = y_i[15]; //[m^3/s]
	double extra_flux = 0.0;
	double ds_p = 0.0;
	double d_ss = 0.0;
	for (i = 0; i < 10; i++)
	{
		s_t[i] = y_i[i + 2];
		if (s_t[i] > theta_s)
		{
			d_ss = (s_t[i] - theta_s) * S_L[i] / h;
			// ds_p = d_ss;
			s_t[i] = theta_s;
		}
		else if (s_t[i] < theta_r)
		{
			s_t[i] = theta_r + 0.01;
		}
	}
	double q_rain = forcing_values[0] * c_1; // m/min

	double q_sl = k_3 * s_t[5] * S_L[5]; // s_s > theta_r ? k_3 * (s_s - theta_r)  * (h_b - L_Top) : 0.0;	//[m/min]
	double q_pl = k_2 * s_p;
	for (i = 0; i < 9; i++)
		q_t[i] = flux_inf(s_t[i], s_t[i + 1], K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L[i], S_L[i + 1], K_scheme);
	double et_tot=0.0;
	double Corr = (s_t[0]-theta_r)*S_L[0] + (s_t[1]-theta_r)*S_L[1] + (s_t[2]-theta_r)*S_L[2] + (s_t[3]-theta_r)*S_L[3] + (s_t[4]-theta_r)*S_L[4];
	// if (e_pot > 0.0 && Corr > 1e-2)
	// {
	for (i = 0; i < 5; i++)
	{
		e_t[i] = s_t[i]> theta_r + 0.02? (s_t[i]-theta_r)*S_L[i]*e_pot/Corr:0;
		et_tot += e_t[i];
	}
	// }

	double q_inf = q_t[0];
	double q_pond_inf = flux_inf_pond(theta_s, s_t[0], K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L[0], S_L[1], K_scheme);
	if (q_rain > q_pond_inf)
	{
		dsp = q_rain - q_pond_inf - q_pl;
		ds0 = (q_pond_inf - q_inf - e_t[0]) / S_L[0];
		extra_flux = ds0 * h + s_t[0] > theta_s ? ds0 - (theta_s - s_t[0]) / h : 0.0;
		dsp = q_rain - q_pond_inf - q_pl + extra_flux * S_L[0];
		ds0 = (q_pond_inf - q_inf - e_t[0]) / S_L[0] - extra_flux;
	}
	else
	{
		dsp = -q_pl;
		ds0 = (q_rain - q_inf - e_t[0]) / S_L[0];
		extra_flux = ds0 * h + s_t[0] > theta_s ? ds0 - (theta_s - s_t[0]) / h : 0.0;
		dsp = -q_pl + extra_flux * S_L[0];
		ds0 = (q_rain - q_inf - e_t[0]) / S_L[0] - extra_flux;
	}

	ans[0] = -q + (q_pl + q_sl) * c_2;
	for (i = 0; i < num_parents; i++)
		ans[0] += y_p[i * dim];
	ans[0] = invtau * pow(q, lambda_1) * ans[0];

	//Hillslope
	ans[1] = dsp; // S_P
	ans[2] = ds0; // S_t[0]
	for (i = 1; i < 5; i++)
	{
		ans[i + 2] = (q_t[i - 1] - q_t[i]-e_t[i]) / S_L[i]; // S_t[1:98]
	}
	for (i = 5; i < 9; i++)
	{
		ans[i + 2] = (q_t[i - 1] - q_t[i]) / S_L[i]; // S_t[1:98]
	}
	ans[7] = (q_t[4] - q_t[5]) / S_L[5] - k_3 * s_t[5];

	ans[11] = (4.0 * ans[10] - ans[9]) / 3.0; // S_t[99] good SM
	ans[12] = et_tot;							  // S_s
											  //Additional states
	ans[13] = forcing_values[0] * c_1;
	ans[14] = q_pl; // q_pl;
	ans[15] = q_sl * A_h - q_b * 60.0;
	for (i = 0; i < num_parents; i++)
		ans[15] += y_p[i * dim + 15] * 60.0;
	ans[15] *= v_B / L;
}

// model_uid = 10006
void navid_w_et(double t, const double *const y_i, unsigned int dim, const double *const y_p, unsigned short num_parents, unsigned int max_dim, const double *const global_params, const double *const params, const double *const forcing_values, const QVSData *const qvs, int state, void *user, double *ans, double h)
{
	unsigned short i;
	double s_t[10];
	double q_t[9];
	double dsp, ds0;
	double e_p, e_t[5]; //

	double lambda_1 = global_params[1];
	double k_3 = global_params[4]; //[1/min]
	double h_b = global_params[6]; //[m]

	double A = global_params[8];
	double B = global_params[9];
	double exponent = global_params[10];
	double v_B = global_params[11];
	int K_scheme = 1;

	double e_pot = forcing_values[1] * (1e-3 / (30.0 * 24.0 * 60.0)); //[mm/month] -> [m/min]

	double L = params[1];	   //[m]
	double A_h = params[2];	   //[m^2]
							   //double h_r = params[3];	//[m]
	double invtau = params[3]; //[1/min]
	double k_2 = params[4];	   //[1/min]
	double k_i = params[5];	   //[1/min]
	double c_1 = params[6];
	double c_2 = params[7];

	double K_sat = params[8];
	//double psi_sat = global_params[12];
	//printf("%i", link->h);
	//printf("%s", link_i->ID);
	double theta_s = (int)(params[9] * 1000) / 1000.0;
	double theta_r = params[10];
	double bc_lambda = params[11];
	double psi_sat = params[12];

	double S_L[10] = {0.01, 0.05, 0.05, 0.10, 0.30, 0.50, 1.0, 1.0, 1.0, 1.0}; // layer depths [m]
	double L_Top = 5.0;

	// Initial conditions (or from last iteration)
	double q = y_i[0];	  //[m^3/s]
	double s_p = y_i[1];  //[m]
	double s_s = y_i[12]; //[m]
	double q_b = y_i[15]; //[m^3/s]
	double extra_flux = 0.0;
	double ds_p = 0.0;
	double d_ss = 0.0;
	for (i = 0; i < 10; i++)
	{
		s_t[i] = y_i[i + 2];
		if (s_t[i] > theta_s)
		{
			d_ss = (s_t[i] - theta_s) * S_L[i] / h;
			// ds_p = d_ss;
			s_t[i] = theta_s;
		}
		else if (s_t[i] < theta_r)
		{
			s_t[i] = theta_r + 0.01;
		}
	}
	double q_rain = forcing_values[0] * c_1; // m/min

	double q_sl = k_3 * s_t[5] * S_L[5]; // s_s > theta_r ? k_3 * (s_s - theta_r)  * (h_b - L_Top) : 0.0;	//[m/min]
	double q_pl = k_2 * s_p;
	for (i = 0; i < 9; i++)
		q_t[i] = flux_inf(s_t[i], s_t[i + 1], K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L[i], S_L[i + 1], K_scheme);

	double Corr = (s_t[0]-theta_r)*S_L[0] + (s_t[1]-theta_r)*S_L[1] + (s_t[2]-theta_r)*S_L[2] + (s_t[3]-theta_r)*S_L[3] + (s_t[4]-theta_r)*S_L[4];
	// if (e_pot > 0.0 && Corr > 1e-2)
	// {
	
	double et_tot=0.0;
	double e_pot_rem = e_pot;
	for (i = 0; i < 5; i++)
	{
		double s_lim = s_t[i] - theta_r - 0.05;
		double C_et = s_lim/pow(0.0005 + pow(s_lim,2.0),0.5);
		e_t[i] = s_t[i]-C_et*e_pot_rem*h > theta_r + q_t[i]*h / S_L[i] ? C_et*e_pot_rem : 0.0 ;
		e_pot_rem = e_pot_rem-e_t[i];
		et_tot += e_t[i];
	}
	// }

	double q_inf = q_t[0];
	double q_pond_inf = flux_inf_pond(theta_s, s_t[0], K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L[0], S_L[1], K_scheme);
	if (q_rain > q_pond_inf)
	{
		dsp = q_rain - q_pond_inf - q_pl;
		ds0 = (q_pond_inf - q_inf - e_t[0]) / S_L[0];
		extra_flux = ds0 * h + s_t[0] > theta_s ? ds0 - (theta_s - s_t[0]) / h : 0.0;
		dsp = q_rain - q_pond_inf - q_pl + extra_flux * S_L[0];
		ds0 = (q_pond_inf - q_inf - e_t[0]) / S_L[0] - extra_flux;
	}
	else
	{
		dsp = -q_pl;
		ds0 = (q_rain - q_inf - e_t[0]) / S_L[0];
		extra_flux = ds0 * h + s_t[0] > theta_s ? ds0 - (theta_s - s_t[0]) / h : 0.0;
		dsp = -q_pl + extra_flux * S_L[0];
		ds0 = (q_rain - q_inf - e_t[0]) / S_L[0] - extra_flux;
	}

	ans[0] = -q + (q_pl + q_sl) * c_2;
	for (i = 0; i < num_parents; i++)
		ans[0] += y_p[i * dim];
	ans[0] = invtau * pow(q, lambda_1) * ans[0];

	//Hillslope
	ans[1] = dsp; // S_P
	ans[2] = ds0; // S_t[0]
	for (i = 1; i < 5; i++)
	{
		ans[i + 2] = (q_t[i - 1] - q_t[i]-e_t[i]) / S_L[i]; // S_t[1:98]
	}
	for (i = 5; i < 9; i++)
	{
		ans[i + 2] = (q_t[i - 1] - q_t[i]) / S_L[i]; // S_t[1:98]
	}
	ans[7] = (q_t[4] - q_t[5]) / S_L[5] - k_3 * s_t[5];

	ans[11] = (4.0 * ans[10] - ans[9]) / 3.0; // S_t[99] good SM
	ans[12] = et_tot;							  // S_s
											  //Additional states
	ans[13] = forcing_values[0] * c_1;
	ans[14] = q_pl; // q_pl;
	ans[15] = q_sl * A_h - q_b * 60.0;
	for (i = 0; i < num_parents; i++)
		ans[15] += y_p[i * dim + 15] * 60.0;
	ans[15] *= v_B / L;
}


// model_uid = 10007
void navid_w_et_kt(double t, const double *const y_i, unsigned int dim, const double *const y_p, unsigned short num_parents, unsigned int max_dim, const double *const global_params, const double *const params, const double *const forcing_values, const QVSData *const qvs, int state, void *user, double *ans, double h)
{
	unsigned short i;
	double s_t[10];
	double q_t[9];
	double dsp, ds0;
	double e_p, e_t[5]; //

	double lambda_1 = global_params[1];
	double k_3 = global_params[4]; //[1/min]
	double h_b = global_params[6]; //[m]

	double A = global_params[8];
	double B = global_params[9];
	double exponent = global_params[10];
	double v_B = global_params[11];
	int K_scheme = 1;

	double e_pot = forcing_values[1] * (1e-3 / (30.0 * 24.0 * 60.0)); //[mm/month] -> [m/min]

	double L = params[1];	   //[m]
	double A_h = params[2];	   //[m^2]
							   //double h_r = params[3];	//[m]
	double invtau = params[3]; //[1/min]
	double k_2 = params[4];	   //[1/min]
	double k_i = params[5];	   //[1/min]
	double c_1 = params[6];
	double c_2 = params[7];

	double K_sat = params[8];
	//double psi_sat = global_params[12];
	//printf("%i", link->h);
	//printf("%s", link_i->ID);
	double theta_s = (int)(params[9] * 1000) / 1000.0;
	double theta_r = params[10];
	double bc_lambda = params[11];
	double psi_sat = params[12];

	double S_L[10] = { 0.05, 0.05, 0.05, 0.10, 0.30, 0.50, 1.0, 1.0, 1.0, 1.0 }; // layer depths [m]
	double L_Top = 5.0;

	// Initial conditions (or from last iteration)
	double q = y_i[0];	  //[m^3/s]
	double s_p = y_i[1];  //[m]
	double s_s = y_i[12]; //[m]
	double q_b = y_i[15]; //[m^3/s]
	double extra_flux = 0.0;
	double ds_p = 0.0;
	double d_ss = 0.0;
	for (i = 0; i < 10; i++)
	{
		s_t[i] = y_i[i + 2];
		if (s_t[i] > theta_s)
		{
			d_ss = (s_t[i] - theta_s) * S_L[i] / h;
			// ds_p = d_ss;
			s_t[i] = theta_s;
		}
		else if (s_t[i] < theta_r)
		{
			s_t[i] = theta_r + 0.01;
		}
	}
	double q_rain = forcing_values[0] * c_1; // m/min

	double q_sl = k_3 * s_t[5] * S_L[5]; // s_s > theta_r ? k_3 * (s_s - theta_r)  * (h_b - L_Top) : 0.0;	//[m/min]
	for (i = 0; i < 9; i++)
		q_t[i] = flux_inf(s_t[i], s_t[i + 1], K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L[i], S_L[i + 1], K_scheme);

	double Corr = (s_t[0] - theta_r)*S_L[0] + (s_t[1] - theta_r)*S_L[1] + (s_t[2] - theta_r)*S_L[2] + (s_t[3] - theta_r)*S_L[3] + (s_t[4] - theta_r)*S_L[4];


	double et_tot = 0.0;
	double e_pot_rem = e_pot;
	// for (i = 0; i < 5; i++)
	// {
	// 	double s_lim = s_t[i] - 1.4*theta_r;
	// 	double C_et = (0.5* s_lim / pow(0.000005 + pow(s_lim, 2.0), 0.5)) + 0.5;
	// 	e_t[i] = s_t[i] - q_t[i] * h> theta_r + 0.02 ? C_et*e_pot_rem : 0;
	// 	e_pot_rem = e_pot_rem - e_t[i];
	// 	et_tot += e_t[i];
	// }
	for (i = 0; i < 5; i++)
	{
		double s_lim = s_t[i] - theta_r - 0.05;
		double C_et = s_lim/pow(0.0005 + pow(s_lim,2.0),0.5);
		e_t[i] = C_et > theta_r + q_t[i] / S_L[i] ? C_et*e_pot_rem : 0.0 ;
		e_pot_rem = e_pot_rem-e_t[i];
		et_tot += e_t[i];
	}

	double q_inf = q_t[0];
	double q_pond_inf = flux_inf_pond(theta_s, s_t[0], K_sat, psi_sat, bc_lambda, theta_s, theta_r, S_L[0], S_L[1], K_scheme);
	/*if (q_rain > q_inf)
	{
		dsp = q_rain - q_pond_inf - q_pl;
		ds0 = (q_pond_inf - q_inf - e_t[0]) / S_L[0];
		extra_flux = ds0 * h + s_t[0] > theta_s ? ds0 - (theta_s - s_t[0]) / h : 0.0;
		dsp = q_rain - q_pond_inf - q_pl + extra_flux * S_L[0];
		ds0 = (q_pond_inf - q_inf - e_t[0]) / S_L[0] - extra_flux;
	}
	else
	{
		dsp = -q_pl;
		ds0 = (q_rain - q_inf - e_t[0]) / S_L[0];
		extra_flux = ds0 * h + s_t[0] > theta_s ? ds0 - (theta_s - s_t[0]) / h : 0.0;
		dsp = -q_pl + extra_flux * S_L[0];
		ds0 = (q_rain - q_inf - e_t[0]) / S_L[0] - extra_flux;
	}*/
	double pow_term = (1.0 - (s_t[0]-theta_r)/(theta_s-theta_r) > 0.0) ? pow(1.0 - (s_t[0] - theta_r) / (theta_s - theta_r), exponent) : 0.0;
	double k_t = (A + B * pow_term) * k_2;

	//Fluxes
	double q_pl = k_2 * s_p;
	double q_pt = k_t * s_p;

	ans[0] = -q + (q_pl + q_sl) * c_2;
	for (i = 0; i < num_parents; i++)
		ans[0] += y_p[i * dim];
	ans[0] = invtau * pow(q, lambda_1) * ans[0];

	//Hillslope
	ans[1] = forcing_values[0] * c_1 - q_pl - q_pt; // S_P
	ans[2] = (q_pt - q_inf -e_t[0])/S_L[0]; // S_t[0]
	for (i = 1; i < 5; i++)
	{
		ans[i + 2] = (q_t[i - 1] - q_t[i] - e_t[i]) / S_L[i]; // S_t[1:98]
	}
	for (i = 5; i < 9; i++)
	{
		ans[i + 2] = (q_t[i - 1] - q_t[i]) / S_L[i]; // S_t[1:98]
	}
	ans[7] = (q_t[4] - q_t[5]) / S_L[5] - k_3 * s_t[5];

	ans[11] = (4.0 * ans[10] - ans[9]) / 3.0; // S_t[99] good SM
	ans[12] = et_tot;							  // S_s
												  //Additional states
	ans[13] = forcing_values[0] * c_1;
	ans[14] = q_pl; // q_pl;
	ans[15] = q_sl * A_h - q_b * 60.0;
	for (i = 0; i < num_parents; i++)
		ans[15] += y_p[i * dim + 15] * 60.0;
	ans[15] *= v_B / L;
}
