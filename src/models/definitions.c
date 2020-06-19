#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if defined(HAVE_UNISTD_H)
#include <unistd.h>
#endif

#include <rksteppers.h>
#include <models/definitions.h>
#include <models/equations.h>
#include <models/check_consistency.h>
#include <models/output_constraints.h>
#include <models/check_state.h>


//Sets the various sizes and flags for the model. This method should set the following fields:
//dim:			The number of unknowns in the differential equations (or the number of ODEs at each link).
//diff_start:		The starting index of the differential unknowns.
//str_flag:		1 if reading rainfall data from a .str file, 0 else.
//binrain_flag:		1 if reading rainfall data from binary files, 0 else.
//uses_dam:		1 if dams are compatible with the model given by model_uid.
//params_size:		The total number of parameters (including precalculations) at links with no dam.
//dam_params_size:	The total number of parameters (including precalculations) at links with a dam.
//area_idx:		The entry in params where the upstream area is stored.
//disk_params:		The number of enries in param that are read from DEM data.
//Currently, this program assumes the same number of differential equations at each link.
//UnivVars* GlobalVars:	Contains the global variables for the system.
void SetParamSizes(
	GlobalVars* globals,
	void* external)
{
	unsigned short int model_uid = globals->model_uid;
	unsigned int num_global_params;

	//Set dim and start of differential variables
	switch (model_uid)
	{
	case 1:	num_global_params = 12;
		globals->uses_dam = 0;
		globals->num_params = 8;
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
		globals->num_disk_params = 3;
		globals->convertarea_flag = 0;
		globals->num_forcings = 3;
		globals->min_error_tolerances = 14;
		break;
	case 8:	num_global_params = 12;
		globals->uses_dam = 0;
		globals->num_params = 8;
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
		globals->num_disk_params = 3;
		globals->convertarea_flag = 0;
		globals->num_forcings = 3;
		globals->min_error_tolerances = 14;
		break;
	case 252:	num_global_params = 11;
		globals->uses_dam = 0;
		globals->num_params = 8;
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
		globals->num_disk_params = 3;
		globals->convertarea_flag = 0;
		globals->num_forcings = 2;
		globals->min_error_tolerances = 4;
		break;
	case 254:	num_global_params = 12;
		globals->uses_dam = 0;
		globals->num_params = 8;
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
		globals->num_disk_params = 3;
		globals->convertarea_flag = 0;
		globals->num_forcings = 3;
		globals->min_error_tolerances = 7;
		break;
	case 1004:	num_global_params = 12;
		globals->uses_dam = 0;
		globals->num_params = 12;
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
		globals->num_disk_params = 7;
		globals->convertarea_flag = 0;
		globals->num_forcings = 3;
		globals->min_error_tolerances = 10;
		break;
	case 1005:	num_global_params = 13;
		globals->uses_dam = 0;
		globals->num_params = 12;
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
		globals->num_disk_params = 7;
		globals->convertarea_flag = 0;
		globals->num_forcings = 3;
		globals->min_error_tolerances = 13;
		break;
	case 1006:	num_global_params = 13;
		globals->uses_dam = 0;
		globals->num_params = 12;
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
		globals->num_disk_params = 7;
		globals->convertarea_flag = 0;
		globals->num_forcings = 3;
		globals->min_error_tolerances = 13;
		break;
	case 1007:	num_global_params = 13;
		globals->uses_dam = 0;
		globals->num_params = 12;
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
		globals->num_disk_params = 7;
		globals->convertarea_flag = 0;
		globals->num_forcings = 3;
		globals->min_error_tolerances = 13;
		break;
	case 1008:	num_global_params = 13;
		globals->uses_dam = 0;
		globals->num_params = 12;
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
		globals->num_disk_params = 7;
		globals->convertarea_flag = 0;
		globals->num_forcings = 3;
		globals->min_error_tolerances = 13;
		break;
	case 1009:	num_global_params = 13;
		globals->uses_dam = 0;
		globals->num_params = 12;
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
		globals->num_disk_params = 7;
		globals->convertarea_flag = 0;
		globals->num_forcings = 3;
		globals->min_error_tolerances = 13;
		break;
	case 1010:	num_global_params = 13;
		globals->uses_dam = 0;
		globals->num_params = 12;
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
		globals->num_disk_params = 7;
		globals->convertarea_flag = 0;
		globals->num_forcings = 3;
		globals->min_error_tolerances = 13;
		break;
	case 1011:	num_global_params = 13;
		globals->uses_dam = 0;
		globals->num_params = 12;
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
		globals->num_disk_params = 7;
		globals->convertarea_flag = 0;
		globals->num_forcings = 3;
		globals->min_error_tolerances = 17;
		break;
	case 1012:	num_global_params = 13;
		globals->uses_dam = 0;
		globals->num_params = 12;
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
		globals->num_disk_params = 7;
		globals->convertarea_flag = 0;
		globals->num_forcings = 3;
		globals->min_error_tolerances = 13;
		break;
	case 1013:	num_global_params = 13;
		globals->uses_dam = 0;
		globals->num_params = 12;
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
		globals->num_disk_params = 7;
		globals->convertarea_flag = 0;
		globals->num_forcings = 3;
		globals->min_error_tolerances = 13;
		break;
	case 10000:	num_global_params = 13;
		globals->uses_dam = 0;
		globals->num_params = 12;
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
		globals->num_disk_params = 7;
		globals->convertarea_flag = 0;
		globals->num_forcings = 3;
		globals->min_error_tolerances = 106;
		break;
	case 10001:	num_global_params = 13;
		globals->uses_dam = 0;
		globals->num_params = 13;
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
		globals->num_disk_params = 8;
		globals->convertarea_flag = 0;
		globals->num_forcings = 3;
		globals->min_error_tolerances = 16;
		break;
	case 10002:	num_global_params = 13;
		globals->uses_dam = 0;
		globals->num_params = 13;
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
		globals->num_disk_params = 8;
		globals->convertarea_flag = 0;
		globals->num_forcings = 3;
		globals->min_error_tolerances = 16;
		break;
	case 10005:	num_global_params = 13;
		globals->uses_dam = 0;
		globals->num_params = 13;
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
		globals->num_disk_params = 8;
		globals->convertarea_flag = 0;
		globals->num_forcings = 3;
		globals->min_error_tolerances = 16;
		break;
		//--------------------------------------------------------------------------------------------
	default:	printf("Error: Invalid model_uid (%u) in SetParamSizes.\n", model_uid);
		MPI_Abort(MPI_COMM_WORLD, 1);
		//--------------------------------------------------------------------------------------------
	}

	//Make sure the appropriate number of global parameters are given
	if (globals->num_global_params < num_global_params)
	{
		printf("\nError: Obtained %u parameters from .gbl file. Expected %u for model model_uid %hu.\n", globals->num_global_params, num_global_params, model_uid);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	if (globals->num_global_params > num_global_params)
		printf("\nWarning: Obtained %u parameters from .gbl file. Expected %u for model model_uid %hu.\n", globals->num_global_params, num_global_params, model_uid);
}



//Sets the function to be used when writing outputs. This method should set the following field:
//output_constrains_hdf5
//output_constrains_psql
//output_constrains_rec
void SetOutputConstraints(GlobalVars* globals)
{
    unsigned short int model_uid = globals->model_uid;
    //Set dim and start of differential variables
    switch (model_uid)
    {
        //--------------------------------------------------------------------------------------------
        case 196:
            globals->OutputConstrainsHdf5 = &OutputConstraints_Model196_Hdf5;
            globals->OutputConstrainsPsql = NULL;
            globals->OutputConstrainsRec = NULL;
            break;
        case 254:
            globals->OutputConstrainsHdf5 = &OutputConstraints_Model254_Hdf5;
            globals->OutputConstrainsPsql = NULL;
            globals->OutputConstrainsRec = NULL;
			break;
        case 256:
            globals->OutputConstrainsHdf5 = &OutputConstraints_Model256_Hdf5;
            globals->OutputConstrainsPsql = NULL;
            globals->OutputConstrainsRec = NULL;
            break;
        default:
            globals->OutputConstrainsHdf5 = NULL;
            globals->OutputConstrainsPsql = NULL;
            globals->OutputConstrainsRec = NULL;
            break;
    }
}


//Performs some unit conversions on the data in params. This takes place immediately after reading in the DEM data,
//so these changes are available in all routines of definetype.c, if params is available. Note that dam data
//and precalculations are not available here.
//VEC* params:		Vector of parameters to convert.
//unsigned int model_uid:	The index of the model.
void ConvertParams(
	double *params,
	unsigned int model_uid,
	void* external)
{
	if (model_uid == 19)
	{
		params[1] *= 1000;	//L: km -> m
		params[2] *= 1e6;	//A_h: km^2 -> m^2
	}
	else if (model_uid == 190 || model_uid == 191 || model_uid == 192 || model_uid == 195 || model_uid == 196)
	{
		params[1] *= 1000;	//L: km -> m
		params[2] *= 1e6;	//A_h: km^2 -> m^2
	}
	else if (model_uid == 20)
	{
		//params[0] *= 1e6;	//km^2 -> m^2
		params[1] *= 1000;	//km -> m
		params[2] *= 1e6;	//km^2 -> m^2
	}
	else if (model_uid == 60)
	{
		params[1] *= 1000;	//L: km -> m
		params[2] *= 1e6;	//A_h: km^2 -> m^2
		params[3] *= 1e-3;	//h_b: mm->m
	}
	else if (model_uid == 21)
	{
		//params[0] *= 1e6;	//km^2 -> m^2
		params[1] *= 1000;	//km -> m
		params[2] *= 1e6;	//km^2 -> m^2
	}
	else if (model_uid == 22 || model_uid == 23 || model_uid == 40)
	{
		//params[0] *= 1e6;	//km^2 -> m^2
		params[1] *= 1000;	//km -> m
		params[2] *= 1e6;	//km^2 -> m^2
	}
	else if (model_uid <= 5)
	{
		params[0] *= 1000;	//km -> m
		params[3] *= .001;	//mm -> m
		params[4] *= .001;	//mm -> m
	}
	else if (model_uid == 6)
	{
		params[0] *= 1000;	//km -> m
		params[3] *= .001;	//mm -> m
	}
	else if (model_uid == 15 || model_uid == 315)
	{
		params[0] *= 1000;	//L: km -> m
		params[3] *= .001;	//h_b: mm -> m
		params[4] *= .001;	//h_H: mm -> m
	}
	else if (model_uid == 30)
	{
		params[0] *= 1000;		//L_h:  km -> m
		params[4] *= .001;		//H_h:  mm -> m
		params[5] *= 1000.0; 	//MaxInfRate:  m/hr -> mm/hr
	}
	else if (model_uid == 105)
	{
		params[0] *= 1000;	//km -> m
		params[3] *= .001;	//mm -> m
		params[4] *= .001;	//mm -> m
	}
	else if (model_uid == 200)	//!!!! Did I screw these up on accident?!!!!
	{
		params[0] *= 1000;	//L_h:  km -> m
							/*
							//params[3] *= .001;	//mm -> m
							params[4] *= .001;	//H_h:  mm -> m
							params[5] *= 1000.0; //MaxInfRate:  m/hr -> mm/hr
							*/
	}
	else if (model_uid == 219)
	{
		params[1] *= 1000;	//L: km -> m
		params[2] *= 1e6;	//A_h: km^2 -> m^2
	}
	else if (model_uid == 250)
	{
		params[1] *= 1000;		//L_h: km -> m
		params[2] *= 1e6;		//A_h: km^2 -> m^2
		params[4] *= .001;		//H_h: mm -> m
	}
	else if (model_uid == 10005 || model_uid == 10002 || model_uid == 1011 || model_uid == 10001 || model_uid==1009 || model_uid == 10000 || model_uid == 1008 || model_uid == 1010 || model_uid == 252 || model_uid == 253 || model_uid == 254 || model_uid == 255 || model_uid == 256 || model_uid == 257 || model_uid == 258 || model_uid == 259 || model_uid == 260 || model_uid == 261 || model_uid == 262 || model_uid == 263)
	{
		params[1] *= 1000;		//L_h: km -> m
		params[2] *= 1e6;		//A_h: km^2 -> m^2
	}
	else if (model_uid == 300 || model_uid == 301)
	{
		params[0] *= 1000;	//km -> m
		params[3] *= .001;	//mm -> m
		params[4] *= .001;	//mm -> m
	}
	else if (model_uid == 2000)
	{
		params[0] *= 1000;	//km -> m
		params[3] *= .001;	//mm -> m
		params[4] *= .001;	//mm -> m
	}
}


//Sets the system of ODEs and the Runge-Kutta solver for link. This method MUST set both link->differential
//	and link->solver. The Jacobian of f (link->jacobian) may be set here, if using an
//	implicit solver.
//Link* link: 		The link at which the ODEs and Runge-Kutta solver are selected.
//unsigned int model_uid: 	The index of the model to be set.
//unsigned int exp_imp: 0 if using an explicit solver, 1 if implicit.
//unsigned int dam: 	0 if no dam is present at link, 1 if a dam is present.
void InitRoutines(
	Link* link,
	unsigned int model_uid,
	unsigned int exp_imp,
	unsigned short dam,
	void* external)
{
	//Select appropriate RK Solver for the numerical method (link->solver)
	if ((model_uid == 21 || model_uid == 22 || model_uid == 23 || model_uid == 40 || model_uid == 261 || model_uid == 262) && dam == 1)
		link->solver = &ExplicitRKIndex1SolverDam;
	else if ((model_uid == 21 || model_uid == 22 || model_uid == 23 || model_uid == 40 || model_uid == 261 || model_uid == 262) && dam == 0)
		link->solver = &ExplicitRKIndex1Solver;
	else if (exp_imp == 0)
		link->solver = &ExplicitRKSolver;
	//	else if(link->method->exp_imp == 1)
	//		link->solver = &RadauRKSolver;
	else
	{
		printf("Warning: No solver selected for link ID %u.\n", link->ID);
	}
	//Select the RHS function of the ODE (link->differential, link->jacobian)
	if (model_uid == 252)
	{
		link->dim = 4;
		link->no_ini_start = link->dim;
		link->diff_start = 0;

		link->num_dense = 1;
		link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
		link->dense_indices[0] = 0;

		link->differential = &TopLayerHillslope;
		link->algebraic = NULL;
		link->check_state = NULL;
		link->check_consistency = &CheckConsistency_Nonzero_4States;
	}
	else if (model_uid == 254)
	{
		link->dim = 7;
		link->no_ini_start = 4;
		link->diff_start = 0;

		link->num_dense = 2;
		link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
		link->dense_indices[0] = 0;
		link->dense_indices[1] = 6;

		if (link->has_res)
		{
			link->differential = &TopLayerHillslope_Reservoirs;
			link->solver = &ForcedSolutionSolver;
		}
		else			
			link->differential = &TopLayerHillslope_extras;
		link->algebraic = NULL;
		link->check_state = NULL;
		link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
	}
	else if (model_uid == 1)
	{
		link->dim = 14;
		link->no_ini_start = link->dim;
		link->diff_start = 0;

		link->num_dense = 2;
		link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
		link->dense_indices[0] = 0;
		link->dense_indices[1] = 13;

		if (link->has_res)
		{
			link->differential = &TopLayerHillslope_Reservoirs;
			link->solver = &ForcedSolutionSolver;
		}
		else			
		link->differential = &Soil_moisture_velocity;
		link->algebraic = NULL;
		link->check_state = NULL;
		link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
	}
	else if (model_uid == 8)
	{
		link->dim = 14;
		link->no_ini_start = 11;
		link->diff_start = 0;

		link->num_dense = 2;
		link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
		link->dense_indices[0] = 0;
		link->dense_indices[1] = 13;

		if (link->has_res)
		{
			link->differential = &TopLayerHillslope_Reservoirs;
			link->solver = &ForcedSolutionSolver;
		}
		else
			link->differential = &Eight_TopLayers;
		link->algebraic = NULL;
		link->check_state = NULL;
		link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
	}
	else if (model_uid == 1004)
	{
		link->dim = 10;
		link->no_ini_start = link->dim;
		link->diff_start = 0;

		link->num_dense = 2;
		link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
		link->dense_indices[0] = 0;
		link->dense_indices[1] = 9;

		if (link->has_res)
		{
			link->differential = &TopLayerHillslope_Reservoirs;
			link->solver = &ForcedSolutionSolver;
		}
		else
		link->differential = &Four_TopLayers;
		link->algebraic = NULL;
		link->check_state = NULL;
		link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
	}
	else if (model_uid == 1005)
	{
		link->dim = 16;
		link->no_ini_start = 16;
		link->diff_start = 0;

		link->num_dense = 2;
		link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
		link->dense_indices[0] = 0;
		link->dense_indices[1] = 15;

		if (link->has_res)
		{
			link->differential = &TopLayerHillslope_Reservoirs;
			link->solver = &ForcedSolutionSolver;
		}
		else
			link->differential = &GAMPT_extras;
		link->algebraic = NULL;
		link->check_state = NULL;
		link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
	}
	else if (model_uid == 1006)
	{
		link->dim = 16;
		link->no_ini_start = 16;
		link->diff_start = 0;

		link->num_dense = 2;
		link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
		link->dense_indices[0] = 0;
		link->dense_indices[1] = 15;

		if (link->has_res)
		{
			link->differential = &TopLayerHillslope_Reservoirs;
			link->solver = &ForcedSolutionSolver;
		}
		else
			link->differential = &navid;
		link->algebraic = NULL;
		link->check_state = NULL;
		link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
	}
	else if (model_uid == 1007)
	{
		link->dim = 16;
		link->no_ini_start = 16;
		link->diff_start = 0;

		link->num_dense = 2;
		link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
		link->dense_indices[0] = 0;
		link->dense_indices[1] = 15;

		if (link->has_res)
		{
			link->differential = &TopLayerHillslope_Reservoirs;
			link->solver = &ForcedSolutionSolver;
		}
		else
			link->differential = &navid0;
		link->algebraic = NULL;
		link->check_state = NULL;
		link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
	}
	else if (model_uid == 1008)
	{
		link->dim = 16;
		link->no_ini_start = 16;
		link->diff_start = 0;

		link->num_dense = 2;
		link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
		link->dense_indices[0] = 0;
		link->dense_indices[1] = 15;

		if (link->has_res)
		{
			link->differential = &TopLayerHillslope_Reservoirs;
			link->solver = &ForcedSolutionSolver;
		}
		else
			link->differential = &navid1;
		link->algebraic = NULL;
		link->check_state = NULL;
		link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
	}
	else if (model_uid == 1009)
	{
		link->dim = 16;
		link->no_ini_start = 16;
		link->diff_start = 0;

		link->num_dense = 2;
		link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
		link->dense_indices[0] = 0;
		link->dense_indices[1] = 15;

		if (link->has_res)
		{
			link->differential = &TopLayerHillslope_Reservoirs;
			link->solver = &ForcedSolutionSolver;
		}
		else
			link->differential = &navid2;
		link->algebraic = NULL;
		link->check_state = NULL;
		link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
	}
	else if (model_uid == 1010)
	{
		link->dim = 16;
		link->no_ini_start = 16;
		link->diff_start = 0;

		link->num_dense = 2;
		link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
		link->dense_indices[0] = 0;
		link->dense_indices[1] = 15;

		if (link->has_res)
		{
			link->differential = &TopLayerHillslope_Reservoirs;
			link->solver = &ForcedSolutionSolver;
		}
		else
			link->differential = &navid_gampt;
		link->algebraic = NULL;
		link->check_state = NULL;
		link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
	}
	else if (model_uid == 1011)
	{
		link->dim = 17;
		link->no_ini_start = 17;
		link->diff_start = 0;

		link->num_dense = 2;
		link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
		link->dense_indices[0] = 0;
		link->dense_indices[1] = 16;

		if (link->has_res)
		{
			link->differential = &TopLayerHillslope_Reservoirs;
			link->solver = &ForcedSolutionSolver;
		}
		else
			link->differential = &herrada_gampt;
		link->algebraic = NULL;
		link->check_state = NULL;
		link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
	}
	else if (model_uid == 1012)
	{
		link->dim = 16;
		link->no_ini_start = 16;
		link->diff_start = 0;

		link->num_dense = 2;
		link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
		link->dense_indices[0] = 0;
		link->dense_indices[1] = 15;

		if (link->has_res)
		{
			link->differential = &TopLayerHillslope_Reservoirs;
			link->solver = &ForcedSolutionSolver;
		}
		else
			link->differential = &navid_gampt1012;
		link->algebraic = NULL;
		link->check_state = NULL;
		link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
	}
	else if (model_uid == 1013)
	{
		link->dim = 17;
		link->no_ini_start = 17;
		link->diff_start = 0;

		link->num_dense = 2;
		link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
		link->dense_indices[0] = 0;
		link->dense_indices[1] = 16;

		if (link->has_res)
		{
			link->differential = &TopLayerHillslope_Reservoirs;
			link->solver = &ForcedSolutionSolver;
		}
		else
			link->differential = &herrada_gampt1013;
		link->algebraic = NULL;
		link->check_state = NULL;
		link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
	}
	else if (model_uid == 10000)
	{
		link->dim = 106;
		link->no_ini_start = 106;
		link->diff_start = 0;

		link->num_dense = 2;
		link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
		link->dense_indices[0] = 0;
		link->dense_indices[1] = 105;

		if (link->has_res)
		{
			link->differential = &TopLayerHillslope_Reservoirs;
			link->solver = &ForcedSolutionSolver;
		}
		else
			link->differential = &navid_100layer;
		link->algebraic = NULL;
		link->check_state = NULL;
		link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
	}
	else if (model_uid == 10001)
	{
		link->dim = 16;
		link->no_ini_start = 16;
		link->diff_start = 0;

		link->num_dense = 2;
		link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
		link->dense_indices[0] = 0;
		link->dense_indices[1] = 15;

		if (link->has_res)
		{
			link->differential = &TopLayerHillslope_Reservoirs;
			link->solver = &ForcedSolutionSolver;
		}
		else
			link->differential = &navid_10layer;
		link->algebraic = NULL;
		link->check_state = NULL;
		link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
	}
	else if (model_uid == 10002)
	{
		link->dim = 16;
		link->no_ini_start = 16;
		link->diff_start = 0;

		link->num_dense = 2;
		link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
		link->dense_indices[0] = 0;
		link->dense_indices[1] = 15;

		if (link->has_res)
		{
			link->differential = &TopLayerHillslope_Reservoirs;
			link->solver = &ForcedSolutionSolver;
		}
		else
			link->differential = &navid_10layer_new;
		link->algebraic = NULL;
		link->check_state = NULL;
		link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
	}
	else if (model_uid == 10005)
	{
		link->dim = 16;
		link->no_ini_start = 16;
		link->diff_start = 0;

		link->num_dense = 2;
		link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
		link->dense_indices[0] = 0;
		link->dense_indices[1] = 15;

		if (link->has_res)
		{
			link->differential = &TopLayerHillslope_Reservoirs;
			link->solver = &ForcedSolutionSolver;
		}
		else
			link->differential = &navid_extra_water;
		link->algebraic = NULL;
		link->check_state = NULL;
		link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
	}
	else {
		printf("Warning: No ODE selected for link ID %u.\n", link->ID);
	}
}

	//Perform precalculations needed for the differential equation.  These should be stored in params after the DEM
	//	data and after the dam data (i.e. params[disk_params] is the first precalcuation, params[params_size]
	//	is the first dam datum). This method is run at each link. This method can do nothing, if no precalculations are
	//	expected for the model. This assumes the same number of precalculations regardless if there is a dam or not.
	//VEC* global_params:		Vector of the global parameters of the system. These are already set and are available for use.
	//VEC* params:			Vector of the parameters at a link. The DEM data and dam data are already set (but may be
	//					altered here). Only the entries for precalculations need to be set.
	//unsigned int disk_params:	The first entry of params that should be set here.
	//unsigned int params_size:	First entry of the dam data. Don't change this entry or later unless you want to modify the dam!
	//unsigned int model_uid:		The index of the model.
	void Precalculations(
		Link* link_i,
		double *global_params, unsigned int num_global_params,
		double *params, unsigned int num_disk_params, unsigned int num_params,
		unsigned short dam,
		unsigned int model_uid,
		void* external)
	{
		if (model_uid == 254 || model_uid == 8 )
		{
			//Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
			//The numbering is:	0   1   2    3     4   5   6   7 
			//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent,v_B
			//The numbering is:        0      1        2     3   4     5        6   7  8 9  10       11
			double* vals = params;
			double A_i = params[0];
			double L_i = params[1];
			double A_h = params[2];

			double v_0 = global_params[0];
			double lambda_1 = global_params[1];
			double lambda_2 = global_params[2];
			double v_h = global_params[3];
			double k_i_factor = global_params[5];

			vals[3] = 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
			vals[4] = v_h * L_i / A_h * 60.0;	//[1/min] k_2
			vals[5] = vals[4] * k_i_factor;	//[1/min] k_i
			vals[6] = (0.001 / 60.0);		//(mm/hr->m/min)  c_1
			vals[7] = A_h / 60.0;	//  c_2
		}
		else if (model_uid == 1004) {
			//Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
			//The numbering is:	0   1   2    3     4   5   6   7 
			//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent,v_B
			//The numbering is:        0      1        2     3   4     5        6   7  8 9  10       11
			double* vals = params;
			double A_i = params[0];
			double L_i = params[1];
			double A_h = params[2];
			double kSat = params[3];
			double n_VG = params[4];

			double v_0 = global_params[0];
			double lambda_1 = global_params[1];
			double lambda_2 = global_params[2];
			double v_h = global_params[3];
			double k_i_factor = global_params[5];
			vals[3] = kSat; // K_sat (cm/hr -> m/min) * 0.01 / 60.0
			vals[4] = (1 - 1 / n_VG); // convert Van Genuchten n to m m = 1-1/n
			vals[7] = 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
			vals[8] = v_h * L_i / A_h * 60.0;	//[1/min] k_2
			vals[9] = vals[8] * k_i_factor;	//[1/min] k_i
			vals[10] = (0.001 / 60.0);		//(mm/hr->m/min)  c_1
			vals[11] = A_h / 60.0;	//  c_2
			
			

		}
		else if (model_uid == 1005 || model_uid == 1006 || model_uid == 1007 || model_uid == 1008 || model_uid == 1009) {
			//Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
			//The numbering is:	0   1   2    3     4   5   6   7 
			//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent,v_B
			//The numbering is:        0      1        2     3   4     5        6   7  8 9  10       11
			// Spatially distributed parameters
			double* vals = params;
			double A_i = params[0]; // Upstream Area [m^2]
			double L_i = params[1]; // Link length [m]
			double A_h = params[2]; // Hillslope Area [m^2]
			double kSat = params[3]; // saturated hydraulic conductivity [cm/hr]
			double BC_lambda = params[4]; // Brooks-Corey's lambda [-]
			double theta_s = params[5]; // Saturated water content [-]
			double theta_r = params[6]; // Residual water content [-]

			double v_0 = global_params[0];
			double lambda_1 = global_params[1];
			double lambda_2 = global_params[2];
			double v_h = global_params[3];
			double k_i_factor = global_params[5]; // 
			double psi_Sat = global_params[12]; // TODO: make it spatially variable 


			//vals[3] = pow(10, kSat)  * 0.01 / 60.0; // K_sat (cm/hr -> m/min)
			//vals[3] = 10.0 / 24.0  * 0.01 / 60.0; // K_sat (cm/hr -> m/min)
			vals[3] = 0.5  * 0.01 / 60.0; // K_sat (cm/hr -> m/min)
			vals[4] = 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
			vals[5] = v_h * L_i / A_h * 60.0;	//[1/min] k_2
			vals[6] = vals[5] * k_i_factor;	//[1/min] k_i
			vals[7] = (0.001 / 60.0);		//(mm/hr->m/min)  c_1
			vals[8] = A_h / 60.0;	//  c_2
			vals[9] = theta_s;
			vals[10] = theta_r;
			vals[11] = BC_lambda;
		}
		else if (model_uid == 1010) {
			//Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
			//The numbering is:	0   1   2    3     4   5   6   7 
			//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent,v_B
			//The numbering is:        0      1        2     3   4     5        6   7  8 9  10       11
			// Spatially distributed parameters
			double* vals = params;
			double A_i = params[0]; // Upstream Area [m^2]
			double L_i = params[1]; // Link length [m]
			double A_h = params[2]; // Hillslope Area [m^2]
			double kSat = params[3]; // saturated hydraulic conductivity [cm/hr]
			double BC_lambda = params[4]; // Brooks-Corey's lambda [-]
			double theta_s = params[5]; // Saturated water content [-]
			double theta_r = params[6]; // Residual water content [-]

			double v_0 = global_params[0];
			double lambda_1 = global_params[1];
			double lambda_2 = global_params[2];
			double v_h = global_params[3];
			double k_i_factor = global_params[5]; // 
			double psi_Sat = global_params[12]; // TODO: make it spatially variable 


												//vals[3] = pow(10, kSat)  * 0.01 / 60.0; // K_sat (cm/hr -> m/min)
												//vals[3] = 10.0 / 24.0  * 0.01 / 60.0; // K_sat (cm/hr -> m/min)
			vals[3] = 0.75  * 0.01 / 60.0; // K_sat (cm/hr -> m/min)
			vals[4] = 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
			vals[5] = v_h * L_i / A_h * 60.0;	//[1/min] k_2
			vals[6] = vals[5] * k_i_factor;	//[1/min] k_i
			vals[7] = (0.001 / 60.0);		//(mm/hr->m/min)  c_1
			vals[8] = A_h / 60.0;	//  c_2
			vals[9] = theta_s;
			vals[10] = theta_r;
			vals[11] = BC_lambda;
		}
		else if (model_uid == 1011) {
			//Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
			//The numbering is:	0   1   2    3     4   5   6   7 
			//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent,v_B
			//The numbering is:        0      1        2     3   4     5        6   7  8 9  10       11
			// Spatially distributed parameters
			double* vals = params;
			double A_i = params[0]; // Upstream Area [m^2]
			double L_i = params[1]; // Link length [m]
			double A_h = params[2]; // Hillslope Area [m^2]
			double kSat = params[3]; // saturated hydraulic conductivity [cm/hr]
			double BC_lambda = params[4]; // Brooks-Corey's lambda [-]
			double theta_s = 0.42;// params[5]; // Saturated water content [-]
			double theta_r = 0.1;// params[6]; // Residual water content [-]

			double v_0 = global_params[0];
			double lambda_1 = global_params[1];
			double lambda_2 = global_params[2];
			double v_h = global_params[3];
			double k_i_factor = global_params[5]; // 
			double psi_Sat = global_params[12]; // TODO: make it spatially variable 


												//vals[3] = pow(10, kSat)  * 0.01 / 60.0; // K_sat (cm/hr -> m/min)
												//vals[3] = 10.0 / 24.0  * 0.01 / 60.0; // K_sat (cm/hr -> m/min)
			vals[3] = pow(10, kSat)  * 0.01 / 60.0; // K_sat (cm/hr -> m/min)
			vals[4] = 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
			vals[5] = v_h * L_i / A_h * 60.0;	//[1/min] k_2
			vals[6] = vals[5] * k_i_factor;	//[1/min] k_i
			vals[7] = (0.001 / 60.0);		//(mm/hr->m/min)  c_1
			vals[8] = A_h / 60.0;	//  c_2
			vals[9] = theta_s;
			vals[10] = theta_r;
			vals[11] = BC_lambda;
		}
		else if (model_uid == 1012) {
			//Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
			//The numbering is:	0   1   2    3     4   5   6   7 
			//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent,v_B
			//The numbering is:        0      1        2     3   4     5        6   7  8 9  10       11
			// Spatially distributed parameters
			double* vals = params;
			double A_i = params[0]; // Upstream Area [m^2]
			double L_i = params[1]; // Link length [m]
			double A_h = params[2]; // Hillslope Area [m^2]
			double kSat = params[3]; // saturated hydraulic conductivity [cm/hr]
			double BC_lambda = params[4]; // Brooks-Corey's lambda [-]
			double theta_s = params[5]; // Saturated water content [-]
			double theta_r = params[6]; // Residual water content [-]

			double v_0 = global_params[0];
			double lambda_1 = global_params[1];
			double lambda_2 = global_params[2];
			double v_h = global_params[3];
			double k_i_factor = global_params[5]; // 
			double psi_Sat = global_params[12]; // TODO: make it spatially variable 


												//vals[3] = pow(10, kSat)  * 0.01 / 60.0; // K_sat (cm/hr -> m/min)
												//vals[3] = 10.0 / 24.0  * 0.01 / 60.0; // K_sat (cm/hr -> m/min)
			vals[3] = 0.75  * 0.01 / 60.0; // K_sat (cm/hr -> m/min)
			vals[4] = 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
			vals[5] = v_h * L_i / A_h * 60.0;	//[1/min] k_2
			vals[6] = vals[5] * k_i_factor;	//[1/min] k_i
			vals[7] = (0.001 / 60.0);		//(mm/hr->m/min)  c_1
			vals[8] = A_h / 60.0;	//  c_2
			vals[9] = theta_s;
			vals[10] = theta_r;
			vals[11] = BC_lambda;
		}
		else if (model_uid == 1013) {
			//Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
			//The numbering is:	0   1   2    3     4   5   6   7 
			//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent,v_B
			//The numbering is:        0      1        2     3   4     5        6   7  8 9  10       11
			// Spatially distributed parameters
			double* vals = params;
			double A_i = params[0]; // Upstream Area [m^2]
			double L_i = params[1]; // Link length [m]
			double A_h = params[2]; // Hillslope Area [m^2]
			double kSat = params[3]; // saturated hydraulic conductivity [cm/hr]
			double BC_lambda = params[4]; // Brooks-Corey's lambda [-]
			double theta_s = params[5]; // Saturated water content [-]
			double theta_r = params[6]; // Residual water content [-]

			double v_0 = global_params[0];
			double lambda_1 = global_params[1];
			double lambda_2 = global_params[2];
			double v_h = global_params[3];
			double k_i_factor = global_params[5]; // 
			double psi_Sat = global_params[12]; // TODO: make it spatially variable 


												//vals[3] = pow(10, kSat)  * 0.01 / 60.0; // K_sat (cm/hr -> m/min)
												//vals[3] = 10.0 / 24.0  * 0.01 / 60.0; // K_sat (cm/hr -> m/min)
			vals[3] = 0.75  * 0.01 / 60.0; // K_sat (cm/hr -> m/min)
			vals[4] = 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
			vals[5] = v_h * L_i / A_h * 60.0;	//[1/min] k_2
			vals[6] = vals[5] * k_i_factor;	//[1/min] k_i
			vals[7] = (0.001 / 60.0);		//(mm/hr->m/min)  c_1
			vals[8] = A_h / 60.0;	//  c_2
			vals[9] = theta_s;
			vals[10] = theta_r;
			vals[11] = BC_lambda;
		}
		else if (model_uid == 10000) {
			//Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
			//The numbering is:	0   1   2    3     4   5   6   7 
			//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent,v_B
			//The numbering is:        0      1        2     3   4     5        6   7  8 9  10       11
			// Spatially distributed parameters
			double* vals = params;
			double A_i = params[0]; // Upstream Area [m^2]
			double L_i = params[1]; // Link length [m]
			double A_h = params[2]; // Hillslope Area [m^2]
			double kSat = params[3]; // saturated hydraulic conductivity [cm/hr]
			double BC_lambda = params[4]; // Brooks-Corey's lambda [-]
			double theta_s = params[5]; // Saturated water content [-]
			double theta_r = params[6]; // Residual water content [-]

			double v_0 = global_params[0];
			double lambda_1 = global_params[1];
			double lambda_2 = global_params[2];
			double v_h = global_params[3];
			double k_i_factor = global_params[5]; // 
			double psi_Sat = global_params[12]; // TODO: make it spatially variable 


												//vals[3] = pow(10, kSat)  * 0.01 / 60.0; // K_sat (cm/hr -> m/min)
												//vals[3] = 10.0 / 24.0  * 0.01 / 60.0; // K_sat (cm/hr -> m/min)
			vals[3] = pow(10.0, kSat)  * 0.01 / 60.0; // K_sat (cm/hr -> m/min)
			vals[4] = 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
			vals[5] = v_h * L_i / A_h * 60.0;	//[1/min] k_2
			vals[6] = vals[5] * k_i_factor;	//[1/min] k_i
			vals[7] = (0.001 / 60.0);		//(mm/hr->m/min)  c_1
			vals[8] = A_h / 60.0;	//  c_2
			vals[9] = theta_s;
			vals[10] = theta_r;
			vals[11] = BC_lambda;
		}
		else if (model_uid == 10001) {
			double* vals = params;
			double A_i = params[0];
			double L_i = params[1];
			double A_h = params[2];

			double kSat = params[3]; // saturated hydraulic conductivity [cm/hr]
			double BC_lambda = params[4]; // Brooks-Corey's lambda [-]
			double theta_s = params[5]; // Saturated water content [-]
			double theta_r = params[6]; // Residual water content [-]
			double psi_sat = params[7]; // Residual water content [-]

			double v_0 = global_params[0];
			double lambda_1 = global_params[1];
			double lambda_2 = global_params[2];
			double v_h = global_params[3];
			double k_i_factor = global_params[5];

			vals[3] = 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	 // [1/min]  invtau
			vals[4] = v_h * L_i / A_h * 60.0;                                // [1/min] k_2
			vals[5] = vals[4] * k_i_factor;                                  // [1/min] k_i
			vals[6] = (0.001 / 60.0);                                        // (mm/hr->m/min)  c_1
			vals[7] = A_h / 60.0;                                            // c_2
			vals[8] = pow(10.0, kSat)  * 0.01 / 60.0; // K_sat (cm/hr -> m/min)
			vals[9] = theta_s;
			vals[10] = theta_r;
			vals[11] = BC_lambda;
			vals[12] = -1 *  pow(10, psi_sat) * 10.0 * 0.01;


		}
		else if (model_uid == 10002 || model_uid ==10005) {
			double* vals = params;
			double A_i = params[0];
			double L_i = params[1];
			double A_h = params[2];

			double kSat = params[3]; // saturated hydraulic conductivity [cm/hr]
			double BC_lambda = params[4]; // Brooks-Corey's lambda [-]
			double theta_s = params[5];// < 0.35 ? 0.35 : params[5]; // Saturated water content [-]
			theta_s = theta_s < 0.25 ? 0.40 : theta_s;
			double theta_r = params[6] > 0.15 ? 0.15 : params[6]; // Residual water content [-]
			double psi_sat = params[7]; // Residual water content [-]

			double v_0 = global_params[0];
			double lambda_1 = global_params[1];
			double lambda_2 = global_params[2];
			double v_h = global_params[3];
			double k_i_factor = global_params[5];

			vals[3] = 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	 // [1/min]  invtau
			vals[4] = v_h * L_i / A_h * 60.0;                                // [1/min] k_2
			vals[5] = vals[4] * k_i_factor;                                  // [1/min] k_i
			vals[6] = (0.001 / 60.0);                                        // (mm/hr->m/min)  c_1
			vals[7] = A_h / 60.0;                                            // c_2
			vals[8] = pow(10.0, kSat)  * 0.01 / 60.0; // K_sat (cm/hr -> m/min)
			vals[9] = theta_s;
			vals[10] = theta_r;
			vals[11] = BC_lambda;
			vals[12] = -1 * pow(10, psi_sat) * 10.0 * 0.01;
		}
	}


	//Set the initial condition for the differential equations a link. This method will be called once for each link. The differential
	//	variables from any .ini, .uini, or database are already set here. Precalculations have also been made previously. All
	//	algebraic variables MUST be set here. Other initial conditions may be set here as well, depending on the model.
	//VEC* global_params:	The vector of global parameters for the system.
	//VEC* params:		The parameters for this link. DEM, dam parameters, and precalculations are available.
	//unsigned int dam:	1 if a dam is present at this link, 0 if no dam is present.
	//VEC* y_0:		The initial condition vector. Store the initial data here.
	//unsigned int model_uid:	The index of the model.
	//Returns the state of the solution (use 0 if state discontinuities are not a concern). !!!! Should the return value be an int? !!!!
int ReadInitData(
	double *global_params, unsigned int num_global_params,
	double *params, unsigned int num_params,
	QVSData* qvs,
	unsigned short int dam,
	double *y_0, unsigned int dim,
	unsigned int model_uid,
	unsigned int diff_start, unsigned int no_init_start,
	void* user,
	void* external)
{

	if (model_uid == 195)
	{
		//For this model_uid, the extra states need to be set (3)
		y_0[3] = 0.0;
	}
	else if (model_uid == 254)
	{
		//For this model_uid, the extra states need to be set (4,5,6)
		y_0[4] = 0.0;
		y_0[5] = 0.0;
		y_0[6] = y_0[0];
	}
	else if (model_uid == 8 || model_uid == 1)
	{
		//For this model_uid, the extra states need to be set (4,5,6)
		y_0[11] = 0.0;
		y_0[12] = 0.0;
		y_0[13] = y_0[0];
	} 
	else if (model_uid == 1004)
	{
		y_0[7] = 0.0;
		y_0[8] = 0.0;
		y_0[9] = y_0[0];
	}
	else if (model_uid == 1005 || model_uid == 1006 || model_uid == 1007 || model_uid == 1009 || model_uid == 1010 || model_uid == 1012)
	{
		y_0[13] = 0.0;
		y_0[14] = 0.0;
		y_0[15] = y_0[0];
	}
	else if (model_uid == 1011)
	{
		y_0[13] = 0.0;
		y_0[14] = 0.0;
		y_0[15] = y_0[0];
	}
	else if (model_uid == 1013)
	{
		y_0[13] = 0.0;
		y_0[14] = 0.0;
		y_0[15] = y_0[0];
	}
	else if (model_uid == 10000)
	{
		y_0[103] = 0.0;
		y_0[104] = 0.0;
		y_0[105] = y_0[0];
	}
	else if (model_uid == 10001 || model_uid == 10002 || model_uid==10005)
	{
		y_0[13] = 0.0;
		y_0[14] = 0.0;
		y_0[15] = y_0[0];
	}
	else
	{
		//If not using algebraic variables, then everything is already set
		return 0;
	}

	return 0;
}
