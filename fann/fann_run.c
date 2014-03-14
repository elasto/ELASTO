//USEFORTEST avgcond
#include "header_fann.h"
#include <stdlib.h>

void run_fann_with_scales_(float *fieldIn,
                           float *fieldOut,
                           char nameOfAnnFile[],int *numData,
                           int *numInput, int *numOutput)
{
  struct fann *ann;
  int i,j, nb;
  float* reorder, *result;
  
  /* Allocate an array for ordering datas coming from a fortran array(numData,numInput)
     because there is a stride on the first dimension */
  reorder=(float*) malloc(*numInput*sizeof(float));
  if (!reorder) {
     nb= *numInput * sizeof(float);
     printf("[ERROR] unable to allocate %d bytes\n",nb);
     exit(1);
     }
  result=(float*) malloc(*numOutput*sizeof(float));
  if (!result) {
     nb= *numOutput * sizeof(float);
     printf("[ERROR] unable to allocate %d bytes\n",nb);
     exit(1);
     }
  
  ann = fann_create_from_file ( nameOfAnnFile );
  if (! ann) {
     printf("\n[ERROR] failed to create ann structure from file [%s]\n",nameOfAnnFile);
     exit(1);
     }

  //quelques vérifications nombre entrées/sorties par rapport aux tailles des champs
  if ( ann->num_input != *numInput ) {
     printf("\n[ERROR] the number of input do not match with ann loaded\n");
     exit(1);
     }
  if ( ann->num_output != *numOutput) {
     printf("\n[ERROR] the number of Output do not match with ann loaded\n");
     exit(1);
     }

  for (i=0;i<*numData;i++)
  {
    /* put datas contiguous in a C 1D array */
    for (j=0;j<*numInput;j++) reorder[j]=fieldIn[i+j*(*numData)];
    /* Works on the data */
    result= fann_run( ann , reorder );
    /* Copy back the data in the fortran array*/
    for (j=0;j<*numOutput;j++) fieldOut[i+j*(*numData)]=result[j];
    //printf("%d inputs %f %f output %f\n",i,reorder[0],reorder[1],result[0]);
  }
  free(reorder);
  free(result);
  fann_destroy(ann);
}

void fann_destroy(struct fann *ann)
{
	if(ann == NULL) return;
	fann_safe_free(ann->weights);
	fann_safe_free(ann->connections);
	fann_safe_free(ann->first_layer->first_neuron);
	fann_safe_free(ann->first_layer);
//Pas alloue dans le cas present
//        printf("Taille de ann->output= %u\n",sizeof(ann->output));
//        free(ann->output);
//        fann_safe_free(ann->output);
	fann_safe_free(ann->train_errors);
	fann_safe_free(ann->train_slopes);
	fann_safe_free(ann->prev_train_slopes);
	fann_safe_free(ann->prev_steps);
	fann_safe_free(ann->prev_weights_deltas);
	fann_safe_free(ann->errstr);
	fann_safe_free(ann->cascade_activation_functions);
	fann_safe_free(ann->cascade_activation_steepnesses);
	fann_safe_free( ann->scale_mean_in );
	fann_safe_free( ann->scale_deviation_in );
	fann_safe_free( ann->scale_new_min_in );
	fann_safe_free( ann->scale_factor_in );
	fann_safe_free( ann->scale_mean_out );
	fann_safe_free( ann->scale_deviation_out );
	fann_safe_free( ann->scale_new_min_out );
	fann_safe_free( ann->scale_factor_out );
	fann_safe_free(ann);
}


struct fann *fann_allocate_structure(unsigned int num_layers)
{
	struct fann *ann;

	if(num_layers < 2)
	{
		return NULL;
	}

	/* allocate and initialize the main network structure */
	ann = (struct fann *) malloc(sizeof(struct fann));
	if(ann == NULL)
	{
	        printf("[ERROR] in fann N");
		return NULL;
	}

	//ann->errno_f = FANN_E_NO_ERROR;
	//ann->error_log = fann_default_error_log;
	ann->errstr = NULL;
	ann->learning_rate = 0.7f;
	ann->learning_momentum = 0.0;
	ann->total_neurons = 0;
	ann->total_connections = 0;
	ann->num_input = 0;
	ann->num_output = 0;
	ann->train_errors = NULL;
	ann->train_slopes = NULL;
	ann->prev_steps = NULL;
	ann->prev_train_slopes = NULL;
	ann->prev_weights_deltas = NULL;
	ann->training_algorithm = FANN_TRAIN_RPROP;
	ann->num_MSE = 0;
	ann->MSE_value = 0;
	ann->num_bit_fail = 0;
	ann->bit_fail_limit = (float)0.35;
	ann->network_type = FANN_NETTYPE_LAYER;
	//ann->train_error_function = FANN_ERRORFUNC_TANH;
	//ann->train_stop_function = FANN_STOPFUNC_MSE;
	//ann->callback = NULL;
        //ann->user_data = NULL; /* User is responsible for deallocation */
	ann->weights = NULL;
	ann->connections = NULL;
	ann->output = NULL;
	ann->scale_mean_in = NULL;
	ann->scale_deviation_in = NULL;
	ann->scale_new_min_in = NULL;
	ann->scale_factor_in = NULL;
	ann->scale_mean_out = NULL;
	ann->scale_deviation_out = NULL;
	ann->scale_new_min_out = NULL;
	ann->scale_factor_out = NULL;
	/* variables used for cascade correlation (reasonable defaults) */
	ann->cascade_output_change_fraction = 0.01f;
	ann->cascade_candidate_change_fraction = 0.01f;
	ann->cascade_output_stagnation_epochs = 12;
	ann->cascade_candidate_stagnation_epochs = 12;
	ann->cascade_num_candidate_groups = 2;
	ann->cascade_weight_multiplier = (float)0.4;
	ann->cascade_candidate_limit = (float)1000.0;
	ann->cascade_max_out_epochs = 150;
	ann->cascade_max_cand_epochs = 150;
	ann->cascade_candidate_scores = NULL;
	ann->cascade_activation_functions_count = 10;
	ann->cascade_activation_functions =(enum fann_activationfunc_enum *)calloc(ann->cascade_activation_functions_count, sizeof(enum fann_activationfunc_enum));
	if(ann->cascade_activation_functions == NULL)
	{
//		fann_error(NULL, FANN_E_CANT_ALLOCATE_MEM);
		free(ann);
		return NULL;
	}
							   
	ann->cascade_activation_functions[0] = FANN_SIGMOID;
	ann->cascade_activation_functions[1] = FANN_SIGMOID_SYMMETRIC;
	ann->cascade_activation_functions[2] = FANN_GAUSSIAN;
	ann->cascade_activation_functions[3] = FANN_GAUSSIAN_SYMMETRIC;
	ann->cascade_activation_functions[4] = FANN_ELLIOT;
	ann->cascade_activation_functions[5] = FANN_ELLIOT_SYMMETRIC;
	ann->cascade_activation_functions[6] = FANN_SIN_SYMMETRIC;
	ann->cascade_activation_functions[7] = FANN_COS_SYMMETRIC;
	ann->cascade_activation_functions[8] = FANN_SIN;
	ann->cascade_activation_functions[9] = FANN_COS;

	ann->cascade_activation_steepnesses_count = 4;
	ann->cascade_activation_steepnesses = 
		(float *)calloc(ann->cascade_activation_steepnesses_count, 
							   sizeof(float));
	if(ann->cascade_activation_steepnesses == NULL)
	{
		fann_safe_free(ann->cascade_activation_functions);
	        printf("[ERROR] in fann M");
		free(ann);
		return NULL;
	}
	
	ann->cascade_activation_steepnesses[0] = (float)0.25;
	ann->cascade_activation_steepnesses[1] = (float)0.5;
	ann->cascade_activation_steepnesses[2] = (float)0.75;
	ann->cascade_activation_steepnesses[3] = (float)1.0;

	/* Variables for use with with Quickprop training (reasonable defaults) */
	ann->quickprop_decay = (float) -0.0001;
	ann->quickprop_mu = 1.75;

	/* Variables for use with with RPROP training (reasonable defaults) */
	ann->rprop_increase_factor = (float) 1.2;
	ann->rprop_decrease_factor = 0.5;
	ann->rprop_delta_min = 0.0;
	ann->rprop_delta_max = 50.0;
	ann->rprop_delta_zero = 0.1;
	//fann_init_error_data((struct fann_error *) ann);

	/* allocate room for the layers */
	ann->first_layer = (struct fann_layer *) calloc(num_layers, sizeof(struct fann_layer));
	if(ann->first_layer == NULL)
	{
	        printf("[ERROR] in fann L");
		free(ann);
		return NULL;
	}

	ann->last_layer = ann->first_layer + num_layers;
	return ann;
}


void fann_set_activation_function_hidden(struct fann *ann,enum fann_activationfunc_enum activation_function)
{
	struct fann_neuron *last_neuron, *neuron_it;
	struct fann_layer *layer_it;
	struct fann_layer *last_layer = ann->last_layer - 1;	/* -1 to not update the output layer */

	for(layer_it = ann->first_layer + 1; layer_it != last_layer; layer_it++)
	{
		last_neuron = layer_it->last_neuron;
		for(neuron_it = layer_it->first_neuron; neuron_it != last_neuron; neuron_it++)
		{
			neuron_it->activation_function = activation_function;
		}
	}
}


void fann_set_activation_function_output(struct fann *ann,enum fann_activationfunc_enum activation_function)
{
	struct fann_neuron *last_neuron, *neuron_it;
	struct fann_layer *last_layer = ann->last_layer - 1;

	last_neuron = last_layer->last_neuron;
	for(neuron_it = last_layer->first_neuron; neuron_it != last_neuron; neuron_it++)
	{
		neuron_it->activation_function = activation_function;
	}
}


void fann_set_activation_steepness_output(struct fann *ann,float steepness)
{
	struct fann_neuron *last_neuron, *neuron_it;
	struct fann_layer *last_layer = ann->last_layer - 1;

	last_neuron = last_layer->last_neuron;
	for(neuron_it = last_layer->first_neuron; neuron_it != last_neuron; neuron_it++)
	{
		neuron_it->activation_steepness = steepness;
	}
}

void fann_set_activation_steepness_hidden(struct fann *ann,float steepness)
{
	struct fann_neuron *last_neuron, *neuron_it;
	struct fann_layer *layer_it;
	struct fann_layer *last_layer = ann->last_layer - 1;	/* -1 to not update the output layer */

	for(layer_it = ann->first_layer + 1; layer_it != last_layer; layer_it++)
	{
		last_neuron = layer_it->last_neuron;
		for(neuron_it = layer_it->first_neuron; neuron_it != last_neuron; neuron_it++)
		{
			neuron_it->activation_steepness = steepness;
		}
	}
}

struct fann * fann_create_from_file(const char *configuration_file)
{
	struct fann *ann;
	FILE *conf = fopen(configuration_file, "r");

	if(!conf)
	{
	        printf("[ERROR] not able to open ANN from file %s",configuration_file);
		return NULL;
	}
	//ann = fann_create_from_fd(conf, configuration_file);
	ann = fann_create_from_fd(conf);
	fclose(conf);
	return ann;
}


float * fann_run(struct fann * ann, float * input)
{
	struct fann_neuron *neuron_it, *last_neuron, *neurons, **neuron_pointers;
	unsigned int i, num_connections, num_input, num_output;
	float neuron_sum, *output;
	float *weights;
	struct fann_layer *layer_it, *last_layer;
	unsigned int activation_function;
	float steepness;

	/* store some variables local for fast access */
	struct fann_neuron *first_neuron = ann->first_layer->first_neuron;

	float max_sum;	

	/* first set the input */
	num_input = ann->num_input;
	for(i = 0; i != num_input; i++)
	{
		first_neuron[i].value = input[i];
	}
	/* Set the bias neuron in the input layer */
	(ann->first_layer->last_neuron - 1)->value = 1;
	last_layer = ann->last_layer;
	for(layer_it = ann->first_layer + 1; layer_it != last_layer; layer_it++)
	{
		last_neuron = layer_it->last_neuron;
		for(neuron_it = layer_it->first_neuron; neuron_it != last_neuron; neuron_it++)
		{
			if(neuron_it->first_con == neuron_it->last_con)
			{
				/* bias neurons */
				neuron_it->value = 1;
				continue;
			}

			activation_function = neuron_it->activation_function;
			steepness = neuron_it->activation_steepness;
			neuron_sum = 0;
			num_connections = neuron_it->last_con - neuron_it->first_con;
			weights = ann->weights + neuron_it->first_con;
			if(ann->connection_rate >= 1)
			{
				if(ann->network_type == FANN_NETTYPE_SHORTCUT)
				{
					neurons = ann->first_layer->first_neuron;
				}
				else
				{
					neurons = (layer_it - 1)->first_neuron;
				}
				/* unrolled loop start */
				i = num_connections & 3;	/* same as modulo 4 */
				switch (i)
				{
					case 3:
						neuron_sum += fann_mult(weights[2], neurons[2].value);
					case 2:
						neuron_sum += fann_mult(weights[1], neurons[1].value);
					case 1:
						neuron_sum += fann_mult(weights[0], neurons[0].value);
					case 0:
						break;
				}

				for(; i != num_connections; i += 4)
				{
					neuron_sum +=
						fann_mult(weights[i], neurons[i].value) +
						fann_mult(weights[i + 1], neurons[i + 1].value) +
						fann_mult(weights[i + 2], neurons[i + 2].value) +
						fann_mult(weights[i + 3], neurons[i + 3].value);
				}
				/* unrolled loop end */

				/*
				 * for(i = 0;i != num_connections; i++){
				 * printf("%f += %f*%f, ", neuron_sum, weights[i], neurons[i].value);
				 * neuron_sum += fann_mult(weights[i], neurons[i].value);
				 * }
				 */
			}
			else
			{
				neuron_pointers = ann->connections + neuron_it->first_con;

				i = num_connections & 3;	/* same as modulo 4 */
				switch (i)
				{
					case 3:
						neuron_sum += fann_mult(weights[2], neuron_pointers[2]->value);
					case 2:
						neuron_sum += fann_mult(weights[1], neuron_pointers[1]->value);
					case 1:
						neuron_sum += fann_mult(weights[0], neuron_pointers[0]->value);
					case 0:
						break;
				}

				for(; i != num_connections; i += 4)
				{
					neuron_sum +=
						fann_mult(weights[i], neuron_pointers[i]->value) +
						fann_mult(weights[i + 1], neuron_pointers[i + 1]->value) +
						fann_mult(weights[i + 2], neuron_pointers[i + 2]->value) +
						fann_mult(weights[i + 3], neuron_pointers[i + 3]->value);
				}
			}

			neuron_sum = fann_mult(steepness, neuron_sum);
			
			max_sum = 150/steepness;
			if(neuron_sum > max_sum)
				neuron_sum = max_sum;
			else if(neuron_sum < -max_sum)
				neuron_sum = -max_sum;
			
			neuron_it->sum = neuron_sum;

			fann_activation_switch(activation_function, neuron_sum, neuron_it->value);
		}
	}

	/* set the output */
	output = ann->output;
	num_output = ann->num_output;
	neurons = (ann->last_layer - 1)->first_neuron;
	for(i = 0; i != num_output; i++)
	{
		output[i] = neurons[i].value;
	}
	return ann->output;
}


struct fann *fann_create_from_fd(FILE * conf)
{
	unsigned int num_layers, layer_size, input_neuron, i, num_connections;
	unsigned int scale_included,trash;
	struct fann_neuron *first_neuron, *neuron_it, *last_neuron, **connected_neurons;
	float *weights;
	struct fann_layer *layer_it;
	struct fann *ann = NULL;

	char *read_version;

	read_version = (char *) calloc(strlen(FANN_CONF_VERSION "\n"), 1);
	if(read_version == NULL)
	{
		//fann_error(NULL, FANN_E_CANT_ALLOCATE_MEM);
		return NULL;
	}

	fread(read_version, 1, strlen(FANN_CONF_VERSION "\n"), conf);	/* reads version */

	/* compares the version information */
	if(strncmp(read_version, FANN_CONF_VERSION "\n", strlen(FANN_CONF_VERSION "\n")) != 0)
	{
		if(strncmp(read_version, "FANN_FLO_1.1\n", strlen("FANN_FLO_1.1\n")) == 0)
		{
			free(read_version);
			//return fann_create_from_fd_1_1(conf, configuration_file);
			return fann_create_from_fd_1_1(conf);
		}

		/* Maintain compatibility with 2.0 version that doesnt have scale parameters. */
		if(strncmp(read_version, "FANN_FLO_2.0\n", strlen("FANN_FLO_2.0\n")) != 0 &&
		   strncmp(read_version, "FANN_FLO_2.1\n", strlen("FANN_FLO_2.1\n")) != 0)
		{
			free(read_version);
	                printf("[ERROR] in fann J");
			return NULL;
		}
	}

	free(read_version);

        fann_scanf("%u", "num_layers", &num_layers);

	ann = fann_allocate_structure(num_layers);
	if(ann == NULL)
	{
		return NULL;
	}
        fann_scanf("%f", "learning_rate", &ann->learning_rate);
        fann_scanf("%f", "connection_rate", &ann->connection_rate);
        fann_scanf("%u", "network_type", (unsigned int *)&ann->network_type);
	fann_scanf("%f", "learning_momentum", &ann->learning_momentum);
	fann_scanf("%u", "training_algorithm", (unsigned int *)&ann->training_algorithm);
	//fann_scanf("%u", "train_error_function", (unsigned int *)&ann->train_error_function);
	//fann_scanf("%u", "train_stop_function", (unsigned int *)&ann->train_stop_function);
	fann_scanf("%u", "train_error_function", &trash);
	fann_scanf("%u", "train_stop_function",  &trash);
	fann_scanf("%f", "cascade_output_change_fraction", &ann->cascade_output_change_fraction);
	fann_scanf("%f", "quickprop_decay", &ann->quickprop_decay);
	fann_scanf("%f", "quickprop_mu", &ann->quickprop_mu);
	fann_scanf("%f", "rprop_increase_factor", &ann->rprop_increase_factor);
	fann_scanf("%f", "rprop_decrease_factor", &ann->rprop_decrease_factor);
	fann_scanf("%f", "rprop_delta_min", &ann->rprop_delta_min);
	fann_scanf("%f", "rprop_delta_max", &ann->rprop_delta_max);
	fann_scanf("%f", "rprop_delta_zero", &ann->rprop_delta_zero);
	fann_scanf("%u", "cascade_output_stagnation_epochs", &ann->cascade_output_stagnation_epochs);
	fann_scanf("%f", "cascade_candidate_change_fraction", &ann->cascade_candidate_change_fraction);
	fann_scanf("%u", "cascade_candidate_stagnation_epochs", &ann->cascade_candidate_stagnation_epochs);
	fann_scanf("%u", "cascade_max_out_epochs", &ann->cascade_max_out_epochs);
	fann_scanf("%u", "cascade_max_cand_epochs", &ann->cascade_max_cand_epochs);	
	fann_scanf("%u", "cascade_num_candidate_groups", &ann->cascade_num_candidate_groups);
	fann_scanf("%f", "bit_fail_limit", &ann->bit_fail_limit);
	fann_scanf("%f", "cascade_candidate_limit", &ann->cascade_candidate_limit);
	fann_scanf("%f", "cascade_weight_multiplier", &ann->cascade_weight_multiplier);
	fann_scanf("%u", "cascade_activation_functions_count", &ann->cascade_activation_functions_count);

	/* reallocate mem */
	ann->cascade_activation_functions = 
		(enum fann_activationfunc_enum *)realloc(ann->cascade_activation_functions, 
		ann->cascade_activation_functions_count * sizeof(enum fann_activationfunc_enum));
	if(ann->cascade_activation_functions == NULL)
	{
	        printf("[ERROR] in fann I");
		fann_destroy(ann);
		return NULL;
	}

	fscanf(conf, "cascade_activation_functions=");
	for(i = 0; i < ann->cascade_activation_functions_count; i++)
        {
		fscanf(conf, "%u ", (unsigned int *)&ann->cascade_activation_functions[i]);
	}
	fann_scanf("%u", "cascade_activation_steepnesses_count", &ann->cascade_activation_steepnesses_count);
	/* reallocate mem */
	ann->cascade_activation_steepnesses = 
		(float *)realloc(ann->cascade_activation_steepnesses, 
		ann->cascade_activation_steepnesses_count * sizeof(float));
	if(ann->cascade_activation_steepnesses == NULL)
	{
	        printf("[ERROR] in fann H");
		fann_destroy(ann);
		return NULL;
	}

	fscanf(conf, "cascade_activation_steepnesses=");
	for(i = 0; i < ann->cascade_activation_steepnesses_count; i++)
        {
		fscanf(conf, "%f ", &ann->cascade_activation_steepnesses[i]);
        }
	fscanf(conf, "layer_sizes=");
	/* determine how many neurons there should be in each layer */
	for(layer_it = ann->first_layer; layer_it != ann->last_layer; layer_it++)
	{
		if(fscanf(conf, "%u ", &layer_size) != 1)
		{
	                printf("[ERROR] in fann G");
			fann_destroy(ann);
			return NULL;
		}
		layer_it->first_neuron = NULL;
		layer_it->last_neuron = layer_it->first_neuron + layer_size;
		ann->total_neurons += layer_size;
	}

	ann->num_input = ann->first_layer->last_neuron - ann->first_layer->first_neuron - 1;
	ann->num_output = ((ann->last_layer - 1)->last_neuron - (ann->last_layer - 1)->first_neuron);
	if(ann->network_type == FANN_NETTYPE_LAYER)
	{
		ann->num_output--;
	}

#define SCALE_LOAD( what, where )\
	fscanf( conf, #what "_" #where "=" );\
	for(i = 0; i < ann->num_##where##put; i++)\
	{					\
		if(fscanf( conf, "%f ", (float *)&ann->what##_##where[ i ] ) != 1) \
		{								\
			/*fann_error((struct fann_error *) ann, FANN_E_CANT_READ_CONFIG, #what "_" #where, configuration_file);*/ \
			printf("[ERROR] in define SCALE_LOAD of routine fann\n");\
			fann_destroy(ann);\
			return NULL;\
		}\
	}
	
	if(fscanf(conf, "scale_included=%u\n", &scale_included) == 1 && scale_included == 1)
	{
		fann_allocate_scale(ann);
		SCALE_LOAD( scale_mean,			in )
		SCALE_LOAD( scale_deviation,	in )
		SCALE_LOAD( scale_new_min,		in )
		SCALE_LOAD( scale_factor,		in )
	
		SCALE_LOAD( scale_mean,			out )
		SCALE_LOAD( scale_deviation,	out )
		SCALE_LOAD( scale_new_min,		out )
		SCALE_LOAD( scale_factor,		out )
	}
#undef SCALE_LOAD
	
	fann_allocate_neurons(ann);
	last_neuron = (ann->last_layer - 1)->last_neuron;
	fscanf(conf, "neurons (num_inputs, activation_function, activation_steepness)=");
	for(neuron_it = ann->first_layer->first_neuron; neuron_it != last_neuron; neuron_it++)
	{
		if(fscanf(conf, "(%u, %u, %f) ", &num_connections, (unsigned int *)&neuron_it->activation_function,&neuron_it->activation_steepness) != 3)
		{
	                printf("[ERROR] in fann F");
			fann_destroy(ann);
			return NULL;
		}
 		neuron_it->first_con = ann->total_connections;
		ann->total_connections += num_connections;
		neuron_it->last_con = ann->total_connections;
	}

	fann_allocate_connections(ann);
	connected_neurons = ann->connections;
	weights = ann->weights;
	first_neuron = ann->first_layer->first_neuron;

	fscanf(conf, "connections (connected_to_neuron, weight)=");
	for(i = 0; i < ann->total_connections; i++)
	{
		if(fscanf(conf, "(%u, " "%f" ") ", &input_neuron, &weights[i]) != 2)
		{
	                printf("[ERROR] in fann E");
			fann_destroy(ann);
			return NULL;
		}
		connected_neurons[i] = first_neuron + input_neuron;
	}

	return ann;
}


struct fann *fann_create_from_fd_1_1(FILE * conf)
{
	unsigned int num_layers, layer_size, input_neuron, i, network_type, num_connections;
	unsigned int activation_function_hidden, activation_function_output;
	float activation_steepness_hidden, activation_steepness_output;
	float learning_rate, connection_rate;
	struct fann_neuron *first_neuron, *neuron_it, *last_neuron, **connected_neurons;
	float *weights;
	struct fann_layer *layer_it;
	struct fann *ann;

	if(fscanf(conf, "%u %f %f %u %u %u " "%f"" " "%f" "\n", &num_layers, &learning_rate,
		&connection_rate, &network_type, &activation_function_hidden,
		&activation_function_output, &activation_steepness_hidden,
		&activation_steepness_output) != 8)
	{
	        printf("[ERROR] in fann D");
		return NULL;
	}

	ann = fann_allocate_structure(num_layers);
	if(ann == NULL)
	{
		return NULL;
	}
	ann->connection_rate = connection_rate;
	ann->network_type = (enum fann_nettype_enum)network_type;
	ann->learning_rate = learning_rate;

	/* determine how many neurons there should be in each layer */
	for(layer_it = ann->first_layer; layer_it != ann->last_layer; layer_it++)
	{
		if(fscanf(conf, "%u ", &layer_size) != 1)
		{
	                printf("[ERROR] in fann A");
			fann_destroy(ann);
			return NULL;
		}
		layer_it->first_neuron = NULL;
		layer_it->last_neuron = layer_it->first_neuron + layer_size;
		ann->total_neurons += layer_size;
	}

	ann->num_input = ann->first_layer->last_neuron - ann->first_layer->first_neuron - 1;
	ann->num_output = ((ann->last_layer - 1)->last_neuron - (ann->last_layer - 1)->first_neuron);
	if(ann->network_type == FANN_NETTYPE_LAYER)
	{
		ann->num_output--;
	}

	fann_allocate_neurons(ann);
	last_neuron = (ann->last_layer - 1)->last_neuron;
	for(neuron_it = ann->first_layer->first_neuron; neuron_it != last_neuron; neuron_it++)
	{
		if(fscanf(conf, "%u ", &num_connections) != 1)
		{
	                printf("[ERROR] in fann C");
			fann_destroy(ann);
			return NULL;
		}
		neuron_it->first_con = ann->total_connections;
		ann->total_connections += num_connections;
		neuron_it->last_con = ann->total_connections;
	}

	fann_allocate_connections(ann);
	connected_neurons = ann->connections;
	weights = ann->weights;
	first_neuron = ann->first_layer->first_neuron;

	for(i = 0; i < ann->total_connections; i++)
	{
		if(fscanf(conf, "(%u " "%f" ") ", &input_neuron, &weights[i]) != 2)
		{
	                printf("[ERROR] in fann B");
			fann_destroy(ann);
			return NULL;
		}
		connected_neurons[i] = first_neuron + input_neuron;
	}

	fann_set_activation_steepness_hidden(ann, activation_steepness_hidden);
	fann_set_activation_steepness_output(ann, activation_steepness_output);
	fann_set_activation_function_hidden(ann, (enum fann_activationfunc_enum)activation_function_hidden);
	fann_set_activation_function_output(ann, (enum fann_activationfunc_enum)activation_function_output);

	return ann;
}

void fann_allocate_neurons(struct fann *ann)
{
	struct fann_layer *layer_it;
	struct fann_neuron *neurons;
	unsigned int num_neurons_so_far = 0;
	unsigned int num_neurons = 0;

	/* all the neurons is allocated in one long array (calloc clears mem) */
	neurons = (struct fann_neuron *) calloc(ann->total_neurons, sizeof(struct fann_neuron));
	ann->total_neurons_allocated = ann->total_neurons;

	if(neurons == NULL)
	{
	        printf("[ERROR] in fann_allocate_neurons");
		return;
	}

	for(layer_it = ann->first_layer; layer_it != ann->last_layer; layer_it++)
	{
		num_neurons = layer_it->last_neuron - layer_it->first_neuron;
		layer_it->first_neuron = neurons + num_neurons_so_far;
		layer_it->last_neuron = layer_it->first_neuron + num_neurons;
		num_neurons_so_far += num_neurons;
	}

	ann->output = (float *) calloc(num_neurons, sizeof(float));
	if(ann->output == NULL)
	{
	        printf("[ERROR] in fann_allocate_neurons");
		return;
	}
}

void fann_allocate_connections(struct fann *ann)
{
	ann->weights = (float *) calloc(ann->total_connections, sizeof(float));
	if(ann->weights == NULL)
	{
	        printf("[ERROR] in fann_allocate_connections");
		return;
	}
	ann->total_connections_allocated = ann->total_connections;
	ann->connections = (struct fann_neuron **) calloc(ann->total_connections_allocated,sizeof(struct fann_neuron *));
	if(ann->connections == NULL)
	{
	        printf("[ERROR] in fann_allocate_connections");
		return;
	}
}

float **initTabFloatDim2( unsigned int Dim1 , unsigned int Dim2)
{
    unsigned int i;
    float **tab = (float **) calloc ( Dim1 , sizeof(float *));
    for (i=0 ; i<Dim1 ; i++)
    {
  tab[i] = (float *) calloc (Dim2 , sizeof (float) );
    }
    return tab;
}

void freeTabFloatDim2 ( float ** tab ,unsigned int dim1 )
{
    unsigned int i=0;
    for ( i=0 ; i < dim1 ; i++)
    {
      free(tab[ i ]);
    }
    free(tab);
}

int fann_allocate_scale(struct fann *ann)
{
	/* todo this should only be allocated when needed */
	unsigned int i = 0;
#define SCALE_ALLOCATE( what, where, default_value )		    			\
		ann->what##_##where = (float *)calloc(								\
			ann->num_##where##put,											\
			sizeof( float )													\
			);																\
		if( ann->what##_##where == NULL )									\
		{																	\
		/*	fann_error( NULL, FANN_E_CANT_ALLOCATE_MEM );				*/	\
			fann_destroy( ann );                            				\
			return 1;														\
		}																	\
		for( i = 0; i < ann->num_##where##put; i++ )						\
			ann->what##_##where[ i ] = ( default_value );

	SCALE_ALLOCATE( scale_mean,		in,		0.0 )
	SCALE_ALLOCATE( scale_deviation,	in,		1.0 )
	SCALE_ALLOCATE( scale_new_min,	in,		-1.0 )
	SCALE_ALLOCATE( scale_factor,		in,		1.0 )

	SCALE_ALLOCATE( scale_mean,		out,	0.0 )
	SCALE_ALLOCATE( scale_deviation,	out,	1.0 )
	SCALE_ALLOCATE( scale_new_min,	out,	-1.0 )
	SCALE_ALLOCATE( scale_factor,		out,	1.0 )
#undef SCALE_ALLOCATE
	return 0;
}
