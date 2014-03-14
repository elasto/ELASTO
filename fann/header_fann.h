//USEFORTEST avgcond
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#define FANN_CONF_VERSION "FANN_FLO_2.1" 

#define fann_scanf(type, name, val) \
{ \
	if(fscanf(conf, name"="type"\n", val) != 1) \
	{ \
	/*	fann_error(NULL, FANN_E_CANT_READ_CONFIG, name, configuration_file);*/ \
		fann_destroy(ann); \
		return NULL; \
	} \
}

#define fann_safe_free(x) {if(x) { free(x); x = NULL; }}
#define fann_abs(value) (((value) > 0) ? (value) : -(value))
#define fann_mult(x,y) (x*y)
#define fann_linear_func(v1, r1, v2, r2, sum) (((((r2)-(r1)) * ((sum)-(v1)))/((v2)-(v1))) + (r1))
#define fann_stepwise(v1, v2, v3, v4, v5, v6, r1, r2, r3, r4, r5, r6, min, max, sum) (sum < v5 ? (sum < v3 ? (sum < v2 ? (sum < v1 ? min : fann_linear_func(v1, r1, v2, r2, sum)) : fann_linear_func(v2, r2, v3, r3, sum)) : (sum < v4 ? fann_linear_func(v3, r3, v4, r4, sum) : fann_linear_func(v4, r4, v5, r5, sum))) : (sum < v6 ? fann_linear_func(v5, r5, v6, r6, sum) : max))
/* FANN_LINEAR */
/*#define fann_linear(steepness, sum) fann_mult(steepness, sum) */
#define fann_linear_derive(steepness, value) (steepness)
/* FANN_SIGMOID */
//MODIF///
/*#define fann_sigmoid(steepness, sum) (1.0f/(1.0f + exp(-2.0f * steepness * sum)))  */
#define fann_sigmoid_real(sum) (1.0f/(1.0f + exp(-2.0f * sum)))
#define fann_sigmoid_derive(steepness, value) (2.0f * steepness * value * (1.0f - value))
/* FANN_SIGMOID_SYMMETRIC */
//MODIF///
/* #define fann_sigmoid_symmetric(steepness, sum) (2.0f/(1.0f + exp(-2.0f * steepness * sum)) - 1.0f) */
#define fann_sigmoid_symmetric_real(sum) (2.0f/(1.0f + exp(-2.0f * sum)) - 1.0f)
#define fann_sigmoid_symmetric_derive(steepness, value) steepness * (1.0f - (value*value))
/* FANN_GAUSSIAN */
//MODIF///
/* #define fann_gaussian(steepness, sum) (exp(-sum * steepness * sum * steepness)) */
#define fann_gaussian_real(sum) (exp(-sum * sum))
#define fann_gaussian_derive(steepness, value, sum) (-2.0f * sum * value * steepness * steepness)
/* FANN_GAUSSIAN_SYMMETRIC */
//MODIF///
/* #define fann_gaussian_symmetric(steepness, sum) ((exp(-sum * steepness * sum * steepness)*2.0)-1.0) */
#define fann_gaussian_symmetric_real(sum) ((exp(-sum * sum)*2.0f)-1.0f)
#define fann_gaussian_symmetric_derive(steepness, value, sum) (-2.0f * sum * (value+1.0f) * steepness * steepness)
/* FANN_ELLIOT */
//MODIF///
/* #define fann_elliot(steepness, sum) (((sum * steepness) / 2.0f) / (1.0f + fann_abs(sum * steepness)) + 0.5f) */
#define fann_elliot_real(sum) (((sum) / 2.0f) / (1.0f + fann_abs(sum)) + 0.5f)
#define fann_elliot_derive(steepness, value, sum) (steepness * 1.0f / (2.0f * (1.0f + fann_abs(sum)) * (1.0f + fann_abs(sum))))
/* FANN_ELLIOT_SYMMETRIC */
/*#define fann_elliot_symmetric(steepness, sum) ((sum * steepness) / (1.0f + fann_abs(sum * steepness))) */
#define fann_elliot_symmetric_real(sum) ((sum) / (1.0f + fann_abs(sum)))
#define fann_elliot_symmetric_derive(steepness, value, sum) (steepness * 1.0f / ((1.0f + fann_abs(sum)) * (1.0f + fann_abs(sum))))
/* FANN_SIN_SYMMETRIC */
/*#define fann_sin_symmetric_real(steepness,sum) (sin(steepness*sum))*/
#define fann_sin_symmetric_real(sum) (sin(sum))
#define fann_sin_symmetric_derive(steepness, sum) (steepness*cos(steepness*sum))
/* FANN_COS_SYMMETRIC */
/*#define fann_cos_symmetric_real(steepness,sum) (cos(steepness*sum))*/
#define fann_cos_symmetric_real(sum) (cos(sum))
#define fann_cos_symmetric_derive(steepness, sum) (steepness*-sin(steepness*sum))
/* FANN_SIN */
/*#define fann_sin_real(steepness,sum) (sin(steepness*sum)/2.0f+0.5f)*/
#define fann_sin_real(sum) (sin(sum)/2.0f+0.5f)
#define fann_sin_derive(steepness, sum) (steepness*cos(steepness*sum)/2.0f)
/* FANN_COS */
/*#define fann_cos_real(steepness,sum) (cos(steepness*sum)/2.0f+0.5f)*/
#define fann_cos_real(sum) (cos(sum)/2.0f+0.5f)
#define fann_cos_derive(steepness, sum) (steepness*-sin(steepness*sum)/2.0f)
#define fann_activation_switch(activation_function, value, result) \
switch(activation_function) \
{ \
	case FANN_LINEAR: \
		result = (float)value; \
        break; \
	case FANN_LINEAR_PIECE: \
		result = (float)((value < 0) ? 0 : (value > 1) ? 1 : value); \
        break; \
	case FANN_LINEAR_PIECE_SYMMETRIC: \
		result = (float)((value < -1) ? -1 : (value > 1) ? 1 : value); \
        break; \
	case FANN_SIGMOID: \
		result = (float)fann_sigmoid_real(value); \
        break; \
	case FANN_SIGMOID_SYMMETRIC: \
		result = (float)fann_sigmoid_symmetric_real(value); \
        break; \
	case FANN_SIGMOID_SYMMETRIC_STEPWISE: \
		result = (float)fann_stepwise(-2.64665293693542480469e+00, -1.47221934795379638672e+00, -5.49306154251098632812e-01, 5.49306154251098632812e-01, 1.47221934795379638672e+00, 2.64665293693542480469e+00, -9.90000009536743164062e-01, -8.99999976158142089844e-01, -5.00000000000000000000e-01, 5.00000000000000000000e-01, 8.99999976158142089844e-01, 9.90000009536743164062e-01, -1, 1, value); \
        break; \
	case FANN_SIGMOID_STEPWISE: \
		result = (float)fann_stepwise(-2.64665246009826660156e+00, -1.47221946716308593750e+00, -5.49306154251098632812e-01, 5.49306154251098632812e-01, 1.47221934795379638672e+00, 2.64665293693542480469e+00, 4.99999988824129104614e-03, 5.00000007450580596924e-02, 2.50000000000000000000e-01, 7.50000000000000000000e-01, 9.49999988079071044922e-01, 9.95000004768371582031e-01, 0, 1, value); \
        break; \
	case FANN_THRESHOLD: \
		result = (float)((value < 0) ? 0 : 1); \
        break; \
	case FANN_THRESHOLD_SYMMETRIC: \
		result = (float)((value < 0) ? -1 : 1); \
        break; \
	case FANN_GAUSSIAN: \
		result = (float)fann_gaussian_real(value); \
        break; \
	case FANN_GAUSSIAN_SYMMETRIC: \
		result = (float)fann_gaussian_symmetric_real(value); \
        break; \
	case FANN_ELLIOT: \
		result = (float)fann_elliot_real(value); \
	    break; \
	case FANN_ELLIOT_SYMMETRIC: \
		result = (float)fann_elliot_symmetric_real(value); \
        break; \
	case FANN_SIN_SYMMETRIC: \
		result = (float)fann_sin_symmetric_real(value); \
        break; \
	case FANN_COS_SYMMETRIC: \
		result = (float)fann_cos_symmetric_real(value); \
        break; \
	case FANN_SIN: \
		result = (float)fann_sin_real(value); \
        break; \
	case FANN_COS: \
		result = (float)fann_cos_real(value); \
        break; \
	case FANN_GAUSSIAN_STEPWISE: \
        result = 0; \
        break; \
}


/*-------------------------------------------
DECLARATION DES TYPES ET ENUMERATIONS 
-------------------------------------------*/
enum fann_train_enum
{
	FANN_TRAIN_INCREMENTAL = 0,
	FANN_TRAIN_BATCH,
	FANN_TRAIN_RPROP,
	FANN_TRAIN_QUICKPROP
};

static char const *const FANN_TRAIN_NAMES[] = {
	"FANN_TRAIN_INCREMENTAL",
	"FANN_TRAIN_BATCH",
	"FANN_TRAIN_RPROP",
	"FANN_TRAIN_QUICKPROP"
};

enum fann_activationfunc_enum
{
	FANN_LINEAR = 0,
	FANN_THRESHOLD,
	FANN_THRESHOLD_SYMMETRIC,
	FANN_SIGMOID,
	FANN_SIGMOID_STEPWISE,
	FANN_SIGMOID_SYMMETRIC,
	FANN_SIGMOID_SYMMETRIC_STEPWISE,
	FANN_GAUSSIAN,
	FANN_GAUSSIAN_SYMMETRIC,
	FANN_GAUSSIAN_STEPWISE,
	FANN_ELLIOT,
	FANN_ELLIOT_SYMMETRIC,
	FANN_LINEAR_PIECE,
	FANN_LINEAR_PIECE_SYMMETRIC,
	FANN_SIN_SYMMETRIC,
	FANN_COS_SYMMETRIC,
	FANN_SIN,
	FANN_COS
};

static char const *const FANN_ACTIVATIONFUNC_NAMES[] = {
	"FANN_LINEAR",
	"FANN_THRESHOLD",
	"FANN_THRESHOLD_SYMMETRIC",
	"FANN_SIGMOID",
	"FANN_SIGMOID_STEPWISE",
	"FANN_SIGMOID_SYMMETRIC",
	"FANN_SIGMOID_SYMMETRIC_STEPWISE",
	"FANN_GAUSSIAN",
	"FANN_GAUSSIAN_SYMMETRIC",
	"FANN_GAUSSIAN_STEPWISE",
	"FANN_ELLIOT",
	"FANN_ELLIOT_SYMMETRIC",
	"FANN_LINEAR_PIECE",
	"FANN_LINEAR_PIECE_SYMMETRIC",
	"FANN_SIN_SYMMETRIC",
	"FANN_COS_SYMMETRIC",
	"FANN_SIN",
	"FANN_COS"
};

enum fann_stopfunc_enum
{
	FANN_STOPFUNC_MSE = 0,
	FANN_STOPFUNC_BIT
};

static char const *const FANN_STOPFUNC_NAMES[] = {
	"FANN_STOPFUNC_MSE",
	"FANN_STOPFUNC_BIT"
};

enum fann_nettype_enum
{
    FANN_NETTYPE_LAYER = 0, /* Each layer only has connections to the next layer */
    FANN_NETTYPE_SHORTCUT /* Each layer has connections to all following layers */
};

static char const *const FANN_NETTYPE_NAMES[] = {
	"FANN_NETTYPE_LAYER",
	"FANN_NETTYPE_SHORTCUT"
};


/* forward declarations for use with the callback */
//struct fann;
//struct fann_train_data;
//FANN_EXTERNAL typedef int (FANN_API * fann_callback_type) (struct fann *ann, struct fann_train_data *train,unsigned int max_epochs,unsigned int epochs_between_reports,float desired_error, unsigned int epochs);


/* ----- Data structures -----
 * No data within these structures should be altered directly by the user.
 */

struct fann_neuron
{
	/* Index to the first and last connection
	 * (actually the last is a past end index)
	 */
	unsigned int first_con;
	unsigned int last_con;
	/* The sum of the inputs multiplied with the weights */
	float sum;
	/* The value of the activation function applied to the sum */
	float value;
	/* The steepness of the activation function */
	float activation_steepness;
	/* Used to choose which activation function to use */
	enum fann_activationfunc_enum activation_function;
};

/* A single layer in the neural network.
 */
struct fann_layer
{
	struct fann_neuron *first_neuron;

	/* A pointer to the neuron past the last neuron in the layer */
	/* the number of neurons is last_neuron - first_neuron */
	struct fann_neuron *last_neuron;
};
/*
struct fann_error
{
	enum fann_errno_enum errno_f;
	FILE *error_log;
	char *errstr;
};
*/

struct fann
{
	//enum fann_errno_enum errno_f;
	FILE *error_log;
	char *errstr;
	float learning_rate;
	float learning_momentum;
	float connection_rate;
	enum fann_nettype_enum network_type;
	struct fann_layer *first_layer;
	struct fann_layer *last_layer;
	unsigned int total_neurons;
	unsigned int num_input;
	unsigned int num_output;
	float *weights;
	struct fann_neuron **connections;
	float *train_errors;
	enum fann_train_enum training_algorithm;
	unsigned int total_connections;
	float *output;
	unsigned int num_MSE;
	float MSE_value;
	unsigned int num_bit_fail;
	float bit_fail_limit;
//	enum fann_errorfunc_enum train_error_function;
//	enum fann_stopfunc_enum train_stop_function;
//	fann_callback_type callback;
//        void *user_data;
	float cascade_output_change_fraction;
	unsigned int cascade_output_stagnation_epochs;
	float cascade_candidate_change_fraction;
	unsigned int cascade_candidate_stagnation_epochs;
	unsigned int cascade_best_candidate;
	float cascade_candidate_limit;
	float cascade_weight_multiplier;
	unsigned int cascade_max_out_epochs;
	unsigned int cascade_max_cand_epochs;	
	enum fann_activationfunc_enum *cascade_activation_functions;
	unsigned int cascade_activation_functions_count;
	float *cascade_activation_steepnesses;
	unsigned int cascade_activation_steepnesses_count;
	unsigned int cascade_num_candidate_groups;
	float *cascade_candidate_scores;
	unsigned int total_neurons_allocated;
	unsigned int total_connections_allocated;
	float quickprop_decay;
	float quickprop_mu;
	float rprop_increase_factor;
	float rprop_decrease_factor;
	float rprop_delta_min;
	float rprop_delta_max;
	float rprop_delta_zero;
	float *train_slopes;
	float *prev_steps;
	float *prev_train_slopes;
	float *prev_weights_deltas;
	float *scale_mean_in;
	float *scale_deviation_in;
	float *scale_new_min_in;
	float *scale_factor_in;
	float *scale_mean_out;
	float *scale_deviation_out;
	float *scale_new_min_out;
	float *scale_factor_out;
};

struct fann_connection
{
    long int from_neuron;
    unsigned int to_neuron;
    float weight;
};


/*-------------------------------------------
DECLARATION DES ROUTINES
-------------------------------------------*/
//void fann_init_error_data(struct fann_error *errdat);
void fann_destroy(struct fann *ann);
struct fann *fann_allocate_structure(unsigned int num_layers);
void fann_set_activation_function_hidden(struct fann *ann,enum fann_activationfunc_enum activation_function);
void fann_set_activation_function_output(struct fann *ann,enum fann_activationfunc_enum activation_function);
void fann_set_activation_steepness_output(struct fann *ann,float steepness);
void fann_set_activation_steepness_hidden(struct fann *ann,float steepness);
void fann_allocate_neurons(struct fann *ann);
struct fann *fann_create_from_file(const char *configuration_file);
float *fann_run(struct fann * ann, float * input);
struct fann *fann_create_from_fd(FILE * conf);
struct fann *fann_create_from_fd_1_1(FILE * conf);
void fann_allocate_connections(struct fann *ann);
int fann_allocate_scale(struct fann *ann);
float **initTabFloatDim2( unsigned int Dim1 , unsigned int Dim2);
void freeTabFloatDim2 ( float ** tab ,unsigned int dim1 );
