// defs.h - Type definitions and structures
// 
// Author: Tuomas Aaltonen : tuomas.aaltonen@tuni.fi
//

#ifndef DEFS_H
#define DEFS_H

#include "ac_fixed.h"
//#include "ovf_ac_fixed.h"
//#include "ovf_ac_fixed_fns.h"
#include "ac_math.h"
#include "ac_float.h"
#include "ac_complex.h"
#include "ac_std_float.h"
#include "ac_channel.h"

#include <vector>
#include <complex>
#include <iomanip>
#include <cmath>
#include <fstream>

//precisions for std::fixed
const int PRECI = 15;

//Number of sample data instances
//const int SAMPLE_NUM = 910;
const int SAMPLE_NUM = 10;

//RxData Input sample dimensions
const int RX_ROWS = 4;
const int RX_COLS = 256;

//Slice of input data fed into DUT
const int RX_ROWS_SLICE = 4;
const int RX_COLS_SLICE = 1; //one column at a time

//Jacobi method
const int RXX_SIZE = 4;
const int MAX_ITER = 16;

const int IJM_SEARCH_NUM = ((RXX_SIZE-1)*RXX_SIZE)/2; //Number of off-diagonal elements
const int SEARCH_IDX = ac::nbits<IJM_SEARCH_NUM>::val;
typedef ac_int<SEARCH_IDX,false> ijm_search_idx_t;

const double ERROR_TOLERANCE = 0.005; //Jacobi termination

//Array manifold dimensions
const int MANIFOLD_TP_ROWS = 181;
const int MANIFOLD_TP_COLS = 4;

const int MANIFOLD_ROWS = 4;
const int MANIFOLD_COLS = 181;

//Noise subspace
const int UN_ROWS = 4;
const int UN_COLS = 3;

const int RXX_IDX = ac::nbits<RXX_SIZE-1>::val;
typedef ac_int<RXX_IDX,false> ijm_rxx_idx_t;

//Peak search
const int MAX_ANGLE = 90;
const int MIN_ANGLE = -90;

const int ANGLE_DIFF = ac::nbits<MAX_ANGLE+(-MIN_ANGLE)>::val;
typedef ac_int<ANGLE_DIFF,true> sms_out_t;

const int ANGLE_INDEX = ac::nbits<MANIFOLD_TP_ROWS>::val;
typedef ac_int<ANGLE_INDEX,false> sms_index_t;

//Debug flag
//#define MUSIC_DEBUG


//Type definitions
#define FIXED
#ifndef FIXED

//SCM
typedef ac_ieee_float64 scm_in_t;
typedef ac_ieee_float64 scm_out_t;
typedef ac_ieee_float64 scm_res_t;
typedef ac_ieee_float64 scm_acc_t;

//IJM
typedef ac_ieee_float64 ijm_out_val_t;
typedef ac_ieee_float64 ijm_out_vec_t;

//MATRIX MULTIPLICATIONs
typedef ac_ieee_float64 ijm_vec_accu_t;
typedef ac_ieee_float64 ijm_val_accu_t;
typedef ac_ieee_float64 ijm_r_t;
typedef ac_ieee_float64 ijm_tmp_accu_t;
typedef ac_ieee_float64 ijm_tmp_t;

//DIAG SEARCH
typedef ac_ieee_float64 ijm_diag_max_t;

//ROTATION MATRIX
typedef ac_ieee_float64 ijm_a_t;
typedef ac_ieee_float64 ijm_b_t;
typedef ac_ieee_float64 ijm_c_nom_t;
typedef ac_ieee_float64 ijm_s_nom_t;
typedef ac_ieee_float64 ijm_cs_den_t;
typedef ac_ieee_float64 ijm_cs_reci_t;
typedef ac_ieee_float64 ijm_c_arg_t;
typedef ac_ieee_float64 ijm_s_arg_t;

typedef ac_ieee_float64 ijm_s_t;
typedef ac_ieee_float64 ijm_c_t;

typedef ac_ieee_float64 ijm_r_diff_t;
typedef ac_ieee_float64 ijm_norm_max_t;

typedef ac_ieee_float64 ijm_denom2t_norm_t;
typedef ac_ieee_float64 ijm_denom2t_diff_t;
typedef ac_ieee_float64 ijm_denom2t_t;
typedef ac_ieee_float64 ijm_reci2t_t;
typedef ac_ieee_float64 ijm_sqrt_denom2t_t;

//SIN2T
typedef ac_ieee_float64 ijm_densin2t_t;
typedef ac_ieee_float64 ijm_sqrt_den_sin2t_t;
typedef ac_ieee_float64 ijm_nom_sin2t_t;
typedef ac_ieee_float64 ijm_sin2t_t;
typedef ac_ieee_float64 ijm_sin2tsqrt_reci_t;

//COS2T
typedef ac_ieee_float64 ijm_nom_cos2t_t;
typedef ac_ieee_float64 ijm_cos2t_t;

//COST
typedef ac_ieee_float64 ijm_nom_cost_t ;
typedef ac_ieee_float64 ijm_arg_cost_t;
typedef ac_ieee_float64 ijm_cost_t;

//SINT
typedef ac_ieee_float64 ijm_nom_sint_t;
typedef ac_ieee_float64 ijm_den_sint_t;
typedef ac_ieee_float64 ijm_sint_t;
typedef ac_ieee_float64 ijm_den_sint_reci_t;

const ac_ieee_float64 CONSTANT_FOUR = 4.0;
const ac_ieee_float64 CONSTANT_TWO = 2.0;
const ac_ieee_float64 CONSTANT_ONE = 1.0;

//NSS
//typedef ac_ieee_float64 nss_out_t;
typedef ac_ieee_float64 nss_max_val_t;

//SMS
typedef ac_ieee_float64 sms_manifold_accu_t;
typedef ac_ieee_float64 sms_res_t;

typedef ac_ieee_float64 sms_sum_t;
typedef ac_ieee_float64 sms_abs_t;
typedef ac_ieee_float64 sms_tmp_t;
typedef ac_ieee_float64 sms_min_t;
typedef ac_ieee_float64 sms_manifold_t;

//typedef int sms_out_t;

//Fixed-point type definitions
#else

//SCM
typedef ac_fixed<16,1,true> scm_in_t;
typedef ac_fixed<17,7,true> scm_out_t;
typedef ac_fixed<16,1,true> scm_res_t;
typedef ac_fixed<16,1,true> scm_acc_t;

//IJM
typedef ac_fixed<18,8,true> ijm_out_val_t;
typedef ac_fixed<12,2,true> ijm_out_vec_t;

//MATRIX MULTIPLICATIONs
typedef ac_fixed<12,2,true> ijm_vec_accu_t;
typedef ac_fixed<18,8,true> ijm_val_accu_t;
typedef ac_fixed<12,2,true> ijm_r_t;
typedef ac_fixed<18,8,true> ijm_tmp_accu_t;
typedef ac_fixed<18,8,true> ijm_tmp_t;

//DIAG SEARCH
typedef ac_fixed<19,9,false> ijm_diag_max_t;

//ROTATION MATRIX
typedef ac_fixed<16,6,true> ijm_a_t;
typedef ac_fixed<16,6,true> ijm_b_t;

typedef ac_fixed<15,5,false> ijm_cs_sqrt_t;
typedef ac_fixed<14,4,true> ijm_cs_reci_t;

typedef ac_fixed<12,2,true> ijm_c_arg_t;
typedef ac_fixed<12,2,true> ijm_s_arg_t;

typedef ac_fixed<12,2,true> ijm_s_t;
typedef ac_fixed<12,2,true> ijm_c_t;

typedef ac_fixed<17,7,true> ijm_r_diff_t;
typedef ac_fixed<15,5,false> ijm_norm_max_t;

typedef ac_fixed<22,12,false> ijm_denom2t_norm_t;
typedef ac_fixed<22,12,false> ijm_denom2t_diff_t;
typedef ac_fixed<22,12,false> ijm_denom2t_t;
typedef ac_fixed<16,6,false> ijm_sqrt_denom2t_t;
typedef ac_fixed<11,1,false> ijm_reci2t_t;

//SIN2T
typedef ac_fixed<16,6,false> ijm_nom_sin2t_t;
typedef ac_fixed<12,2,true> ijm_sin2t_t;

//COS2T
typedef ac_fixed<16,6,false> ijm_nom_cos2t_t;
typedef ac_fixed<11,1,false> ijm_cos2t_t;

//COST
typedef ac_fixed<12,2,false> ijm_nom_cost_t;
typedef ac_fixed<11,1,false> ijm_arg_cost_t;
typedef ac_fixed<11,1,false> ijm_cost_t;

//SINT
typedef ac_fixed<16,6,true> ijm_nom_sint_t;
typedef ac_fixed<12,2,false> ijm_den_sint_t;
typedef ac_fixed<13,3,false> ijm_den_sint_reci_t;
typedef ac_fixed<12,2,true> ijm_sint_t;

const int CONSTANT_FOUR = 4;
const int CONSTANT_TWO = 2;
const int CONSTANT_ONE = 1;

//NSS
//typedef ac_fixed<64,32,true> nss_out_t;
typedef ac_fixed<18,8,true> nss_max_val_t;

//SMS
typedef ac_fixed<22,2,true> sms_manifold_accu_t;
typedef ac_fixed<22,2,true> sms_res_t;

//typedef ac_fixed<24,3,true> sms_sum_t;
typedef ac_fixed<23,3,false> sms_sum_t;
typedef ac_fixed<24,4,false> sms_min_t;
typedef ac_fixed<18,2,true> sms_manifold_t;


#endif

//Type independent complex values
//SCM
typedef ac_complex<scm_in_t> scm_in_cpx_t;
typedef ac_complex<scm_out_t> scm_out_cpx_t;
typedef ac_complex<scm_res_t> scm_res_cpx_t;
typedef ac_complex<scm_acc_t> scm_acc_cpx_t;

//IJM
typedef ac_complex<ijm_out_val_t> ijm_out_val_cpx_t;
typedef ac_complex<ijm_out_vec_t> ijm_out_vec_cpx_t;

typedef ac_complex<ijm_vec_accu_t> ijm_vec_accu_cpx_t;
typedef ac_complex<ijm_val_accu_t> ijm_val_accu_cpx_t;
typedef ac_complex<ijm_r_t> ijm_r_cpx_t;
typedef ac_complex<ijm_tmp_accu_t> ijm_tmp_accu_cpx_t;
typedef ac_complex<ijm_tmp_t> ijm_tmp_cpx_t;

typedef ac_complex<ijm_s_t> ijm_s_cpx_t;
typedef ac_complex<ijm_c_t> ijm_c_cpx_t;

//NSS
//typedef ac_complex<nss_out_t> nss_out_cpx_t;

//SMS
typedef ac_complex<sms_sum_t> sms_sum_cpx_t;
typedef ac_complex<sms_manifold_t> manifold_cpx_t;
typedef ac_complex<sms_manifold_accu_t> sms_manifold_accu_cpx_t;
typedef ac_complex<sms_res_t> sms_res_cpx_t;

//Ac_channel interconnect data structs
struct inStruct_t{
    scm_in_cpx_t data[RX_ROWS_SLICE][RX_COLS_SLICE];
    bool rdy_flag;
};

struct rxxStruct_t{
    scm_out_cpx_t data[RXX_SIZE][RXX_SIZE];
};

struct nssStruct_t{
    ijm_out_vec_cpx_t data[UN_ROWS][UN_COLS];
};


struct eigenStruct_t{
    ijm_out_vec_cpx_t vec_data[RXX_SIZE][RXX_SIZE];
    ijm_out_val_t val_data[RXX_SIZE];
};


//Reference data and file read datatypes
typedef std::vector<std::vector<std::complex<double>>> matrix;
typedef std::vector<std::complex<double>> matrix_row;

#endif // DEFS_H
