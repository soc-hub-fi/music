// ijm_class.h - Iterative Jacobi Method (IJM) hierarchical block
//
// Author: Tuomas Aaltonen : tuomas.aaltonen@tuni.fi
// Tampere University - Unit of Computing Sciences - 2024
//

#ifndef IJM_CLASS_H
#define IJM_CLASS_H

#include "ijm_ccore_pt1.h"
#include "ijm_ccore_pt2.h"
#include "Music_Matrix_Class.h"

#include "max_search.h"
#include "defs.h"
#include <mc_scverify.h>

#pragma design
class IJM_class
{
    private:

        //CCORE instances
        ijm_ccore_pt1 ijm_ccore_inst1;
        ijm_ccore_pt2 ijm_ccore_inst2;

        //Matrix multiplication instances
        unitary_left_mult_class<ijm_r_cpx_t,ijm_out_val_cpx_t,ijm_tmp_accu_cpx_t,ijm_tmp_cpx_t,
                                RXX_SIZE> Left_mult_val_inst;

        unitary_right_mult_class<ijm_tmp_cpx_t,ijm_r_cpx_t,ijm_val_accu_cpx_t,ijm_out_val_cpx_t, 
                                RXX_SIZE> Right_mult_val_inst;

        unitary_right_mult_class<ijm_out_vec_cpx_t,ijm_r_cpx_t,ijm_vec_accu_cpx_t,ijm_out_vec_cpx_t,
                                RXX_SIZE> Right_mult_vec_inst;

    public:
        IJM_class() {};

#ifdef MUSIC_DEBUG
    #pragma hls_design interface
    void CCS_BLOCK(run)(ac_channel<rxxStruct_t> &dataIn_ch,
                        ac_channel<eigenStruct_t> &debugOut_ch,
                        ac_channel<eigenStruct_t> &dataOut_ch)
#else
    #pragma hls_design interface
    void CCS_BLOCK(run)(ac_channel<rxxStruct_t> &dataIn_ch,
                        ac_channel<eigenStruct_t> &dataOut_ch)
#endif //MUSIC_DEBUG
    {

        if(dataIn_ch.available(1))
        {
            //Channel read
            rxxStruct_t inData = dataIn_ch.read();

            ijm_out_val_cpx_t EIGVAL_MTX[RXX_SIZE][RXX_SIZE];
            ijm_out_vec_cpx_t EIGVEC_MTX[RXX_SIZE][RXX_SIZE];

            //Init eigenvalue matrix to covariance matrix
            #pragma hls_unroll yes
            UNPACK_ROWS: for(int i=0; i<RXX_SIZE; i++){
                #pragma hls_unroll yes
                UNPACK_COLS: for(int j=0; j<RXX_SIZE; j++){
                    EIGVAL_MTX[i][j] = inData.data[i][j];
                }
            }

            //Init eigenvector matrix to identity
            #pragma hls_unroll yes
            INIT_VEC_ROWS: for(int i=0; i<RXX_SIZE; i++){
                #pragma hls_unroll yes
                INIT_VEC_COLS: for(int j=0; j<RXX_SIZE; j++){
                    EIGVEC_MTX[i][j] = 0;

                    if(i==j)
                        EIGVEC_MTX[i][j].set_r(1.0);

                }
            }

            //Iterative jacobi method
            JACOBI_LOOP: for(int iter=0; iter<MAX_ITER; iter++){

                ijm_rxx_idx_t row_idx;
                ijm_rxx_idx_t col_idx;
                ijm_diag_max_t max_val;

                bool is_converged = false;

                //Off-diagonal max. search, convergence
                ijm_ccore_inst1.run(EIGVAL_MTX, row_idx, col_idx,
                                   max_val, is_converged);

                //Convergence condition
                if(is_converged == true && iter != 0){
                    //std::cout << "DUT    Converged at iteration: " << iter << std::endl;
                    break;
                }

                //Rotation matrix
                ijm_r_cpx_t R_MTX[RXX_SIZE][RXX_SIZE];

                //Rotation matrix creation
                ijm_ccore_inst2.run(EIGVAL_MTX,R_MTX,row_idx,col_idx,max_val);

                //Intermediate matrix mult. result
                ijm_tmp_cpx_t TMP_MTX[RXX_SIZE][RXX_SIZE];

                //APPLY ROTATION TRANSFORMATIONS
                //EIGVAL = R x EIGVAL x R'
                Left_mult_val_inst.left_mult(row_idx,col_idx,R_MTX,EIGVAL_MTX,TMP_MTX);
                Right_mult_val_inst.right_mult(row_idx,col_idx,TMP_MTX,R_MTX,EIGVAL_MTX);

                //Intermediate matrix mult. result
                ijm_out_vec_cpx_t EIGVEC_TMP[RXX_SIZE][RXX_SIZE];
                //EIGVEC = EIGVEC x R'
                Right_mult_vec_inst.right_mult(row_idx,col_idx,EIGVEC_MTX,R_MTX,EIGVEC_TMP);

                //Assign values from temporary matrix
                #pragma hls_unroll yes
                UPDATE_ROWS: for(int i=0; i<RXX_SIZE; i++){
                    #pragma hls_unroll yes
                    UPDATE_COLS: for(int j=0; j<RXX_SIZE; j++){
                        EIGVEC_MTX[i][j] = EIGVEC_TMP[i][j];
                    }
                }
            }

            //Pack outputs into struct, send data further
            eigenStruct_t outStruct;

            #pragma hls_unroll yes
            PACK_VEC_ROWS: for(int i=0; i<RXX_SIZE; i++){
                #pragma hls_unroll yes
                PACK_VEC_COLS: for(int j=0; j<RXX_SIZE; j++){
                    outStruct.vec_data[i][j] = EIGVEC_MTX[i][j];
                }
            }

            #pragma hls_unroll yes
            PACK_EIGVAL: for(int i=0; i<RXX_SIZE; i++){
                outStruct.val_data[i] = EIGVAL_MTX[i][i].real();
            }

            dataOut_ch.write(outStruct);
#ifdef MUSIC_DEBUG
            debugOut_ch.write(outStruct);

#endif
        }//Data available

    }//void run

};//ijm class

#endif // IJM_CLASS_H
