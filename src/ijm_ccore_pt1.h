// ijm_ccore_pt1.h - CCORE for off-diagonal max. search
//
// Author: Tuomas Aaltonen : tuomas.aaltonen@tuni.fi
//

#ifndef IJM_CCORE_PT1_H
#define IJM_CCORE_PT1_H

#include "max_search.h"
#include "defs.h"
#include <mc_scverify.h>

class ijm_ccore_pt1
{
    private:


        //Off-diagonal search values
        void searchArrayCreate(ijm_out_val_cpx_t EIGVAL_MTX[RXX_SIZE][RXX_SIZE],
                            element_struct SEARCH_ARRAY[IJM_SEARCH_NUM])
        {
            ijm_search_idx_t k=0;
            #pragma hls_unroll yes
            SEARCH_ROWS: for(int i=0; i<RXX_SIZE-1; i++){
                #pragma hls_unroll yes
                SEARCH_COLS: for(int j=RXX_SIZE-1; j>i; j--){

                    element_struct tmp_s;

                    ijm_a_t tmp_r = EIGVAL_MTX[i][j].real();
                    ijm_b_t tmp_i = EIGVAL_MTX[i][j].imag();

                    ijm_diag_max_t tmp_max = 0;
                    tmp_max = tmp_r*tmp_r+tmp_i*tmp_i;

                    tmp_s.val = tmp_r*tmp_r+tmp_i*tmp_i;
                    tmp_s.r = i;
                    tmp_s.c = j;

                    SEARCH_ARRAY[k] = tmp_s;
                    k++;
                }
            }
        }

        //Off-diagonal max. search
        //Value and indices
        void offDiagMax(element_struct SEARCH_ARRAY[IJM_SEARCH_NUM],
                        ijm_diag_max_t &max_val,ijm_rxx_idx_t &r_idx,
                        ijm_rxx_idx_t &c_idx)
        {
            element_struct max_struct = max<IJM_SEARCH_NUM>(SEARCH_ARRAY);

            max_val = max_struct.val;
            r_idx = max_struct.r;
            c_idx = max_struct.c;
        }

    public:
        ijm_ccore_pt1() {};

        #pragma hls_design interface ccore
        void CCS_BLOCK(run)(ijm_out_val_cpx_t EIGVAL_MTX[RXX_SIZE][RXX_SIZE],ijm_rxx_idx_t &row_idx,
                            ijm_rxx_idx_t &col_idx, ijm_diag_max_t &max_val, bool &is_converged)
        {
            element_struct SEARCH_ARRAY[IJM_SEARCH_NUM];

            row_idx = 0;
            col_idx = 1;
            max_val = 0.0;

            //Create search array
            searchArrayCreate(EIGVAL_MTX,SEARCH_ARRAY);
            //Find off-diagonal max and indices
            offDiagMax(SEARCH_ARRAY,max_val,row_idx,col_idx);

            //Check convergence
            if(max_val < ERROR_TOLERANCE)
                is_converged = true;

        }
};

#endif // IJM_CCORE_PT1_H
