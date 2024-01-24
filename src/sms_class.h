//sms_class.h - Spectral Maximum Search (SMS) hierarchical block
//
// Author: Tuomas Aaltonen : tuomas.aaltonen@tuni.fi
// Tampere University - Unit of Computing Sciences - 2024
//

#ifndef SMS_CLASS_H
#define SMS_CLASS_H


#include "defs.h"
#include "Vector_Matrix_Ops_class.h"

#include <mc_scverify.h>

class SMS_class
{
    private:

        //Array manifold matrix
        manifold_cpx_t MANI_MTX[MANIFOLD_TP_ROWS][MANIFOLD_TP_COLS];

        //Matrix multiplication instance
        matrix_x_matrix_multiply_class<manifold_cpx_t,ijm_out_vec_cpx_t,sms_manifold_accu_cpx_t,
                                    sms_res_cpx_t,MANIFOLD_TP_ROWS,MANIFOLD_TP_COLS,UN_ROWS,UN_COLS,MANIFOLD_TP_ROWS,UN_COLS,
                                    false,false> Mult_res_inst;

	//
	#include "matrix_mani.h"

    public:

        SMS_class()
        {

            //Scan manifold into local array
            for(int i=0; i<MANIFOLD_TP_ROWS; i++){
                for(int j=0; j<MANIFOLD_TP_COLS; j++){
                    MANI_MTX[i][j] = MTX_MANIFOLD[i][j];
                }
            }
        };

        #pragma hls_design interface
        void CCS_BLOCK(run)(ac_channel<nssStruct_t> &dataIn_ch,
                ac_channel<sms_out_t> &dataOut_ch)
        {

            if(dataIn_ch.available(1))
            {

                //Channel read
                nssStruct_t inData = dataIn_ch.read();

                //Matrix multiplication
                sms_res_cpx_t RES_MTX[MANIFOLD_TP_ROWS][UN_COLS];
                Mult_res_inst.Product(MANI_MTX,inData.data,RES_MTX);

                //First result of min. search
                sms_sum_t col_sum = 0;
                sms_min_t min_val = 0;
                #pragma hls_unroll yes
                FIRST_SUM: for(int i=0; i<UN_COLS; i++){

                    sms_res_cpx_t tmp_val = RES_MTX[0][i];
                    min_val = tmp_val.real()*tmp_val.real()+tmp_val.imag()*tmp_val.imag();
                    col_sum += min_val;
                }

                sms_min_t min_sum = col_sum;
                sms_index_t index = 0;

                //Find smallest value index of remaining data
                ROW_SCAN: for(int i=1; i<MANIFOLD_TP_ROWS; i++){

                    col_sum = 0;

                    //Columnwise summation
                    #pragma hls_unroll yes
                    COL_SUM:for(int j=0; j<UN_COLS; j++){

                        sms_res_cpx_t tmp_val = RES_MTX[i][j];
                        min_val = tmp_val.real()*tmp_val.real()+tmp_val.imag()*tmp_val.imag();
                        col_sum += min_val;
                    }
                    
                    //Set new minimum
                    if(col_sum < min_sum){
                        min_sum = col_sum;
                        index = i;
                    }

                }

                //Adjust index to angle search range
                sms_out_t angle = index - 90;

                //Output write
                dataOut_ch.write(angle);

            }//if dataIn_ch available

        }//run
};

#endif // SMS_CLASS_H
