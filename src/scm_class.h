// scm_class.h - Sensor Covariance Matrix (SCM) hierarchical block
// 
// Author: Tuomas Aaltonen : tuomas.aaltonen@tuni.fi
//

#ifndef SCM_CLASS_H
#define SCM_CLASS_H

#include "defs.h"

#include "Vector_Matrix_Ops_class.h"

#pragma design
class SCM_class
{

private:

    //Store accumulated covariance matrix
    rxxStruct_t outStruct;

    //Matrix multiply instance
    matrix_x_matrix_multiply_class<scm_in_cpx_t,scm_in_cpx_t,scm_acc_cpx_t,
                                    scm_res_cpx_t,RX_ROWS_SLICE,RX_COLS_SLICE,
                                    RX_ROWS_SLICE,RX_COLS_SLICE,RXX_SIZE,RXX_SIZE,
                                    false,true> Mult_inst;

public:

    SCM_class()
    {
        //Init covariance matrix to zero
        ac::init_array<AC_VAL_0>(&outStruct.data[0][0],RXX_SIZE*RXX_SIZE);
    };

#ifdef MUSIC_DEBUG
    #pragma hls_design interface
    void    CCS_BLOCK(run)(ac_channel<inStruct_t> &dataIn_ch,
                ac_channel<rxxStruct_t> &dataOut_ch,
                ac_channel<rxxStruct_t> &debugOut_ch)
#else
    #pragma hls_design interface
    void    CCS_BLOCK(run)(ac_channel<inStruct_t> &dataIn_ch,
                ac_channel<rxxStruct_t> &dataOut_ch)
#endif
    {
        //Read data if available
        if(dataIn_ch.available(1)){

            //Channel data read
            inStruct_t inData = dataIn_ch.read();

            //Matrix multiplication result
            scm_res_cpx_t RXX_MTX[RXX_SIZE][RXX_SIZE];

            //Matrix multiplication
            Mult_inst.Product(inData.data,inData.data,RXX_MTX);

            //Store accumulated result
            #pragma hls_unroll yes
            ACCUM_ROWS: for(int i=0; i<RXX_SIZE; i++){
                #pragma hls_unroll yes
                ACCUM_COLS: for(int j=0; j<RXX_SIZE; j++){
                    outStruct.data[i][j] += RXX_MTX[i][j];
                }
            }

            //If end of data instance send data further
            if(inData.rdy_flag == 1)
            {
                dataOut_ch.write(outStruct);
#ifdef MUSIC_DEBUG
                debugOut_ch.write(outStruct);
#endif //MUSIC_DEBUG

                #pragma hls_unroll yes
                ZERO_ROWS: for(int i=0; i<RXX_SIZE; i++){
                    #pragma hls_unroll yes
                    ZERO_COLS: for(int j=0; j<RXX_SIZE; j++){
                        outStruct.data[i][j] = 0;
                    }
                }

            }//rd_flag

        }//dataIn_available()

    }//run
};

#endif // SCM_CLASS_H
