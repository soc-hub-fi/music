// nss_class.h - Noise Subspace Selection (NSS) hierarchical block
// 
// Author: Tuomas Aaltonen : tuomas.aaltonen@tuni.fi
//

#ifndef NSS_CLASS_H
#define NSS_CLASS_H

#include "defs.h"

#include <mc_scverify.h>

#pragma design
class NSS_class
{
    private:

    public:

        NSS_class() {};

#ifdef MUSIC_DEBUG
        #pragma hls_design interface
        void CCS_BLOCK(run)(ac_channel<eigenStruct_t> &dataIn_ch,
                ac_channel<nssStruct_t> &debugOut_ch,
                ac_channel<nssStruct_t> &dataOut_ch)
#else
        #pragma hls_design interface
        void CCS_BLOCK(run)(ac_channel<eigenStruct_t> &dataIn_ch,
                ac_channel<nssStruct_t> &dataOut_ch)
#endif
        {

            if(dataIn_ch.available(1))
            {
                //Channel data read
                eigenStruct_t inData = dataIn_ch.read();

                //Init first index, value
                ijm_rxx_idx_t max_index = 0;
                nss_max_val_t max_value = inData.val_data[0];

                //Find index of largest eigenvalue
                #pragma hls_unroll yes
                EIGEN_SCAN: for(int i=1; i<RXX_SIZE; i++){
                    if(inData.val_data[i] > max_value){
                        max_index = i;
                        max_value = inData.val_data[i];
                    }
                }

                //Noise subspace selection
                //Omit column corresponding to largest eigenvalue
                nssStruct_t outData;
                #pragma hls_unroll yes
                SEL_ROWS: for(int i=0; i<UN_ROWS; i++){
                    #pragma hls_unroll yes
                    SEL_COLS: for(int j=0, k=0; k<UN_COLS; j++,k++){

                        if(j==max_index){
                            j++;
                        }

                        outData.data[i][k] = inData.vec_data[i][j];
                    }
                }

                //Channel write
                dataOut_ch.write(outData);

#ifdef MUSIC_DEBUG
                debugOut_ch.write(outData);
#endif

            }//data.available()

        }//run()

};

#endif //NSS_CLASS_H
