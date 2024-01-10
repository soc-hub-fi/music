// tb_music_doa.cpp - Testbench for MUSIC
// 
// Author: Tuomas Aaltonen : tuomas.aaltonen@tuni.fi
//

#include "defs.h"
#include "music_top_class.h"

#include <mc_scverify.h>

//Test data
matrix loadMatrixFile(const char* filename, int rows, int cols,double scale);

CCS_MAIN(int argc, char *argv[]){

    //Interconnect channels
    static ac_channel<inStruct_t> inChan;
    static ac_channel<sms_out_t> outChan;

    //Interconnect data
    inStruct_t inData;
    sms_out_t outData;

    //Debug channels
    static ac_channel<rxxStruct_t> SCMDbgChan;
    static ac_channel<nssStruct_t> NSSDbgChan;
    static ac_channel<eigenStruct_t> IJMDbgChan;

    //Debug data
    rxxStruct_t SCMDbgData;
    eigenStruct_t IJMDbgData;
    nssStruct_t NSSDbgData;

    //Dut, golden reference instances
    music_top_class DUT;

    //Input data
    std::string name_s = "data256/testi";
    std::string name_e = ".txt";

    //File for output data
    /*
    FILE *doaFile;
    doaFile = fopen("doa_outdata.txt","w");
    if(doaFile == NULL)
    {
        std::cout << "Could not open doa output file" << std::endl;
    }
    */

    //Iterate through sample datasets
    for(int i=1; i<SAMPLE_NUM+1; i++)
    {

    //Read signal data from file
    std::string tmp_s = name_s+std::to_string(i)+name_e;
    const char *c = tmp_s.c_str();
    matrix data = loadMatrixFile(c,RX_ROWS,RX_COLS,0.1);

    //Feed data into DUT
    scm_in_cpx_t val;
    for(int j=0; j<RX_COLS; j+=RX_COLS_SLICE){

        for(int k=0; k<RX_COLS_SLICE; k++){

            for(int i=0; i<RX_ROWS; i++){
                val.set_r(data[i][j+k].real());
                val.set_i(data[i][j+k].imag());
                inData.data[i][k] = val;
            }
        }

        if(j==RX_COLS-RX_COLS_SLICE){
            //rd_flag = 1;
            inData.rdy_flag = true;
        }else{
            inData.rdy_flag = false;
        }

        inChan.write(inData);

//Ouput/Debug data
#ifdef MUSIC_DEBUG
        //GET DOA ESTIMATE
        DUT.run(inChan,SCMDbgChan,IJMDbgChan,NSSDbgChan,outChan);

#else
        DUT.run(inChan,outChan);

#endif //MUSIC_DEBUG


        if(SCMDbgChan.available(1)){
            SCMDbgData = SCMDbgChan.read();

            std::cout << "DUT // (GOLDEN)" << std::endl;

            std::cout << "RXX: " << std::endl;
            for(int i=0; i<RXX_SIZE; i++){
                for(int j=0; j<RXX_SIZE; j++){
                    std::cout << std::fixed << std::setprecision(PRECI) << "Index [" << i << j << "]: " <<
                    SCMDbgData.data[i][j].r().to_double() << " " << SCMDbgData.data[i][j].i().to_double() << std::endl;
                }
            }
        }

        if(IJMDbgChan.available(1)){

            IJMDbgData = IJMDbgChan.read();

            std::cout << "EIGVAL" << std::endl;
            for(int i=0; i<RXX_SIZE; i++){
                std::cout << std::fixed << std::setprecision(PRECI) << "Index [" << i << "]: " <<
                IJMDbgData.val_data[i].to_double()  << std::endl;
            }

            std::cout << "EIGVEC" << std::endl;
            for(int i=0; i<RXX_SIZE; i++){
                for(int j=0; j<RXX_SIZE; j++){
                    std::cout << std::fixed << std::setprecision(PRECI) << "Index [" << i << j << "]: " <<
                    IJMDbgData.vec_data[i][j].r().to_double() << " " << IJMDbgData.vec_data[i][j].i().to_double() << std::endl;
                }
            }
        }

        if(NSSDbgChan.available(1)){

            NSSDbgData = NSSDbgChan.read();

            std::cout << "UN" << std::endl;
            for(int i=0; i<UN_ROWS; i++){
                for(int j=0; j<UN_COLS; j++){
                    std::cout << std::fixed << std::setprecision(PRECI) << "Index [" << i << j << "]: " <<
                    NSSDbgData.data[i][j].r().to_double() << std::endl;
                }
            }

        }

        if(outChan.available(1)){

            outData = outChan.read();

            std::cout << "DOA: " << outData << std::endl;

            //Write output into file
            //fprintf(doaFile,"%d\n",outData);

        }//End if

    } //End COLS

    }//End SAMPLE_NUM

    //fclose(doaFile); //DoA estimates
    CCS_RETURN(0);

}//END main

//Data read from file, scale to appropriate data range
matrix loadMatrixFile(const char* filename, int rows, int cols,double scale){

    matrix mat;

    std::ifstream inFile;
    std::string line;
    inFile.open(filename);

    if(!inFile){
        std::cout << "File " << filename << "did not open" << std::endl;
    }

    double rl;
    double cpx;

    std::vector<std::complex<double>> row;
    std::complex<double> val;

    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++){

            std::getline(inFile,line);
            std::stringstream ss(line);
            ss >> rl >> cpx;

            val.real(rl*scale);
            val.imag(cpx*scale);

            row.push_back(val);
        }

        mat.push_back(row);
        row.erase(row.begin(),row.end());
    }

    return mat;
}
