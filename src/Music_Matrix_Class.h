// Music_Matrix_Class.h - Matrix multiplication for symmetry and sparsity
//
// Author: Tuomas Aaltonen : tuomas.aaltonen@tuni.fi
//

#ifndef MUSIC_MATRIX_CLASS_H
#define MUSIC_MATRIX_CLASS_H

#include "defs.h"
#include "Vector_Matrix_Ops_class.h"

//Exploit conjugate symmetry in covariance matrix calculation
template<typename inA_T, typename inB_T, typename accu_T, typename out_T, int inA_rows, int inA_cols>
class sym_covar_mult_class
{
    
private:
    multiply_accumulate_class<inA_T,inB_T,accu_T,out_T,inA_cols,false,true> MAC;

public:

    sym_covar_mult_class() {};

void CCS_BLOCK(sym_mult)(inA_T inDataA[inA_rows][inA_cols],
                    inB_T inDataB[inA_rows][inA_cols],
                    out_T outData[inA_rows][inA_rows])
{

    inA_T Arow[inA_cols];
    inB_T Bcol[inA_rows];

    #pragma hls_unroll yes
    SCAN_ROW: for(int i=0; i<inA_rows; ++i)
    {
        #pragma hls_unroll yes
        COPY_ROW: for(int k=0; k<inA_cols; ++k)
        {
            Arow[k] = inDataA[i][k];
        }

        //Compute diagonal and upper triangle
        #pragma hls_unroll yes
        SCAN_COL: for(int j=i; j<inA_rows; ++j)
        {   
            #pragma hls_unroll yes
            COPY_COL: for(int l=0; l<inA_cols; ++l)
            {
                Bcol[l] = inDataB[j][l];
            }

            outData[i][j] = MAC.MAC(Arow,Bcol);
        }
    }

    //Mirror lower triangle over diagonal
    #pragma hls_unroll yes
    ASSIGN_ROW: for(int i=0; i<inA_rows-1; ++i)
    {
        #pragma hls_unroll yes
        ASSIGN_COL: for(int j=inA_rows-1; j>i; --j)
        {
            out_T tmp = outData[i][j];
            outData[j][i] = tmp.conj();
        }
    }


}//sym_mult

};//class


//Exploit sparsity in rotation transformations
template<typename inA_T, typename inB_T, typename accu_T, typename out_T, int Size>
class unitary_left_mult_class
{

private:
    
    multiply_accumulate_class<inA_T,inB_T,accu_T,out_T,2,false,false> MAC;

public:

    unitary_left_mult_class() {};

void CCS_BLOCK(left_mult)(ijm_rxx_idx_t ROW_INDEX,
                          ijm_rxx_idx_t COL_INDEX,
                          inA_T inDataA[Size][Size],
                          inB_T inDataB[Size][Size],
                          out_T outData[Size][Size])
{

    inA_T Arow[2];
    inB_T Bcol[2];

    //Init outdata
    #pragma hls_unroll yes
    INIT_OUTROW: for(int i=0; i<Size; ++i)
    {
        #pragma hls_unroll yes
        INIT_OUTCOL: for(int j=0; j<Size; ++j)
        {
            outData[i][j] = inDataB[i][j];
        }
    }

    //First row
    Arow[0] = inDataA[ROW_INDEX][ROW_INDEX];
    Arow[1] = inDataA[ROW_INDEX][COL_INDEX];

    #pragma hls_unroll yes
    FIRST_ROW: for(int i=0; i<Size; ++i)
    {
        Bcol[0] = inDataB[ROW_INDEX][i];
        Bcol[1] = inDataB[COL_INDEX][i];

        //First row according to ROW_INDEX
        outData[ROW_INDEX][i] = MAC.MAC(Arow,Bcol);
    }

    //Second row
    Arow[0] = inDataA[COL_INDEX][ROW_INDEX];
    Arow[1] = inDataA[COL_INDEX][COL_INDEX];

    #pragma hls_unroll yes
    SECOND_ROW: for(int i=0; i<Size; ++i)
    {
        Bcol[0] = inDataB[ROW_INDEX][i];
        Bcol[1] = inDataB[COL_INDEX][i];

        //Second row according to COL_INDEX
        outData[COL_INDEX][i] = MAC.MAC(Arow,Bcol);
    }

}//left_mult

};//class

//Exploit sparsity in rotation transformation
template<typename inA_T, typename inB_T, typename accu_T, typename out_T, int Size>
class unitary_right_mult_class
{

private:
    
    multiply_accumulate_class<inA_T,inB_T,accu_T,out_T,2,false,true> MAC;

public:

    unitary_right_mult_class() {};

void CCS_BLOCK(right_mult)(ijm_rxx_idx_t ROW_INDEX,
                          ijm_rxx_idx_t COL_INDEX,
                          inA_T inDataA[Size][Size],
                          inB_T inDataB[Size][Size],
                          out_T outData[Size][Size])
{

    inA_T Arow[2];
    inB_T Bcol[2];

    //Init outdata
    #pragma hls_unroll yes
    INIT_OUTROW: for(int i=0; i<Size; ++i)
    {
        #pragma hls_unroll yes
        INIT_OUTCOL: for(int j=0; j<Size; ++j)
        {
            outData[i][j] = inDataA[i][j];
        }
    }

    //First column
    Bcol[0] = inDataB[ROW_INDEX][ROW_INDEX];
    Bcol[1] = inDataB[ROW_INDEX][COL_INDEX];

    #pragma hls_unroll yes
    FIRST_ROW: for(int i=0; i<Size; ++i)
    {
        Arow[0] = inDataA[i][ROW_INDEX];
        Arow[1] = inDataA[i][COL_INDEX];

        //Column according to ROW_INDEX
        outData[i][ROW_INDEX] = MAC.MAC(Arow,Bcol);
    }

    //Second column
    //transpose of inDataB
    Bcol[0] = inDataB[COL_INDEX][ROW_INDEX];
    Bcol[1] = inDataB[COL_INDEX][COL_INDEX];

    #pragma hls_unroll yes
    SECOND_ROW: for(int i=0; i<Size; ++i)
    {
        Arow[0] = inDataA[i][ROW_INDEX];
        Arow[1] = inDataA[i][COL_INDEX];

        //Column according to COL_INDEX
        outData[i][COL_INDEX] = MAC.MAC(Arow,Bcol);
    }

}//left_mult

};//class

#endif //MUSIC_MATRIX_CLASS_H
