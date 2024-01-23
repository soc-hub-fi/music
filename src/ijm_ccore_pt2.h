// ijm_ccore_pt2.h - CCORE for rotation matrix creation
//
// Author: Tuomas Aaltonen : tuomas.aaltonen@tuni.fi
//

#ifndef IJM_CCORE_PT2_H
#define IJM_CCORE_PT2_H

#include "defs.h"
#include <mc_scverify.h>

class ijm_ccore_pt2
{
    private:

        //Coefficients c and s for rotation matrix
        void findCS(ijm_out_val_cpx_t max_ele, ijm_diag_max_t max_val, ijm_c_cpx_t &c, ijm_s_cpx_t &s)
        {
            ijm_a_t a = 0;
            ijm_b_t b = 0;
            ijm_cs_sqrt_t cs_sqrt = 0;
            ijm_cs_reci_t cs_reci = 0;
            ijm_c_arg_t arg_c = 0;
            ijm_s_arg_t arg_s = 0;

            a = max_ele.real();
            b = max_ele.imag();

            ac_math::ac_sqrt_pwl(max_val,cs_sqrt);
            cs_reci = ac_math::ac_reciprocal_pwl<ijm_cs_reci_t>(cs_sqrt);

            arg_c = a*cs_reci;
            arg_s = b*cs_reci;

            c.set_r(arg_c);
            s.set_i(arg_s);
        }

        //Common trigonometric reciprocal
        void find2TReci(ijm_out_val_cpx_t r_ii,ijm_out_val_cpx_t r_jj,
                        ijm_diag_max_t max_val,ijm_r_diff_t &r_diff, ijm_norm_max_t &norm_max, ijm_reci2t_t &reci2T)
        {
                ijm_denom2t_norm_t den2T_norm = 0;
                ijm_denom2t_diff_t den2T_diff = 0;
                ijm_denom2t_t den2T = 0;
                ijm_sqrt_denom2t_t sqrt_den2T = 0;

                r_diff = r_ii.real()-r_jj.real();
                ac_math::ac_sqrt_pwl(max_val,norm_max);

                den2T_norm = CONSTANT_FOUR*norm_max*norm_max;
                den2T_diff = r_diff*r_diff;
                den2T = den2T_norm+den2T_diff;

                ac_math::ac_sqrt_pwl(den2T,sqrt_den2T);
                reci2T = ac_math::ac_reciprocal_pwl<ijm_reci2t_t>(sqrt_den2T);
        }

        //Sin(2T)
        void findSin2T(ijm_reci2t_t reci2T,ijm_norm_max_t norm_max,
                       ijm_r_diff_t r_diff,ijm_sin2t_t &sin2T)
        {
            ijm_nom_sin2t_t nom_sin2T = 0;
            ijm_sin2t_t sin2T_tmp = 0;

            nom_sin2T = CONSTANT_TWO*norm_max;
            sin2T = nom_sin2T*reci2T;
            sin2T_tmp = sin2T;

            //If r_diff < 0, negate sin(2T)
            if(r_diff[r_diff.length()-1] == 1){

                #pragma hls_unroll yes
                NEG_LOOP: for(int i=0; i<sin2T.length(); i++){
                    sin2T[i] = sin2T[i]^1;
                }
                    sin2T++;
            }
        }

        //Cos(T)
        void findCosT(ijm_reci2t_t reci2T,ijm_r_diff_t r_diff, ijm_cost_t &cosT)
        {
            //Cos(2T)
            ijm_nom_cos2t_t nom_cos2T = 0;
            ijm_cos2t_t cos2T = 0;
            ijm_nom_cost_t nom_cosT = 0;
            ijm_arg_cost_t arg_cosT = 0;

            ac_math::ac_abs(r_diff,nom_cos2T);
            cos2T = nom_cos2T*reci2T;

            //Cos(T)
            nom_cosT = CONSTANT_ONE+cos2T;
            arg_cosT = nom_cosT/CONSTANT_TWO;

            ac_math::ac_sqrt_pwl(arg_cosT,cosT);
        }

        //Sin(T)
        void findSinT(ijm_nom_sint_t nom_sinT, ijm_cost_t cosT, ijm_sint_t &sinT)
        {
            ijm_den_sint_t den_sinT = 0;
            ijm_den_sint_reci_t den_sinT_reci = 0;

            den_sinT = CONSTANT_TWO*cosT;
            den_sinT_reci = ac_math::ac_reciprocal_pwl<ijm_den_sint_reci_t>(den_sinT);
            sinT = nom_sinT*den_sinT_reci;
        }

    public:
        ijm_ccore_pt2() {};

    #pragma hls_design interface ccore
    void CCS_BLOCK(run)(ijm_out_val_cpx_t EIGVAL_MTX[RXX_SIZE][RXX_SIZE], ijm_r_cpx_t R_MTX[RXX_SIZE][RXX_SIZE],
                        ijm_rxx_idx_t row_idx, ijm_rxx_idx_t col_idx, ijm_diag_max_t max_val)
    {
        //Find c,s
        ijm_s_cpx_t s = 0;
        ijm_c_cpx_t c = 0;
        ijm_out_val_cpx_t max_ele = EIGVAL_MTX[row_idx][col_idx];

        findCS(max_ele,max_val,c,s);

        //Find common reciprocal
        ijm_r_diff_t r_diff = 0.0;
        ijm_norm_max_t norm_max = 0.0;
        ijm_reci2t_t reci2T = 0.0;
        find2TReci(EIGVAL_MTX[row_idx][row_idx],EIGVAL_MTX[col_idx][col_idx],
                    max_val,r_diff,norm_max,reci2T);

        //Sin2(T)
        ijm_sin2t_t sin2T = 0.0;
        findSin2T(reci2T,norm_max,r_diff,sin2T);

        //Cos(T)
        ijm_cost_t cosT = 0.0;
        findCosT(reci2T,r_diff,cosT);

        //Sin(T)
        ijm_sint_t sinT = 0.0;
        findSinT(sin2T,cosT,sinT);

        //Init rotation matrix
        #pragma hls_unroll yes
        SET_R_ROWS: for(int i=0; i<RXX_SIZE; i++){
            #pragma hls_unroll yes
            SET_R_COLS: for(int j=0; j<RXX_SIZE; j++){

                R_MTX[i][j] = 0.0;

                if(i==j)
                    R_MTX[i][j].set_r(1.0);
            }
        }

        //Rotation matrix pivot points
        ijm_r_cpx_t sinTc = 0;
        ijm_r_cpx_t sinTs = 0;

        sinTc = sinT*c;
        sinTs = sinT*s;

        R_MTX[row_idx][row_idx] = cosT;
        R_MTX[row_idx][col_idx] = sinTc + sinTs;
        R_MTX[col_idx][row_idx] = -sinTc + sinTs;
        R_MTX[col_idx][col_idx] = cosT;
    }

};

#endif // IJM_CCORE_PT2_H
