// nsr_class.h - Noise Subspace Recovery (NSR) hierarchical block
//
// Author: tuomas.aaltonen@tuni.fi
//

#ifndef NSR_CLASS_H
#define NSR_CLASS_H

#include "defs.h"

#include "ijm_class.h"
#include "nss_class.h"

#include <mc_scverify.h>

class nsr_top_class
{

    private:

        //Sub-block class instances
        IJM_class ijm_inst;
        NSS_class nss_inst;

        //Interconnect channel
        ac_channel<eigenStruct_t> Eigen_chan;

    public:

        nsr_top_class() {};


#ifdef MUSIC_DEBUG
    #pragma hls_design interface
    void CCS_BLOCK(run)(ac_channel<rxxStruct_t> &IJMIn,
            ac_channel<eigenStruct_t> &IJM_Out_Dbg_chan,
            ac_channel<nssStruct_t> &NSROut,
            ac_channel<nssStruct_t> &NSR_Out_Dbg_chan)
#else
    #pragma hls_design interface
    void CCS_BLOCK(run)(ac_channel<rxxStruct_t> &IJMIn,
            ac_channel<nssStruct_t> &NSROut)
#endif // MUSIC_DEBUG
    {
#ifdef MUSIC_DEBUG
        ijm_inst.run(IJMIn,IJM_Out_Dbg_chan,Eigen_chan);
        nss_inst.run(Eigen_chan,NSR_Out_Dbg_chan,NSROut);
#else
        ijm_inst.run(IJMIn,Eigen_chan);
        nss_inst.run(Eigen_chan,NSROut);

#endif // MUSIC_DEBUG
    }   //run
};//nsr_top_class

#endif // NSR_CLASS_H
