// music_top_class.h - Top level hierarchical block
//
// Author: Tuomas Aaltonen : tuomas.aaltonen@tuni.fi
// Tampere University - Unit of Computing Sciences - 2024
//

#ifndef MUSIC_TOP_CLASS_H
#define MUSIC_TOP_CLASS_H

#include "defs.h"

#include "scm_class.h"
#include "ijm_class.h"
#include "sms_class.h"
#include "nss_class.h"
#include "nsr_class.h"

#include <mc_scverify.h>

#pragma hls_design top
class music_top_class
{

private:

    //Sub-block class instances
    SCM_class scm_inst;
    IJM_class ijm_inst;
    SMS_class sms_inst;
    nsr_top_class nsr_inst; //Hierarchy level

    //Interconnect channels
    ac_channel<rxxStruct_t> Rxx_chan;
    ac_channel<nssStruct_t> Un_chan;

public:

    music_top_class() {};

#ifdef MUSIC_DEBUG
    #pragma hls_design interface
    void CCS_BLOCK(run)(ac_channel<inStruct_t> &SCMIn,
            ac_channel<rxxStruct_t> &SCM_Out_Dbg_chan,
            ac_channel<eigenStruct_t> &IJM_Out_Dbg_chan,
            ac_channel<nssStruct_t> &NSR_Out_Dbg_chan,
            ac_channel<sms_out_t> &SMSOut)
#else
    #pragma hls_design interface
    void CCS_BLOCK(run)(ac_channel<inStruct_t> &SCMIn,
            ac_channel<sms_out_t> &SMSOut)
#endif //MUSIC_DEBUG
    {
#ifdef MUSIC_DEBUG
        scm_inst.run(SCMIn,Rxx_chan,SCM_Out_Dbg_chan);
        nsr_inst.run(Rxx_chan,IJM_Out_Dbg_chan,Un_chan,NSR_Out_Dbg_chan);
        sms_inst.run(Un_chan,SMSOut);
#else
        scm_inst.run(SCMIn,Rxx_chan);
        nsr_inst.run(Rxx_chan,Un_chan);
        sms_inst.run(Un_chan,SMSOut);
#endif //MUSIC_DEBUG
    }// run

}; //music_top_class

#endif //MUSIC_TOP_CLASS
