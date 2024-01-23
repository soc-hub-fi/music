#ifndef MAX_SEARCH_H_INCLUDED
#define MAX_SEARCH_H_INCLUDED

#include "defs.h"

//Store max. value, row, and column indices
struct element_struct{
    ijm_diag_max_t val;
    ijm_rxx_idx_t r;
    ijm_rxx_idx_t c;
};

//Recuresive template search
//High-level synthesis blue book Ex. 6-22
template<int N> struct max_s{

    static element_struct max(element_struct *a){
        element_struct m0 = max_s<N/2>::max(a);
        element_struct m1 = max_s<N-N/2>::max(a+N/2);
        return m0.val > m1.val ? m0 : m1;
    }
};

template<> struct max_s<1>{
    static element_struct max(element_struct *a){
        return a[0];
    }
};

template<int N> element_struct max(element_struct *a){
    return max_s<N>::max(a);
}

#endif // MAX_SEARCH_H_INCLUDED
