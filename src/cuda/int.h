#include "include/h_all_files.h"
#include "include/h_all_subroutines.h"

__device__ void vertical_spdf(int I, int J, int K, int L,
                         QUICKDouble* YVerticalTemp, QUICKDouble* store,
                         QUICKDouble Ptempx, QUICKDouble Ptempy, QUICKDouble Ptempz,  \
                         QUICKDouble WPtempx,QUICKDouble WPtempy,QUICKDouble WPtempz, \
                         QUICKDouble Qtempx, QUICKDouble Qtempy, QUICKDouble Qtempz,  \
                         QUICKDouble WQtempx,QUICKDouble WQtempy,QUICKDouble WQtempz, \
                         QUICKDouble ABCDtemp,QUICKDouble ABtemp, \
                         QUICKDouble CDtemp, QUICKDouble ABcom, QUICKDouble CDcom)
{
    
    if ((I+J) >=  0 && (K+L) >= 5) {
        h_0_5(YVerticalTemp, store, \
              Ptempx, Ptempy, Ptempz, WPtempx, WPtempy, WPtempz, Qtempx, Qtempy, Qtempz,  \
              WQtempx, WQtempy, WQtempz, ABCDtemp, ABtemp, CDtemp,  ABcom, CDcom);
    }
    if ((I+J) >=  0 && (K+L) >= 6) {
        h_0_6(YVerticalTemp, store, \
              Ptempx, Ptempy, Ptempz, WPtempx, WPtempy, WPtempz, Qtempx, Qtempy, Qtempz,  \
              WQtempx, WQtempy, WQtempz, ABCDtemp, ABtemp, CDtemp,  ABcom, CDcom);
    }
    
    if ((I+J) >=  1 && (K+L) >= 5) {
        h_1_5(YVerticalTemp, store, \
              Ptempx, Ptempy, Ptempz, WPtempx, WPtempy, WPtempz, Qtempx, Qtempy, Qtempz,  \
              WQtempx, WQtempy, WQtempz, ABCDtemp, ABtemp, CDtemp,  ABcom, CDcom);
    }
    if ((I+J) >=  1 && (K+L) >= 6) {
        h_1_6(YVerticalTemp, store, \
              Ptempx, Ptempy, Ptempz, WPtempx, WPtempy, WPtempz, Qtempx, Qtempy, Qtempz,  \
              WQtempx, WQtempy, WQtempz, ABCDtemp, ABtemp, CDtemp,  ABcom, CDcom);
    }
    
    if ((I+J) >=  2 && (K+L) >= 5) {
        h_2_5(YVerticalTemp, store, \
              Ptempx, Ptempy, Ptempz, WPtempx, WPtempy, WPtempz, Qtempx, Qtempy, Qtempz,  \
              WQtempx, WQtempy, WQtempz, ABCDtemp, ABtemp, CDtemp,  ABcom, CDcom);
    }
    if ((I+J) >=  2 && (K+L) >= 6) {
        h_2_6(YVerticalTemp, store, \
              Ptempx, Ptempy, Ptempz, WPtempx, WPtempy, WPtempz, Qtempx, Qtempy, Qtempz,  \
              WQtempx, WQtempy, WQtempz, ABCDtemp, ABtemp, CDtemp,  ABcom, CDcom);
    }
    
    if ((I+J) >=  3 && (K+L) >= 5) {
        h_3_5(YVerticalTemp, store, \
              Ptempx, Ptempy, Ptempz, WPtempx, WPtempy, WPtempz, Qtempx, Qtempy, Qtempz,  \
              WQtempx, WQtempy, WQtempz, ABCDtemp, ABtemp, CDtemp,  ABcom, CDcom);
    }
    if ((I+J) >=  3 && (K+L) >= 6) {
        h_3_6(YVerticalTemp, store, \
              Ptempx, Ptempy, Ptempz, WPtempx, WPtempy, WPtempz, Qtempx, Qtempy, Qtempz,  \
              WQtempx, WQtempy, WQtempz, ABCDtemp, ABtemp, CDtemp,  ABcom, CDcom);
    }
    
    if ((I+J) >=  4 && (K+L) >= 5) {
        h_4_5(YVerticalTemp, store, \
              Ptempx, Ptempy, Ptempz, WPtempx, WPtempy, WPtempz, Qtempx, Qtempy, Qtempz,  \
              WQtempx, WQtempy, WQtempz, ABCDtemp, ABtemp, CDtemp,  ABcom, CDcom);
    }
    if ((I+J) >=  4 && (K+L) >= 6) {
        h_4_6(YVerticalTemp, store, \
              Ptempx, Ptempy, Ptempz, WPtempx, WPtempy, WPtempz, Qtempx, Qtempy, Qtempz,  \
              WQtempx, WQtempy, WQtempz, ABCDtemp, ABtemp, CDtemp,  ABcom, CDcom);
    }
}






__device__ void vertical_spdf2(int I, int J, int K, int L,
                              QUICKDouble* YVerticalTemp, QUICKDouble* store,
                              QUICKDouble Ptempx, QUICKDouble Ptempy, QUICKDouble Ptempz,  \
                              QUICKDouble WPtempx,QUICKDouble WPtempy,QUICKDouble WPtempz, \
                              QUICKDouble Qtempx, QUICKDouble Qtempy, QUICKDouble Qtempz,  \
                              QUICKDouble WQtempx,QUICKDouble WQtempy,QUICKDouble WQtempz, \
                              QUICKDouble ABCDtemp,QUICKDouble ABtemp, \
                               QUICKDouble CDtemp, QUICKDouble ABcom, QUICKDouble CDcom)
{
    if ((I+J) >=  5 && (K+L) >= 0) {
        h_5_0(YVerticalTemp, store, \
              Ptempx, Ptempy, Ptempz, WPtempx, WPtempy, WPtempz, Qtempx, Qtempy, Qtempz,  \
              WQtempx, WQtempy, WQtempz, ABCDtemp, ABtemp, CDtemp,  ABcom, CDcom);
    }
    if ((I+J) >=  5 && (K+L) >= 1) {
        h_5_1(YVerticalTemp, store, \
              Ptempx, Ptempy, Ptempz, WPtempx, WPtempy, WPtempz, Qtempx, Qtempy, Qtempz,  \
              WQtempx, WQtempy, WQtempz, ABCDtemp, ABtemp, CDtemp,  ABcom, CDcom);
    }
    if ((I+J) >=  5 && (K+L) >= 2) {
        h_5_2(YVerticalTemp, store, \
              Ptempx, Ptempy, Ptempz, WPtempx, WPtempy, WPtempz, Qtempx, Qtempy, Qtempz,  \
              WQtempx, WQtempy, WQtempz, ABCDtemp, ABtemp, CDtemp,  ABcom, CDcom);
    }
    if ((I+J) >=  5 && (K+L) >= 3) {
        h_5_3(YVerticalTemp, store, \
              Ptempx, Ptempy, Ptempz, WPtempx, WPtempy, WPtempz, Qtempx, Qtempy, Qtempz,  \
              WQtempx, WQtempy, WQtempz, ABCDtemp, ABtemp, CDtemp,  ABcom, CDcom);
    }
    if ((I+J) >=  5 && (K+L) >= 4) {
        h_5_4(YVerticalTemp, store, \
              Ptempx, Ptempy, Ptempz, WPtempx, WPtempy, WPtempz, Qtempx, Qtempy, Qtempz,  \
              WQtempx, WQtempy, WQtempz, ABCDtemp, ABtemp, CDtemp,  ABcom, CDcom);
    }
    if ((I+J) >=  6 && (K+L) >= 0) {
        h_6_0(YVerticalTemp, store, \
              Ptempx, Ptempy, Ptempz, WPtempx, WPtempy, WPtempz, Qtempx, Qtempy, Qtempz,  \
              WQtempx, WQtempy, WQtempz, ABCDtemp, ABtemp, CDtemp,  ABcom, CDcom);
    }
    if ((I+J) >=  6 && (K+L) >= 1) {
        h_6_1(YVerticalTemp, store, \
              Ptempx, Ptempy, Ptempz, WPtempx, WPtempy, WPtempz, Qtempx, Qtempy, Qtempz,  \
              WQtempx, WQtempy, WQtempz, ABCDtemp, ABtemp, CDtemp,  ABcom, CDcom);
    }
    if ((I+J) >=  6 && (K+L) >= 2) {
        h_6_2(YVerticalTemp, store, \
              Ptempx, Ptempy, Ptempz, WPtempx, WPtempy, WPtempz, Qtempx, Qtempy, Qtempz,  \
              WQtempx, WQtempy, WQtempz, ABCDtemp, ABtemp, CDtemp,  ABcom, CDcom);
    }
    if ((I+J) >=  6 && (K+L) >= 3) {
        h_6_3(YVerticalTemp, store, \
              Ptempx, Ptempy, Ptempz, WPtempx, WPtempy, WPtempz, Qtempx, Qtempy, Qtempz,  \
              WQtempx, WQtempy, WQtempz, ABCDtemp, ABtemp, CDtemp,  ABcom, CDcom);
    }
    if ((I+J) >=  6 && (K+L) >= 4) {
        h_6_4(YVerticalTemp, store, \
              Ptempx, Ptempy, Ptempz, WPtempx, WPtempy, WPtempz, Qtempx, Qtempy, Qtempz,  \
              WQtempx, WQtempy, WQtempz, ABCDtemp, ABtemp, CDtemp,  ABcom, CDcom);
    }
}






__device__ void vertical_spdf3(int I, int J, int K, int L,
                               QUICKDouble* YVerticalTemp, QUICKDouble* store,
                               QUICKDouble Ptempx, QUICKDouble Ptempy, QUICKDouble Ptempz,  \
                               QUICKDouble WPtempx,QUICKDouble WPtempy,QUICKDouble WPtempz, \
                               QUICKDouble Qtempx, QUICKDouble Qtempy, QUICKDouble Qtempz,  \
                               QUICKDouble WQtempx,QUICKDouble WQtempy,QUICKDouble WQtempz, \
                               QUICKDouble ABCDtemp,QUICKDouble ABtemp, \
                               QUICKDouble CDtemp, QUICKDouble ABcom, QUICKDouble CDcom)
{
    
    if ((I+J) >=  5 && (K+L) >= 5) {
        h_5_5(YVerticalTemp, store, \
              Ptempx, Ptempy, Ptempz, WPtempx, WPtempy, WPtempz, Qtempx, Qtempy, Qtempz,  \
              WQtempx, WQtempy, WQtempz, ABCDtemp, ABtemp, CDtemp,  ABcom, CDcom);
    }
    
    if ((I+J) >=  5 && (K+L) >= 6) {
        h_5_6(YVerticalTemp, store, \
              Ptempx, Ptempy, Ptempz, WPtempx, WPtempy, WPtempz, Qtempx, Qtempy, Qtempz,  \
              WQtempx, WQtempy, WQtempz, ABCDtemp, ABtemp, CDtemp,  ABcom, CDcom);
    }
}


__device__ void vertical_spdf4(int I, int J, int K, int L,
                               QUICKDouble* YVerticalTemp, QUICKDouble* store,
                               QUICKDouble Ptempx, QUICKDouble Ptempy, QUICKDouble Ptempz,  \
                               QUICKDouble WPtempx,QUICKDouble WPtempy,QUICKDouble WPtempz, \
                               QUICKDouble Qtempx, QUICKDouble Qtempy, QUICKDouble Qtempz,  \
                               QUICKDouble WQtempx,QUICKDouble WQtempy,QUICKDouble WQtempz, \
                               QUICKDouble ABCDtemp,QUICKDouble ABtemp, \
                               QUICKDouble CDtemp, QUICKDouble ABcom, QUICKDouble CDcom)
{
    
    if ((I+J) >=  6 && (K+L) >= 5) {
        h_6_5(YVerticalTemp, store, \
              Ptempx, Ptempy, Ptempz, WPtempx, WPtempy, WPtempz, Qtempx, Qtempy, Qtempz,  \
              WQtempx, WQtempy, WQtempz, ABCDtemp, ABtemp, CDtemp,  ABcom, CDcom);
    }
    if ((I+J) >=  6 && (K+L) >= 6) {
        h_6_6(YVerticalTemp, store, \
              Ptempx, Ptempy, Ptempz, WPtempx, WPtempy, WPtempz, Qtempx, Qtempy, Qtempz,  \
              WQtempx, WQtempy, WQtempz, ABCDtemp, ABtemp, CDtemp,  ABcom, CDcom);
    }
}




