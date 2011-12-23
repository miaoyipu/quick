/*
 *  spdfgh.h
 *  new_quick
 *
 *  Created by Yipu Miao on 10/6/11.
 *  Copyright 2011 University of Florida. All rights reserved.
 *
 */

#include "gpu.h"
#include <cuda.h>


#ifndef TEST
__device__
#endif
void PSSS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble Ptempx, QUICKDouble Ptempy,QUICKDouble Ptempz, \
          QUICKDouble WPtempx, QUICKDouble WPtempy, QUICKDouble WPtempz)
{
    VY( 1, 0, mtemp) = Ptempx * VY( 0, 0, mtemp) + WPtempx * VY( 0, 0, mtemp+1);
    VY( 2, 0, mtemp) = Ptempy * VY( 0, 0, mtemp) + WPtempy * VY( 0, 0, mtemp+1);
    VY( 3, 0, mtemp) = Ptempz * VY( 0, 0, mtemp) + WPtempz * VY( 0, 0, mtemp+1);
	return;
}

#ifndef TEST
__device__
#endif
void SSPS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble Qtempx, QUICKDouble Qtempy,QUICKDouble Qtempz, \
          QUICKDouble WQtempx, QUICKDouble WQtempy, QUICKDouble WQtempz)
{
    
    VY( 0, 1, mtemp) = Qtempx * VY( 0, 0, mtemp) + WQtempx * VY( 0, 0, mtemp+1);
    VY( 0, 2, mtemp) = Qtempy * VY( 0, 0, mtemp) + WQtempy * VY( 0, 0, mtemp+1);
    VY( 0, 3, mtemp) = Qtempz * VY( 0, 0, mtemp) + WQtempz * VY( 0, 0, mtemp+1);
    
	return;
}

#ifndef TEST
__device__
#endif
void PSPS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble Ptempx, QUICKDouble Ptempy,QUICKDouble Ptempz, \
          QUICKDouble WPtempx, QUICKDouble WPtempy, QUICKDouble WPtempz, QUICKDouble ABCDtemp)
{
    
    VY( 1, 1, mtemp) = Ptempx * VY( 0, 1, mtemp) + WPtempx * VY( 0, 1, mtemp+1) + ABCDtemp * VY( 0, 0, mtemp+1);
    VY( 2, 1, mtemp) = Ptempy * VY( 0, 1, mtemp) + WPtempy * VY( 0, 1, mtemp+1);
    VY( 3, 1, mtemp) = Ptempz * VY( 0, 1, mtemp) + WPtempz * VY( 0, 1, mtemp+1);
    
    VY( 1, 2, mtemp) = Ptempx * VY( 0, 2, mtemp) + WPtempx * VY( 0, 2, mtemp+1);
    VY( 2, 2, mtemp) = Ptempy * VY( 0, 2, mtemp) + WPtempy * VY( 0, 2, mtemp+1) + ABCDtemp * VY( 0, 0, mtemp+1);
    VY( 3, 2, mtemp) = Ptempz * VY( 0, 2, mtemp) + WPtempz * VY( 0, 2, mtemp+1);
    
    VY( 1, 3, mtemp) = Ptempx * VY( 0, 3, mtemp) + WPtempx * VY( 0, 3, mtemp+1);
    VY( 2, 3, mtemp) = Ptempy * VY( 0, 3, mtemp) + WPtempy * VY( 0, 3, mtemp+1);
    VY( 3, 3, mtemp) = Ptempz * VY( 0, 3, mtemp) + WPtempz * VY( 0, 3, mtemp+1) + ABCDtemp * VY( 0, 0, mtemp+1);
    
	return;
}

#ifndef TEST
__device__
#endif
void DSSS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble Ptempx, QUICKDouble Ptempy,QUICKDouble Ptempz, \
          QUICKDouble WPtempx, QUICKDouble WPtempy, QUICKDouble WPtempz, QUICKDouble ABtemp, QUICKDouble CDcom)
{
    VY( 4, 0, mtemp) = Ptempx * VY( 2, 0, mtemp) + WPtempx * VY( 2, 0, mtemp+1);
    VY( 5, 0, mtemp) = Ptempy * VY( 3, 0, mtemp) + WPtempy * VY( 3, 0, mtemp+1);
    VY( 6, 0, mtemp) = Ptempx * VY( 3, 0, mtemp) + WPtempx * VY( 3, 0, mtemp+1);
    
    VY( 7, 0, mtemp) = Ptempx * VY( 1, 0, mtemp) + WPtempx * VY( 1, 0, mtemp+1)+ ABtemp*(VY( 0, 0, mtemp) - CDcom * VY( 0, 0, mtemp+1));
    VY( 8, 0, mtemp) = Ptempy * VY( 2, 0, mtemp) + WPtempy * VY( 2, 0, mtemp+1)+ ABtemp*(VY( 0, 0, mtemp) - CDcom * VY( 0, 0, mtemp+1));
    VY( 9, 0, mtemp) = Ptempz * VY( 3, 0, mtemp) + WPtempz * VY( 3, 0, mtemp+1)+ ABtemp*(VY( 0, 0, mtemp) - CDcom * VY( 0, 0, mtemp+1));
    
	return;
}

#ifndef TEST
__device__
#endif
void SSDS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble Qtempx, QUICKDouble Qtempy,QUICKDouble Qtempz, \
          QUICKDouble WQtempx, QUICKDouble WQtempy, QUICKDouble WQtempz, QUICKDouble CDtemp, QUICKDouble ABcom)
{
    VY( 0, 4, mtemp) = Qtempx * VY( 0, 2, mtemp) + WQtempx * VY( 0, 2, mtemp+1);
    VY( 0, 5, mtemp) = Qtempy * VY( 0, 3, mtemp) + WQtempy * VY( 0, 3, mtemp+1);
    VY( 0, 6, mtemp) = Qtempx * VY( 0, 3, mtemp) + WQtempx * VY( 0, 3, mtemp+1);
    
    VY( 0, 7, mtemp) = Qtempx * VY( 0, 1, mtemp) + WQtempx * VY( 0, 1, mtemp+1)+ CDtemp*(VY( 0, 0, mtemp) - ABcom * VY( 0, 0, mtemp+1));
    VY( 0, 8, mtemp) = Qtempy * VY( 0, 2, mtemp) + WQtempy * VY( 0, 2, mtemp+1)+ CDtemp*(VY( 0, 0, mtemp) - ABcom * VY( 0, 0, mtemp+1));
    VY( 0, 9, mtemp) = Qtempz * VY( 0, 3, mtemp) + WQtempz * VY( 0, 3, mtemp+1)+ CDtemp*(VY( 0, 0, mtemp) - ABcom * VY( 0, 0, mtemp+1));
    
	return;
}



#ifndef TEST
__device__
#endif
void DSPS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble Qtempx, QUICKDouble Qtempy,QUICKDouble Qtempz, \
          QUICKDouble WQtempx, QUICKDouble WQtempy, QUICKDouble WQtempz, QUICKDouble ABCDtemp)
{
    
    VY( 4, 1, mtemp) = Qtempx * VY( 4, 0, mtemp) + WQtempx * VY( 4, 0, mtemp + 1) + ABCDtemp * VY( 2, 0, mtemp + 1);
    VY( 4, 2, mtemp) = Qtempy * VY( 4, 0, mtemp) + WQtempy * VY( 4, 0, mtemp + 1) + ABCDtemp * VY( 1, 0, mtemp + 1);
    VY( 4, 3, mtemp) = Qtempz * VY( 4, 0, mtemp) + WQtempz * VY( 4, 0, mtemp + 1);
    
    VY( 5, 1, mtemp) = Qtempx * VY( 5, 0, mtemp) + WQtempx * VY( 5, 0, mtemp + 1);
    VY( 5, 2, mtemp) = Qtempy * VY( 5, 0, mtemp) + WQtempy * VY( 5, 0, mtemp + 1) + ABCDtemp * VY( 3, 0, mtemp + 1);
    VY( 5, 3, mtemp) = Qtempz * VY( 5, 0, mtemp) + WQtempz * VY( 5, 0, mtemp + 1) + ABCDtemp * VY( 2, 0, mtemp + 1);
    
    VY( 6, 1, mtemp) = Qtempx * VY( 6, 0, mtemp) + WQtempx * VY( 6, 0, mtemp + 1) + ABCDtemp * VY( 3, 0, mtemp + 1);
    VY( 6, 2, mtemp) = Qtempy * VY( 6, 0, mtemp) + WQtempy * VY( 6, 0, mtemp + 1);
    VY( 6, 3, mtemp) = Qtempz * VY( 6, 0, mtemp) + WQtempz * VY( 6, 0, mtemp + 1) + ABCDtemp * VY( 1, 0, mtemp + 1);
    
    VY( 7, 1, mtemp) = Qtempx * VY( 7, 0, mtemp) + WQtempx * VY( 7, 0, mtemp + 1) + ABCDtemp * VY( 1, 0, mtemp + 1) * 2;
    VY( 7, 2, mtemp) = Qtempy * VY( 7, 0, mtemp) + WQtempy * VY( 7, 0, mtemp + 1);
    VY( 7, 3, mtemp) = Qtempz * VY( 7, 0, mtemp) + WQtempz * VY( 7, 0, mtemp + 1);
    
    VY( 8, 1, mtemp) = Qtempx * VY( 8, 0, mtemp) + WQtempx * VY( 8, 0, mtemp + 1);
    VY( 8, 2, mtemp) = Qtempy * VY( 8, 0, mtemp) + WQtempy * VY( 8, 0, mtemp + 1) + ABCDtemp * VY( 2, 0, mtemp + 1) * 2;
    VY( 8, 3, mtemp) = Qtempz * VY( 8, 0, mtemp) + WQtempz * VY( 8, 0, mtemp + 1);
    
    VY( 9, 1, mtemp) = Qtempx * VY( 9, 0, mtemp) + WQtempx * VY( 9, 0, mtemp + 1);
    VY( 9, 2, mtemp) = Qtempy * VY( 9, 0, mtemp) + WQtempy * VY( 9, 0, mtemp + 1);
    VY( 9, 3, mtemp) = Qtempz * VY( 9, 0, mtemp) + WQtempz * VY( 9, 0, mtemp + 1) + ABCDtemp * VY( 3, 0, mtemp + 1) * 2;            
    
	return;
}

#ifndef TEST
__device__
#endif
void PSDS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble Ptempx, QUICKDouble Ptempy,QUICKDouble Ptempz, \
          QUICKDouble WPtempx, QUICKDouble WPtempy, QUICKDouble WPtempz, QUICKDouble ABCDtemp)
{
    
    
    VY( 1, 4, mtemp) = Ptempx * VY( 0, 4, mtemp) + WPtempx * VY( 0, 4, mtemp + 1) + ABCDtemp * VY( 0, 2, mtemp + 1);
    VY( 2, 4, mtemp) = Ptempy * VY( 0, 4, mtemp) + WPtempy * VY( 0, 4, mtemp + 1) + ABCDtemp * VY( 0, 1, mtemp + 1);
    VY( 3, 4, mtemp) = Ptempz * VY( 0, 4, mtemp) + WPtempz * VY( 0, 4, mtemp + 1);
    
    VY( 1, 5, mtemp) = Ptempx * VY( 0, 5, mtemp) + WPtempx * VY( 0, 5, mtemp + 1);
    VY( 2, 5, mtemp) = Ptempy * VY( 0, 5, mtemp) + WPtempy * VY( 0, 5, mtemp + 1) + ABCDtemp * VY( 0, 3, mtemp + 1);
    VY( 3, 5, mtemp) = Ptempz * VY( 0, 5, mtemp) + WPtempz * VY( 0, 5, mtemp + 1) + ABCDtemp * VY( 0, 2, mtemp + 1);
    
    VY( 1, 6, mtemp) = Ptempx * VY( 0, 6, mtemp) + WPtempx * VY( 0, 6, mtemp + 1) + ABCDtemp * VY( 0, 3, mtemp + 1);
    VY( 2, 6, mtemp) = Ptempy * VY( 0, 6, mtemp) + WPtempy * VY( 0, 6, mtemp + 1);
    VY( 3, 6, mtemp) = Ptempz * VY( 0, 6, mtemp) + WPtempz * VY( 0, 6, mtemp + 1) + ABCDtemp * VY( 0, 1, mtemp + 1);
    
    VY( 1, 7, mtemp) = Ptempx * VY( 0, 7, mtemp) + WPtempx * VY( 0, 7, mtemp + 1) + ABCDtemp * VY( 0, 1, mtemp + 1) * 2;
    VY( 2, 7, mtemp) = Ptempy * VY( 0, 7, mtemp) + WPtempy * VY( 0, 7, mtemp + 1);
    VY( 3, 7, mtemp) = Ptempz * VY( 0, 7, mtemp) + WPtempz * VY( 0, 7, mtemp + 1);
    
    VY( 1, 8, mtemp) = Ptempx * VY( 0, 8, mtemp) + WPtempx * VY( 0, 8, mtemp + 1);
    VY( 2, 8, mtemp) = Ptempy * VY( 0, 8, mtemp) + WPtempy * VY( 0, 8, mtemp + 1) + ABCDtemp * VY( 0, 2, mtemp + 1) * 2;
    VY( 3, 8, mtemp) = Ptempz * VY( 0, 8, mtemp) + WPtempz * VY( 0, 8, mtemp + 1);
    
    VY( 1, 9, mtemp) = Ptempx * VY( 0, 9, mtemp) + WPtempx * VY( 0, 9, mtemp + 1);
    VY( 2, 9, mtemp) = Ptempy * VY( 0, 9, mtemp) + WPtempy * VY( 0, 9, mtemp + 1);
    VY( 3, 9, mtemp) = Ptempz * VY( 0, 9, mtemp) + WPtempz * VY( 0, 9, mtemp + 1) + ABCDtemp * VY( 0, 3, mtemp + 1) * 2;    
    
	return;
}

#ifndef TEST
__device__
#endif
void DSDS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble Ptempx, QUICKDouble Ptempy,QUICKDouble Ptempz, \
          QUICKDouble WPtempx, QUICKDouble WPtempy, QUICKDouble WPtempz, QUICKDouble ABCDtemp, QUICKDouble ABtemp, QUICKDouble CDcom)
{
    
    VY( 4, 4, mtemp) = Ptempx * VY( 2, 4, mtemp) + WPtempx * VY( 2, 4, mtemp+1) + ABCDtemp * VY( 2, 2, mtemp+1);
    VY( 4, 5, mtemp) = Ptempx * VY( 2, 5, mtemp) + WPtempx * VY( 2, 5, mtemp+1);
    VY( 4, 6, mtemp) = Ptempx * VY( 2, 6, mtemp) + WPtempx * VY( 2, 6, mtemp+1) + ABCDtemp * VY( 2, 3, mtemp+1);
    VY( 4, 7, mtemp) = Ptempx * VY( 2, 7, mtemp) + WPtempx * VY( 2, 7, mtemp+1) + 2 * ABCDtemp * VY( 2, 1, mtemp+1);
    VY( 4, 8, mtemp) = Ptempx * VY( 2, 8, mtemp) + WPtempx * VY( 2, 8, mtemp+1);
    VY( 4, 9, mtemp) = Ptempx * VY( 2, 9, mtemp) + WPtempx * VY( 2, 9, mtemp+1);
    
    VY( 5, 4, mtemp) = Ptempy * VY( 3, 4, mtemp) + WPtempy * VY( 3, 4, mtemp+1) + ABCDtemp * VY( 3, 1, mtemp+1);
    VY( 5, 5, mtemp) = Ptempy * VY( 3, 5, mtemp) + WPtempy * VY( 3, 5, mtemp+1) + ABCDtemp * VY( 3, 3, mtemp+1);
    VY( 5, 6, mtemp) = Ptempy * VY( 3, 6, mtemp) + WPtempy * VY( 3, 6, mtemp+1);
    VY( 5, 7, mtemp) = Ptempy * VY( 3, 7, mtemp) + WPtempy * VY( 3, 7, mtemp+1);
    VY( 5, 8, mtemp) = Ptempy * VY( 3, 8, mtemp) + WPtempy * VY( 3, 8, mtemp+1) + 2 * ABCDtemp * VY( 3, 2, mtemp+1);
    VY( 5, 9, mtemp) = Ptempy * VY( 3, 9, mtemp) + WPtempy * VY( 3, 9, mtemp+1);
    
    VY( 6, 4, mtemp) = Ptempx * VY( 3, 4, mtemp) + WPtempx * VY( 3, 4, mtemp+1) + ABCDtemp * VY( 3, 2, mtemp+1);
    VY( 6, 5, mtemp) = Ptempx * VY( 3, 5, mtemp) + WPtempx * VY( 3, 5, mtemp+1);
    VY( 6, 6, mtemp) = Ptempx * VY( 3, 6, mtemp) + WPtempx * VY( 3, 6, mtemp+1) + ABCDtemp * VY( 3, 3, mtemp+1);
    VY( 6, 7, mtemp) = Ptempx * VY( 3, 7, mtemp) + WPtempx * VY( 3, 7, mtemp+1) + 2 * ABCDtemp * VY( 3, 1, mtemp+1);
    VY( 6, 8, mtemp) = Ptempx * VY( 3, 8, mtemp) + WPtempx * VY( 3, 8, mtemp+1);
    VY( 6, 9, mtemp) = Ptempx * VY( 3, 9, mtemp) + WPtempx * VY( 3, 9, mtemp+1);
    
    VY( 7, 4, mtemp) = Ptempx * VY( 1, 4, mtemp) + WPtempx * VY( 1, 4, mtemp+1) +  ABtemp * (VY( 0, 4,mtemp) - CDcom * VY( 0, 4,mtemp+1)) + ABCDtemp * VY( 1, 2, mtemp+1);
    VY( 7, 5, mtemp) = Ptempx * VY( 1, 5, mtemp) + WPtempx * VY( 1, 5, mtemp+1) +  ABtemp * (VY( 0, 5,mtemp) - CDcom * VY( 0, 5,mtemp+1));
    VY( 7, 6, mtemp) = Ptempx * VY( 1, 6, mtemp) + WPtempx * VY( 1, 6, mtemp+1) +  ABtemp * (VY( 0, 6,mtemp) - CDcom * VY( 0, 6,mtemp+1)) + ABCDtemp * VY( 1, 3, mtemp+1);
    VY( 7, 7, mtemp) = Ptempx * VY( 1, 7, mtemp) + WPtempx * VY( 1, 7, mtemp+1) +  ABtemp * (VY( 0, 7,mtemp) - CDcom * VY( 0, 7,mtemp+1)) + 2 * ABCDtemp * VY( 1, 1, mtemp+1);
    VY( 7, 8, mtemp) = Ptempx * VY( 1, 8, mtemp) + WPtempx * VY( 1, 8, mtemp+1) +  ABtemp * (VY( 0, 8,mtemp) - CDcom * VY( 0, 8,mtemp+1));
    VY( 7, 9, mtemp) = Ptempx * VY( 1, 9, mtemp) + WPtempx * VY( 1, 9, mtemp+1) +  ABtemp * (VY( 0, 9,mtemp) - CDcom * VY( 0, 9,mtemp+1));
    
    
    VY( 8, 4, mtemp) = Ptempy * VY( 2, 4, mtemp) + WPtempy * VY( 2, 4, mtemp+1) +  ABtemp * (VY( 0, 4,mtemp) - CDcom * VY( 0, 4,mtemp+1)) + ABCDtemp * VY( 2, 1, mtemp+1);
    VY( 8, 5, mtemp) = Ptempy * VY( 2, 5, mtemp) + WPtempy * VY( 2, 5, mtemp+1) +  ABtemp * (VY( 0, 5,mtemp) - CDcom * VY( 0, 5,mtemp+1)) + ABCDtemp * VY( 2, 3, mtemp+1);
    VY( 8, 6, mtemp) = Ptempy * VY( 2, 6, mtemp) + WPtempy * VY( 2, 6, mtemp+1) +  ABtemp * (VY( 0, 6,mtemp) - CDcom * VY( 0, 6,mtemp+1));
    VY( 8, 7, mtemp) = Ptempy * VY( 2, 7, mtemp) + WPtempy * VY( 2, 7, mtemp+1) +  ABtemp * (VY( 0, 7,mtemp) - CDcom * VY( 0, 7,mtemp+1));
    VY( 8, 8, mtemp) = Ptempy * VY( 2, 8, mtemp) + WPtempy * VY( 2, 8, mtemp+1) +  ABtemp * (VY( 0, 8,mtemp) - CDcom * VY( 0, 8,mtemp+1)) + 2 * ABCDtemp * VY( 2, 2, mtemp+1);
    VY( 8, 9, mtemp) = Ptempy * VY( 2, 9, mtemp) + WPtempy * VY( 2, 9, mtemp+1) +  ABtemp * (VY( 0, 9,mtemp) - CDcom * VY( 0, 9,mtemp+1));
    
    VY( 9, 4, mtemp) = Ptempz * VY( 3, 4, mtemp) + WPtempz * VY( 3, 4, mtemp+1) +  ABtemp * (VY( 0, 4,mtemp) - CDcom * VY( 0, 4,mtemp+1));
    VY( 9, 5, mtemp) = Ptempz * VY( 3, 5, mtemp) + WPtempz * VY( 3, 5, mtemp+1) +  ABtemp * (VY( 0, 5,mtemp) - CDcom * VY( 0, 5,mtemp+1)) + ABCDtemp * VY( 3, 2, mtemp+1);
    VY( 9, 6, mtemp) = Ptempz * VY( 3, 6, mtemp) + WPtempz * VY( 3, 6, mtemp+1) +  ABtemp * (VY( 0, 6,mtemp) - CDcom * VY( 0, 6,mtemp+1)) + ABCDtemp * VY( 3, 1, mtemp+1);
    VY( 9, 7, mtemp) = Ptempz * VY( 3, 7, mtemp) + WPtempz * VY( 3, 7, mtemp+1) +  ABtemp * (VY( 0, 7,mtemp) - CDcom * VY( 0, 7,mtemp+1));
    VY( 9, 8, mtemp) = Ptempz * VY( 3, 8, mtemp) + WPtempz * VY( 3, 8, mtemp+1) +  ABtemp * (VY( 0, 8,mtemp) - CDcom * VY( 0, 8,mtemp+1));
    VY( 9, 9, mtemp) = Ptempz * VY( 3, 9, mtemp) + WPtempz * VY( 3, 9, mtemp+1) +  ABtemp * (VY( 0, 9,mtemp) - CDcom * VY( 0, 9,mtemp+1)) + 2 * ABCDtemp * VY( 3, 3, mtemp+1);
    
	return;
}

