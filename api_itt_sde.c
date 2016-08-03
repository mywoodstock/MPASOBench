#include "ittnotify.h"

void fortran_sde_start()
{
#ifdef USE_SDE
   __SSC_MARK(0x111);
#else
#warning Not using SDE markers
#endif
}  

void fortran_sde_stop()
{
#ifdef USE_SDE
  __SSC_MARK(0x222);
#endif
}  

void fortran_itt_resume()
{
  __itt_resume();
}

void fortran_itt_pause()
{
  __itt_pause();
}
