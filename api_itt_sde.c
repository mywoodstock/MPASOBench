#ifdef USE_VTUNE
#include "ittnotify.h"
#endif

void fortran_sde_start()
{
#ifdef USE_SDE
   __SSC_MARK(0x111);
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
#ifdef USE_VTUNE
  __itt_resume();
#endif
}

void fortran_itt_pause()
{
#ifdef USE_VTUNE
  __itt_pause();
#endif
}
