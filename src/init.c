#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void calcmixmvdens(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void calcmixtdens(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"calcmixmvdens", (DL_FUNC) &calcmixmvdens,  9},
    {"calcmixtdens",  (DL_FUNC) &calcmixtdens,  10},
    {NULL, NULL, 0}
};

void R_init_iterLap(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

