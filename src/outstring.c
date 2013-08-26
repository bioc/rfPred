#include <R.h>
#include <Rdefines.h>
#include <R_ext/Error.h>

SEXP outstring ( SEXP x, SEXP nfield, SEXP sepfield ) {

int ix,nx=length(x),ny,sep,is,ns,iy,ib;
const char *pts;
char ptb[1000];
SEXP s;
ny=INTEGER(nfield)[0];
sep=CHAR(STRING_ELT(sepfield,0))[0];
PROTECT(s=allocMatrix(STRSXP, nx,ny));
for (ix = 0; ix < nx; ix++) {
  pts=CHAR(STRING_ELT(x,ix));
  ns=strlen(pts);
  iy=0;
  ib=0;
  if (ns==0) {
   ptb[ib]=0;
   SET_STRING_ELT(s,(iy*nx)+ix, mkChar(ptb));
  }
  else {
    for (is=0;is<ns;is++) {
      ptb[ib]=pts[is];
      if(ptb[ib]==sep) {
        ptb[ib]=0;
        SET_STRING_ELT(s,(iy*nx)+ix, mkChar(ptb));
        ib=0;
        iy++;
      }
      else ib++;
    }
  }
  ptb[ib]=0;
  SET_STRING_ELT(s,(iy*nx)+ix, mkChar(ptb));
}
UNPROTECT(1);
return(s);
return R_NilValue;
}
