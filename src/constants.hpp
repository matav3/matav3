#ifndef _ABSDIF_DEFINED
#define _ABSDIF_DEFINED
// C-style macro for type catholicism... can't really write an inline function for this.
#define ABSDIF(a,b)  ((a>=b) ? (a-b) : (b-a))
// safe to use whenever a and b are of the same signed or unsigned integral type
#endif
#ifndef _DEBUG_DEFINED
#define _DEBUG_DEFINED
// global variable 'debug' clashes with Qt stuff, so stick it in a namespace...
namespace bv { 
  const bool debug = false;
}
#endif

