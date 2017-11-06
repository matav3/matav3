/** @file matrix3d.hpp 
    @author deniz
    @description 3d matrix for simple-cubic representation of images
    * x is the fastest varying index, z is the slowest varying index
    * constructor can throw bad_alloc or domain_error exceptions
    * Even for multi-channel images,
    * you can declare your array or struct type, and matrix3d can get it for you.
    * (hence we do NOT assume that codomain is a simple type in constructor)
    * (but we do implement some simpe-type-assuming functions nevertheless) 
    */


namespace bv {
  extern const bool debug;
}

template<class codomain> class matrix3d {
  totSizeType subtotalSize;
  totSizeType totalSize;
  dimSizeType N_z, N_y, N_x;
  // unsigned int N_z, N_y, N_x;
  /* "number of grid points in each dimension"--have to be all POSITIVE */
  dimIndexType min_x, min_y, min_z, max_x, max_y, max_z;
    /* where origin is concerned! (see code conventions) */
  codomain*   x_array;
  codomain*   past_x; // points to location one past the last valid one in x_array
  codomain**  y_array;
  codomain*** z_array;
  codomain*** origin;
  codomain minOfMatrix;
  codomain maxOfMatrix;
public:
  // according to DIG coding conventions of March 2003,
  // place the origin at the center of the 3d lattice.
  // hence, the whole business of specifying hi and lo points is obsolete.

  // can't use valarray, since negative subscripts are prohibited for those.
  // vector is needless, since we don't expect our images to grow

  // constructor can throw bad_alloc or domain_error exceptions
  matrix3d(dimSizeType N_z, dimSizeType N_y, dimSizeType N_x) {
    if(bv::debug) // chose over "bv::debug &&" style to suppress "-Wall" warnings
      std::cout << "Creating a new matrix3d object" << std::endl;
    //  std::cout << sizeof(totSizeType) << " " << sizeof(dimSizeType) << std::endl;
    if(sizeof(totSizeType) < 3 * sizeof(dimSizeType)) {
      std::cerr << "long long needs to be at least 3 times as short.  "
		<< "get a new C++ compiler!" << std::endl;
      throw std::exception();
    }
    // " unsigned long long must be at least 3 times as big as unsigned short.";
    subtotalSize = (static_cast<totSizeType>(N_z)) *
                   (static_cast<totSizeType>(N_y));
    totalSize = (static_cast<totSizeType>(subtotalSize)) * 
                (static_cast<totSizeType>(N_x));
    // was getting errors and worse due to lack of this cast...
    //    if ((N_z<1) || (N_y<1) || (N_x<1)) {
    if (0 == totalSize) {
      throw std::domain_error("Non-positive dimension requested!");
      // std::cerr << "Non-positive dim value passed to matrix3d--quitting" << std::endl;
      //      exit(1); // not the best way to handle this...
    }
    this->N_z = N_z;
    this->N_y = N_y;
    this->N_x = N_x;
    // we call new only 3 times instead of O(N_z*N_y) times
    if(bv::debug) {
      std::cout << "totalSize is " << totalSize << std::endl;
      //      std::cout << "N_z*N_y*N_x is " << (N_z*N_y*N_x) << std::endl;
      std::cout << "x_array needs " << totalSize * sizeof(codomain) 
		<< " 'bytes' (sizeof(char)) or bytes' (sizeof(float))" << std::endl;
    }
    try {
      x_array = new codomain   [totalSize]; // N_z*N_y*N_x
      past_x = x_array + totalSize;
      y_array = new codomain*  [subtotalSize]; // N_z*N_y
      z_array = new codomain** [N_z];
    }
    catch(std::bad_alloc) {
      std::cerr << "Heap memory exhausted by matrix3d::matrix3d()" << std::endl;
      //  if(x_array) delete [] x_array;
      //  if(y_array) delete [] y_array;
      //  if(z_array) delete [] z_array;
      // having the above 3 lines uncommented caused segmentation fault on sleepy
      // std::cerr << "i get here now" << std::endl;
      throw;
    }
    if(bv::debug)
      std::cout << "Memory successfully allocated for matrix3d!" << std::endl;

    // need these casts from dimSizeType to dimIndexType
    min_z = - static_cast<dimIndexType>(N_z/2); 
    min_y = - static_cast<dimIndexType>(N_y/2);
    min_x = - static_cast<dimIndexType>(N_x/2);
    max_z = static_cast<dimIndexType>(N_z-1)/2;
    max_y = static_cast<dimIndexType>(N_y-1)/2;
    max_x = static_cast<dimIndexType>(N_x-1)/2;

    // bv::debug && std::cout << min_z << " " << max_z << std::endl;    
    // the above shall be the inclusive lower- & upper-bounds for our origin
    // now we let everything fall into place (setting up indexing hierarchy)
    codomain** y_ptr = y_array;
    codomain** y_bound = y_array + (N_z*N_y);
    codomain* y_elt = x_array + N_x/2;
    do {
      *y_ptr++ = y_elt;
      y_elt += N_x;
    } while(y_ptr < y_bound);
    codomain*** z_ptr = z_array;
    codomain*** z_bound = z_array + N_z;
    codomain** z_elt = y_array + N_y/2;
    do {
      *z_ptr++ = z_elt;
      z_elt += N_y;
    } while(z_ptr < z_bound);
    origin = z_array + N_z/2;
    if(bv::debug)
      std::cout << "Leaving matrix3d::matrix3d()" << std::endl;
  } // --matrix3d ctor

  codomain** operator[](dimIndexType z) {
    // allows our matrix to be addressed just like an array in C++
    // however, would be (even more) time-expensive to use it exclusively
    // hence we can also return origin directly using matrix3d::getOrigin()
    return origin[z];
  } // --matrix3d::operator[]

  void init_to_zeros() {
    if (!x_array) {
      std::cerr << "3dmatrix::init_to_zeros() : x_array not allocated" << std::endl;
      return;
    }
    codomain *x_bound = x_array + totalSize ; // N_z * N_y * N_x;
    codomain *x_ptr = x_array;
    do {
      *x_ptr++ = 0;
    } while (x_ptr<x_bound);
  } // --matrix3d::init_to_zeros()
  
  bool is_all_zero() {    
    dimIndexType z, y, x;
    codomain*** im = getOrigin();
    for(z=min_z; z<=max_z; z++) {
      for(y=min_y; y<=max_y; y++) {
	for(x=min_x; x<=max_x; x++) {
	  if( 0 != im[z][y][x] )
	    return false;
	} // --for x
      } // --for y
    } // --for z
    return true;
  } // --matrix3d::is_all_zero()


  bool is_all_nonzero() {    
    dimIndexType z, y, x;
    codomain*** im = getOrigin();
    for(z=min_z; z<=max_z; z++) {
      for(y=min_y; y<=max_y; y++) {
				for(x=min_x; x<=max_x; x++) {
	 			 if( 0 == im[z][y][x] )
	    		return false;
	} // --for x
      } // --for y
    } // --for z
    return true;
  } // --matrix3d::is_all_nonzero()

  codomain*** getOrigin() {
    // origin is indexable with z, y, and x
    return origin;
  } // --matrix3d::getOrigin()

  /** return pointer to the (first position of) x_array, useful for
   *  low-level functions that just want to deal with the image vector
   */
  codomain* getFirst() {
    return x_array;
  } // --matrix3d::getFirst()

  /** return pointer to the position past the x_array, useful for
   *  low-level functions that just want to deal with the image vector
   */
  codomain* getPast() {
    return past_x;
  } // --matrix3d::getPast()

  /** copy "sizeof our matrix3d" elements from inbuf to x_array of the matrix3d */
  int copy_from_buffer(codomain *inbuf) {
    codomain *x_bound = x_array + totalSize; // N_z * N_y * N_x;
    codomain *x_ptr = x_array;    
    while(x_ptr != x_bound) {
      *x_ptr++ = *inbuf++;
    }    
    return 0;
  } // --matrix3d::copy_from_buffer()
  
  /** wrapper for reversal of boolean images */
  void reverse_bool() { reverse(1); }
  
  /** 'reverses' an image of discrete domain                  *
   *  converts all zeros to ones and everything else to zeros *
   */
  void reverse_discrete() {
    if (!x_array) {
      std::cerr << "No image to reverse: x_array not allocated" << std::endl;
      return;
    }
    codomain* x_bound = x_array + totalSize; // N_z * N_y * N_x;
    codomain* x_ptr = x_array;
    do {
      if(0 == (*x_ptr))
	*x_ptr = 1;
      else
	*x_ptr = 0;
      x_ptr++;
      //  *x_ptr = maxcodomain - *x_ptr++;
    } while (x_ptr<x_bound);
  } // --matrix3d::reverse_discrete()

  void calcMinMax() {
    // minOfMatrix and maxOfMatrix are data members in matrix3d
    if (!x_array) {
      std::cerr << "No image to reverse: x_array not allocated" << std::endl;
      return;
    }
    codomain *x_bound = x_array + totalSize; // N_z * N_y * N_x;
    codomain *x_ptr = x_array;
    minOfMatrix = *x_ptr;
    maxOfMatrix = *x_ptr;
    x_ptr++;
    while(x_ptr != x_bound) {
      minOfMatrix = (minOfMatrix <= (*x_ptr)) ? minOfMatrix : (*x_ptr);
      maxOfMatrix = (maxOfMatrix >= (*x_ptr)) ? maxOfMatrix : (*x_ptr);
      x_ptr++;
    }
  } // --matrix3d::calcMinMax()

  void calcMinMax(codomain* minOfImageP, codomain* maxOfImageP) {
    calcMinMax();
    *minOfImageP = minOfMatrix;
    *maxOfImageP = maxOfMatrix;
  } // --matrix3d::calcMinMax(codomain*, codomain*)
  
  /* 'reverses' the image by subtracting from maxcodomain */
  void reverse(codomain maxcodomain) {
    // maxcodomain is maximum of image codomain
    if (!x_array) {
      std::cerr << "No image to reverse: x_array not allocated" << std::endl;
      return;
    }
    codomain *x_bound = x_array + totalSize; // N_z * N_y * N_x;
    codomain *x_ptr = x_array;
    do {
      *x_ptr = maxcodomain - *x_ptr++;
    } while (x_ptr<x_bound);
  } // --matrix3d::reverse(maxcodomain)

  // min() and max() sometimes give trouble.
  // but you know what?  we don't even use this anymore.
  void make_rect(dimIndexType cz, dimIndexType cy, dimIndexType cx, 
		 dimIndexType rz, dimIndexType ry, dimIndexType rx, codomain val) {
    if((rz<=0) || (ry<=0) || (rx<=0)) {
      throw std::domain_error("matrix3d::make_rect() got a negative radius argument");
    };
    dimIndexType min_z = max(this->min_z, static_cast<dimIndexType>(cz-rz));
    dimIndexType min_y = max(this->min_y, static_cast<dimIndexType>(cy-ry));
    dimIndexType min_x = max(this->min_x, static_cast<dimIndexType>(cx-rx));
    dimIndexType max_z = min(this->max_z, static_cast<dimIndexType>(cz+rz));		    
    dimIndexType max_y = min(this->max_y, static_cast<dimIndexType>(cy+ry));		    
    dimIndexType max_x = min(this->max_x, static_cast<dimIndexType>(cx+rx));
    dimIndexType z, y, x;
    std::cout << "Filled with " << val << "s--"
	      << "z: " << min_z << ".." << max_z << "; " 
	      << "y: " << min_y << ".." << max_y << "; " 
	      << "x: " << min_x << ".." << max_x << std::endl;    
    for(z = min_z; z <= max_z; z++) {
      for(y = min_y; y <= max_y; y++) {
        for(x = min_x; x <= max_x; x++) {
	  origin[z][y][x]=val;
        }
      }
    } 
  } // --matrix3d::make_rect()
  
  ~matrix3d() {
    if(bv::debug)
      std::cout << "Destroying a matrix3d object" << std::endl;
    if (x_array) delete [] x_array;
    if (y_array) delete [] y_array;
    if (z_array) delete [] z_array;
  } // --matrix3d destructor
}; // --matrix3d<codomain> class definition
