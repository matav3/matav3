/************************************************************************
* @file image3d.cpp                                                     *
* @date June 2003---2005.06.22                                          *
* @author deniz <sarioz at acm dot org>                                 *
* @language ISO C++ 1998, using the C++ Standard Library profusely      *
*                                                                       *
* This class template file is part of Deniz's research package          *
*                                                                       *
* Now just concentrates on history trees                                *
* Removed some features enabling to do betti number-persistence         *
*     in order to to reduce bloat                                       *
*   the file with older functionality is: image3d.cpp.bloated_aug_15    *  
*                                                                       *
* main() is left to be defined somewhere else, this is like a library   *
*                                                                       *
* You have to declare a new image3d from your driver                    *
* You probably want to load a DIGFile using load                        *
*   (although you are entitled to 'manually' edit your image3d in code) *
*   (or passing copy_from_buffer() a buffer of appropriate type)        *
*                                                                       *
* We completely assume Simple Cubic grid for voronoi images and use     *
*   6-connectivity for the 'ones', 26-connectivity for the 'zeros'.     *
* We also assume single-channel (grayscale), non-complex (real) data    *
*                                                                       *
* All matrix3d objects are to be manipulated through image3d objects    *
*   for the sake of safety                                              *
*   (we could define the former as an inner class of the latter)        *
*                                                                       *
* ALL applicable functions and array indices take                       *
*        first z, then y, and then x.                                   *
*                                                                       *
* Author's notes to self: make comments fit C++-doc format better       *
*   DON'T comment out unused functions--they do NOT cause code bloat!   *
*                                                                       *
* Copyrighted under the GNU General Public License by:                  *
*                                                                       *
* Discrete Imaging and Graphics Group (DIG)                             *
* Department of Computer Science                                        *
* The Graduate Center -- City University of New York                    *
* New York, NY, USA                                                     *
* http://www.dig.cs.gc.cuny.edu                                         *
*                                                                       *
************************************************************************/


#include <iostream>
#include <fstream>
#include <iomanip> 
using namespace std;

double kaiser_value(double r, double a, double alpha, int m);

namespace bv {
  extern const bool debug;
}
// const bool debug26 = false; // whether or not to debug 26-connectivity


/** 'wrapper' class for images matrix3d (with some more metadata than matrix3d) 
  * however, note that more than one matrix3d can be created from an image3d.
  * this class deals with I/O
  */

template<class codomain> class image3d {
  bool matrix_set;
  bool minmax_set;
  bettiType b0;
  bettiType b1;
  bettiType b2;
  eulerType E; // I WAS WRONG! the Euler number of a 3d object /can/ be NEGATIVE 
  // (and its calculation can turn thorny, so let's just use long int)
  unsigned smallest; // smallest size for (6) connected components (of 1s)
  codomain maxValInImage, minValInImage, threshold;

  int dym;
  int dym2;
  dsf_for_tree_gen* FGDSF;  // ForeGround DSF for tree generation (related to b0)
  dsf_for_tree_bgnd_gen* BGDSF;  // BackGround DSF for tree generation (related to b2)
  historytree* HTP;
  historytree_bgnd* HTBP;
  typedef std::vector<coord> coordlist;
  typedef std::map<codomain, coordlist> list_of_lists;
  // cannot be defined in BettiVoxel.h due to variability of codomain
  list_of_lists all_voxels;
public:
  dimSizeType N_x, N_y, N_z; // "dimensions" of the image
	codomain minOfImage, maxOfImage;
  dimIndexType min_x, min_y, min_z, max_x, max_y, max_z; // inclusive bounds
  /* "centered" at (0,0,0); s.t. if N_i is even, |min_i| == 1 + max_i */
protected:
  typedef std::pair<codomain, coord> myPairType;
  // cannot be defined in BettiVoxel.h due to variability of codomain
  matrix3d<codomain>* image_data;
  /* to be loaded in from a DIGFile or some other kind of buffer    */
  matrix3d<labelType>* multilabel_matrix;
  /* a new multilabel_matrix can be created at each threshold value * 
   * and freed before the next one...                               *
   * but would be a bit inefficient all those heap memory requests  *
   * so just repopulate for each threshold value                    *
   * free the space when done with all of them.                     *
   *                                                                   *
   * now, it's IMPERATIVE to keep a multilabel_matrix (for each threshold) * 
   * due to filtering out small components by 1s-counter                   *
   *                                                                 *
   * The incremental algorithms suggested by Yung on December 9 2003       *
   *   DEFINITELY require multilabel_matrix to stay there once it's there! *
   */

  void init() { 
    minmax_set = false; 
    matrix_set = false;
    threshold = 0; // contingency: works for boolean image
    smallest = 1; // contingency: components have to be at least size 1 (of course)
    image_data = 0; // as in NULL
    multilabel_matrix = 0; // as in NULL
    // also set all betti numbers and euler number to zero
    b0 = b1 = b2 = 0;
    E = 0;
    //    myDSFpb0 = 0;
    //    myDSFpb2 = 0;
    FGDSF = 0;
    BGDSF = 0;
    HTP = 0;
    HTBP = 0;
    all_voxels.clear();
  } // --image3d::init()



public:

  //    we do NOT want our image_data to be exposed to the outside
  codomain*** getOrigin() {
     if (image_data) {
       return image_data->getOrigin();
     } else {
       // should never get here
       std::cerr << "Origin requested while matrix3d undefined!";
       exit(1);
       // Uhm, throw exception instead?
     }
  }

  image3d() {
    init(); 
    if(bv::debug)
      std::cout << "Creating a new image3d object" << std::endl;
  }

  // make this unpublic?? friend???
  codomain** operator[](dimIndexType z) {
    return (*image_data)[z];
  }

  void init_to_zeros() {
    if (image_data) {
      image_data->init_to_zeros();
    }
    if(bv::debug)
      std::cout << "3d matrix initialized to all zeros" << std::endl;
  }

  void set_dims_no_init(dimSizeType N_z, dimSizeType N_y, dimSizeType N_x) { 
    try {
      image_data = new matrix3d<codomain>(N_z, N_y, N_x);
    }
    catch(...) {
      if (image_data) delete image_data;
      std::cerr << "in image3d::set_dims_no_init(" << N_z << ", " << N_y << ", " 
		<< N_x << "): allocation of a new matrix3d failed." << std::endl;
      throw;
    }
    this->N_z = N_z;
    this->N_y = N_y;
    this->N_x = N_x;
    // need these casts from dimSizeType to dimIndexType
    min_z = - static_cast<dimIndexType>(N_z/2); 
    min_y = - static_cast<dimIndexType>(N_y/2);
    min_x = - static_cast<dimIndexType>(N_x/2);
    max_z = static_cast<dimIndexType>(N_z-1)/2;
    max_y = static_cast<dimIndexType>(N_y-1)/2;
    max_x = static_cast<dimIndexType>(N_x-1)/2;
    matrix_set = true;
  } // --image3d::set_dims_no_init()
  
  void set_dims(dimSizeType N_z, dimSizeType N_y, dimSizeType N_x) { 
    set_dims_no_init(N_z, N_y, N_x);
    init_to_zeros();
  } // --image3d::set_dims()
  
  // pre: image_data must be already initialized
  void copy_from_buffer(codomain *inbuf) {
    if(image_data)
      image_data->copy_from_buffer(inbuf);
    else
      throw std::invalid_argument("image3d::copy_from_buffer requires imagedata to be initialized!");
  } // --image3d::copy_from_buffer()

 
  /** load data given a DIGFile reference                    *
   * WARNING! not to be called without checking sizes first! */
  int load(DIGFile &In) {
    // pre: already existing/allocated image_data matrix3d object
    unsigned int BufferSize = 0;
    In.GetArrayBufferSize(&BufferSize);
    unsigned int Channels = 0;
    In.GetChannels(&Channels);
    if (Channels != 1)
      throw std::runtime_error("Number of channels must be 1 in image3d::load()");
    unsigned int suspected_type_size = BufferSize / N_z / N_y / N_x / Channels;
    if (suspected_type_size != sizeof(codomain))
      throw std::runtime_error("Type size based on buffersize does not match codomain size in image3d::load()");

    if(bv::debug)
      std::cout << "About to load data into matrix from a DIGFile!" << std::endl;
    codomain* bufp = (&((*image_data)[min_z][min_y][min_x]));
	  
		printf(" size[%d][%d][%d]=%d \n",N_x,N_y,N_z,N_x*N_y*N_z);
 
    if(0 != In.GetArrayData(static_cast<void *>(bufp) ) )
      throw std::runtime_error("Error calling InputDIGFile.GetArrayData() from image3d::load()");

    return 0;
  } // --image3d::load()

  /** Read Bruno-type bcc (blob) images from buffer                           *
   *   pre: matrix and buffer MUST have appropriate dimensions wrt each other *
   *   pre: (would be nice if) matrix is initialized to all 0s                */
  void bccread(codomain *buf) {
    codomain*** im = getOrigin();   // image's origin
    // the "boolean" offset affects INITIAL y- and x-indices at each z-layer.
    // it is defined so that z, y, and x have the same "parity" (are equivalent mod 2)
    // might not be the fastest way, but does allows the buffer to be read sequentially
    int y_offset = (min_z-min_y)%2; // min_y + y_offset = min_z (mod 2) so that y = z (mod 2)
    int x_offset = (min_z-min_x)%2; // min_x + x_offset = min_z (mod 2) so that x = z (mod 2)
    for(dimIndexType z = min_z; z <= max_z; z++) {
      for(dimIndexType y = min_y + y_offset; y <= max_y; y+=2) {
				for(dimIndexType x = min_x + x_offset; x <= max_x; x+=2) {
	 			 im[z][y][x] = *buf++;
				}
      }
      x_offset = 1 - x_offset; // NOT A BUG--don't want to change x-offset at each row.
      y_offset = 1 - y_offset;
    }
  } // --image3d::bccread()

  /** Samples the blob based image (1st param) into the voronoi image *this   *
   * Voronoi image *this assumed to be on SC grid.                            *
   * Blob radius is passed in in terms of the DIG grid size                   *
   *    (equal to *HALF* of a side of the sparse elementar cube for bcc!)     *
   *  Pre: codomain (common to *this and *blobIm) is floating point.          *
   *  Pre: Blob image is on DIGGrid BCC1 [look only where "i=j=k (mod 2)"]    *
   * Automagically makes voxel-image 'include' the support of the blob-image. *
   * Note that for blobs on SC grid, this code would require minor change--   *
   *     just don't make it skip based on parity!                             */
  void convert_from_blob(image3d<codomain>* blobIm, double radius = 0, bool normalize = true, double blobnorm = 0.5682131592218539) {
    if(bv::debug)
      std::cout << "executing image3d::convert_from_blob()" << std::endl;
    if(radius == 0) {
      radius = static_cast<codomain>(2*sqrt(2));
    }
    if(radius <= 0) {
      throw std::invalid_argument("image3d::convert_from_blob requires positive blob radius.");
    }
    if(bv::debug) {
      std::cout << "Blob radius is " << radius << " in blob grid units." << std::endl;
    }
    // the players:
    const int m=2; // parameters to the modified Kaiser-Bessel window function
    const double alpha=10.4;
    const dimIndexType& Vmin_z = min_z;
    const dimIndexType& Vmin_y = min_y;
    const dimIndexType& Vmin_x = min_x;
    const dimIndexType  Bmin_z = blobIm->min_z;
    const dimIndexType  Bmin_y = blobIm->min_y;
    const dimIndexType  Bmin_x = blobIm->min_x;
    const dimIndexType& Vmax_z = max_z;
    const dimIndexType& Vmax_y = max_y;
    const dimIndexType& Vmax_x = max_x;
    if((Vmax_z<1) || (Vmax_y<1) || (Vmax_x<1))
      throw std::invalid_argument("In image3d::convert_from_blob(), this voronoi image must contain 3x3x3 'cubes'.");
    const dimIndexType  Bmax_z = blobIm->max_z;
    const dimIndexType  Bmax_y = blobIm->max_y;
    const dimIndexType  Bmax_x = blobIm->max_x;
    codomain***         voxOri = getOrigin();
    codomain***         blobOri = blobIm->getOrigin();
    dimIndexType        Bz, By, Bx, Vz, Vy, Vx, Vzlo, Vzhi, Vylo, Vyhi, Vxlo, Vxhi;
    double Vzdif, Vydif, Vxdif;
    double Vzdifsq, Vydifsq, Vxdifsq;
    double Vdifsqsum, radial;
    codomain blob_coefficient; // float or double expected
    codomain kaisarr[10001]; // set up kaisarr on-the-fly: extra pos for sanity
    std::cout << "Now creating kaiser array" << std::endl;
    for(int i=0; i<=10000; i++) {
      double radial= (radius*i)/10000;
      kaisarr[i]=kaiser_value(radial, radius, alpha, m);
    }
    std::cout << "Done reating kaiser array" << std::endl;
    // although the dimensions of the two images are fixed,
    // it is up to this function to implement spatial correspondence.
    // treat the voxel image (based on the blob radius) as just-fitting the blob image.    
    const bool dofudge = false; // don't use fudge at all if false
    const double fudge = 0.0; // might want to increase (or decrease) this
    // having it 0.0 makes it a requirement for voxOri[1][1][1] to exist.
    // make it 1.5 and you'll get a segmentation fault.    
    // "intuition": Vmax_i+fudge increases => same |Bi| corresponds to greater |Vi|
    // variable names abbreviated from "Blob Grid [coordinate] to Voxel Grid [coordinate]" etc.
    const double Z_bg2vg = dofudge ? ((Vmax_z + fudge) / (radius - Bmin_z)) : (Vmax_z/(radius - Bmin_z));
    const double Z_vg2bg = dofudge ? ((radius - Bmin_z) / (Vmax_z + fudge)) : ((radius - Bmin_z)/Vmax_z);
    const double Y_bg2vg = dofudge ? ((Vmax_y + fudge) / (radius - Bmin_y)) : (Vmax_y/(radius - Bmin_y));
    const double Y_vg2bg = dofudge ? ((radius - Bmin_y) / (Vmax_y + fudge)) : ((radius-Bmin_y)/Vmax_y);
    const double X_bg2vg = dofudge ? ((Vmax_x + fudge) / (radius - Bmin_x)) : (Vmax_x/(radius-Bmin_x));
    const double X_vg2bg = dofudge ? ((radius - Bmin_x) / (Vmax_x + fudge)) : ((radius-Bmin_x)/Vmax_x);
    // meaning: blobOri[Bz][By][Bx] <=>  voxOri[Bz*Z_bg2vg][By*Y_bg2vg][Bx*X_bg2vg]
    //           voxOri[Vz][Vy][Vx] <=> blobOri[Vz*Z_vg2bg][Vy*Y_vg2bg][Vx*X_vg2bg]
    
    // sanity check: (should never happen, even as an exception)
    if (! ( ( (Bmin_z-radius)*Z_bg2vg) >= Vmin_z) &&
	  ( ( (Bmin_y-radius)*Y_bg2vg) >= Vmin_y) &&
	  ( ( (Bmin_x-radius)*X_bg2vg) >= Vmin_x) &&
	  ( ( (Bmax_z+radius)*Z_bg2vg) <= Vmax_z) &&
	  ( ( (Bmax_y+radius)*Y_bg2vg) <= Vmax_y) &&
	  ( ( (Bmax_x+radius)*X_bg2vg) <= Vmax_x) ) {
      std::cerr << "Something in image3d::convert_from_blob very wrong." << std::endl;
      std::cerr << "Returning this with all zeros." << std::endl;
      return;
    }
    // Bruno's code treats SC as really cubic--I allow for stretched cubic grid

    // iterate through all blob centers on BCC1 grid
    int By_offset = (Bmin_z-Bmin_y)%2; // min_y + y_offset = min_z (mod 2) so that y = z (mod 2)
    int Bx_offset = (Bmin_z-Bmin_x)%2; // min_x + x_offset = min_z (mod 2) so that x = z (mod 2)
    std::cout << "z goes from " << Bmin_z << " to " << Bmax_z << std::endl;
    for(Bz = Bmin_z; Bz <= Bmax_z; Bz++) {
      //      std::c << "Bz is: ";
      std::cerr << " " << Bz << " "; // stdout doesn't work as fast, for some reason
      std::cout << " " << Bz << " ";
      // Vzlo and Vzhi are inclusive bounds
      Vzlo = static_cast<dimIndexType> ( ceil ((Bz - radius) * Z_bg2vg) );
      Vzhi = static_cast<dimIndexType> ( floor((Bz + radius) * Z_bg2vg) );
      //      std::cout << "Vzlo: " << Vzlo << " ; Vzhi: " << Vzhi << std::endl;
      // Vzlo > Vzhi implies that the support of this blob does not contain any voxel-sc centers
      if(Vzlo > Vzhi) continue;
      for(By = Bmin_y + By_offset; By <= Bmax_y; By+=2) {
	Vylo = static_cast<dimIndexType> ( ceil ((By - radius) * Y_bg2vg) );
	Vyhi = static_cast<dimIndexType> ( floor((By + radius) * Y_bg2vg) );
	if(Vylo > Vyhi) continue;
	for(Bx = Bmin_x + Bx_offset; Bx <= Bmax_x; Bx+=2) {
	  blob_coefficient = blobOri[Bz][By][Bx];
	  ///	  std::cout << "bc@["<<Bz<<"]["<<By<<"]["<<Bx<<"] = "<<blob_coefficient<<std::endl;
	  if(blob_coefficient == 0) continue;
	  Vxlo = static_cast<dimIndexType> ( ceil ((Bx - radius) * X_bg2vg) );
	  Vxhi = static_cast<dimIndexType> ( floor((Bx + radius) * X_bg2vg) );
	  if(Vxlo > Vxhi) continue;
	  // just calculate value of blob image AT sc grid point (no sophisticated averaging)

	  // have this option of normalizing
	  if(normalize) blob_coefficient *= blobnorm;

	  for(Vz = Vzlo; Vz <= Vzhi; Vz++) {

	    Vzdif = (Vz * Z_vg2bg) - Bz; // 0 plus or minus distance to Bz in blob grid units
	    Vzdifsq = Vzdif * Vzdif; // square distance from Vz to Bz in blob grid units
	    for(Vy = Vylo; Vy <= Vyhi; Vy++) {
	      Vydif = (Vy * Y_vg2bg) - By;
	      Vydifsq = Vydif * Vydif;
	      for(Vx = Vxlo; Vx <= Vxhi; Vx++) {
            Vxdif = (Vx * X_vg2bg) - Bx;
            Vxdifsq = Vxdif * Vxdif;
            Vdifsqsum = Vzdifsq + Vydifsq + Vxdifsq;
            radial = sqrt(Vdifsqsum); // radial distance of voxel center from blob center
            if (radial < radius) {
                voxOri[Vz][Vy][Vx] += blob_coefficient * kaisarr[static_cast<int>(floor(radial*10000/radius))];
            }
	      } // --Vx's in the x-support of this blob (voxel)
	    } // --Vy's in the y-support of this blob
	  } // --Vz's in the z-support of this blob
	} // --Bx's for this By (blob)
      } // --By's for this Bz
      Bx_offset = 1 - Bx_offset; // NOT A BUG--don't want to change x-offset at each y.
      By_offset = 1 - By_offset;
    } // --Bz's
    std::cout << std::endl;
  } // --image3d::convert_from_blob()
  
  /** converts this floating-point image3d into an uchar image3d with N bins [0..N-1] based on param
   *  unsigned char is 1 byte == 8 bits on our machines: can hold 0..255 
   *  @exception invalid_argument, new-related error
   */
  image3d<unsigned char>* float_to_uchar(unsigned int N) {
    if((N < 2) || (N > 256)) {
      throw std::invalid_argument("Argument to float_to_uchar MUST be between 2 and 256, inclusive.");
    }
    image3d<unsigned char>* newImgP;
    try {
      newImgP = new image3d<unsigned char>;
      newImgP->set_dims_no_init(N_z, N_y, N_x);
    } catch(...) {
      std::cerr << "Exception occured in uchar image3d initialization image3d::float_to_uchar()" << std::endl;
      if(newImgP) delete newImgP;
      throw;
    }
    image3d<unsigned char>& newIm = *newImgP;
    if(!minmax_set) calcMinMax();
    //  now minOfImage, maxOfImage are set
    codomain span = maxOfImage - minOfImage; // assume well-behaved-enough image so that this is possible...
    codomain*** orImOr  = getOrigin();   // original image's origin [a bit faster addressing]
    ///    unsigned char*** newImOr = newImgP->getOrigin(); // new image's origin
    // a bit slow but linear up to division so it really doesn't matter much
    dimIndexType z,y,x;
    unsigned int M = N-1;
    unsigned char imax = static_cast<unsigned char> (M);
    codomain gray;
    for(z=min_z; z<=max_z; z++) {
      for(y=min_y; y<=max_y; y++) {
				for(x=min_x; x<=max_x; x++) {
	  			gray = orImOr[z][y][x];
	  			if(gray >= maxOfImage) { // shouldn't ever be greater, but...
	    			newIm[z][y][x] = imax;
	    			///	 std::cerr << imax << " ";
	  			} else {
	    			newIm[z][y][x] = static_cast<unsigned char> ((gray - minOfImage) / span * N);
	    			///	 std::cerr << (static_cast<unsigned int>(newIm[z][y][x])) << std::endl;
	  			}
 				}
      }
    }
    return newImgP;
  } // image3d::float_to_uchar()
		

  image3d<long long>* codomain_to_long_Nbins(unsigned int N) {
    image3d<long long>* newImgP;
    try {
      newImgP = new image3d<long long>;
      newImgP->set_dims_no_init(N_z, N_y, N_x);
    } catch(...) {
      std::cerr << "Exception occured in uchar image3d initialization image3d::float_to_uchar()" << std::endl;
      if(newImgP) delete newImgP;
      throw;
    }
    image3d<long long>& newIm = *newImgP;
    if(!minmax_set) calcMinMax();
    //  now minOfImage, maxOfImage are set
    codomain span =maxOfImage - minOfImage; // assume well-behaved-enough image so that this is possible...
    printf("ORIGINAL maxOfImage =%f minOfImage=%f  \n",maxOfImage,minOfImage);
    newIm.maxOfImage=N;
    newIm.minOfImage=0;
    printf("maxOfImage =%f minOfImage=%f span=%f \n",newIm.maxOfImage,newIm.minOfImage,span);
    codomain*** orImOr  = getOrigin();   // original image's origin [a bit faster addressing]
    ///    unsigned char*** newImOr = newImgP->getOrigin(); // new image's origin
    // a bit slow but linear up to division so it really doesn't matter much
    dimIndexType z,y,x;
    unsigned int M = N-1;
    long long imax = static_cast<long long> (M);    
    codomain gray;			
    for(z=min_z; z<=max_z; z++) {
        for(y=min_y; y<=max_y; y++) {
            for(x=min_x; x<=max_x; x++) {
	  			gray = orImOr[z][y][x];
	  			if(gray >= convert_f2ll(maxOfImage)) { // shouldn't ever be greater, but...
	    			newIm[z][y][x] =  imax;
	  			} else {
                    if (gray==minOfImage)
                        newIm[z][y][x] = 0;
                    else
                        newIm[z][y][x] = static_cast<long long> ((gray - minOfImage) / span * N);
	  			}
            }
        }
    }
    printf(":: Different levels BINS LONG LONG]):  %d  \n",N);
    return newImgP;
  } // image3d::codomain_to_long_Nbins(		

long long convert_f2ll(codomain f){		
		long long ret = 0;	
		if (sizeof(f)<=4){
            ret = static_cast<long long>(f*10000000);
		}else{
			ret = static_cast<long long>(f*100000000000000);			
		}
		return ret;	
}
  
		/** converts this floating-point image3d into an long long image3d with N bins ( N is the number off
		differents gray scale. 8 byte == 64 bits on our machines: can hold 0..255 
   *  @exception invalid_argument, new-related error
   */
  image3d<long long>* codomain_to_long() {

    image3d<long long>* newImgP;
    try {
      newImgP = new image3d<long long>;
      newImgP->set_dims_no_init(N_z, N_y, N_x);
    } catch(...) {
      std::cerr << "Exception occured in uchar image3d initialization image3d::float_to_uchar()" << std::endl;
      if(newImgP) delete newImgP;
      throw;
    }
    image3d<long long>& newIm = *newImgP;
    if(!minmax_set) calcMinMax();

    //  now minOfImage, maxOfImage are set
    codomain span = abs(maxOfImage - minOfImage);
    codomain*** orImOr  = getOrigin();   // original image's origin [a bit faster addressing]
    long long*** ornewImgP  = newImgP->getOrigin();
    // a bit slow but linear up to division so it really doesn't matter much
    dimIndexType z,y,x;
    long long gray;
    //newIm.maxOfImage=convert_f2ll(maxOfImage);
    newImgP->maxOfImage = newIm.maxOfImage;
    //newIm.minOfImage=convert_f2ll(minOfImage);
    newImgP->minOfImage = newIm.minOfImage;
    printf("maxOfImage=%4.7f  minOfImage=%4.8f newImgP->maxOfImage=%lld newImgP->minOfImage=%lld\n",maxOfImage,minOfImage,newImgP->maxOfImage,newImgP->minOfImage );	// assume well-behaved-enough image so that this is possible...
    std::map<codomain, coordlist> teste;		
    std::map<long long, coordlist> testeLL;
    for(z=min_z; z<=max_z; z++) {
        for(y=min_y; y<=max_y; y++) {
            for(x=min_x; x<=max_x; x++){
	  			gray = convert_f2ll(orImOr[z][y][x]);
                //newIm[z][y][x] = gray;
                ornewImgP[z][y][x] = gray;
                coord coo = {z, y, x};
                teste[ orImOr[z][y][x] ].push_back(coo);
                testeLL[ newIm[z][y][x] ].push_back(coo);
            }
        }
    }
    printf(":: Different levels [LONG LONG])= %d  [FLOAT]=%d \n",testeLL.size(),teste.size());
    newImgP->calcMinMax();
    return newImgP;
  } // image3d::float_to_uchar()

  /* saves the image3d into DIGFile format                                            *
   * save takes a REFERENCE (not pointer) to a DIGFile                                *
   * pre: everything (other than copying of the image buffer) is already done to Out  *
   */
    int save(DIGFile &Out) {
    codomain* bufp;
    if(image_data)
      bufp = image_data->getFirst();
    else return -1;
    if(0 != Out.AppendArraySet(
           "PHAN", // ??
           "the only arrayset",
           "here goes array set 1 parameter string",
           "array set 1: this image is not necessarily a 'phantom'"
          )
       ){
			std::cerr << "Error appending an array set to the file" << std::endl;
            return 1;
    }
    if(0 != Out.AppendArray(
            1,
            "only array of the only array set",
            static_cast<void *>(bufp)
            )
       )
      {
        std::cerr << "Error appending the array to the file" << std::endl;
        std::cerr << "Quitting." << std::endl;
        return 1;
      }
    return 0;
  } // --image3d::save()

  /** return pointer to a flip/rotate combination of the image            *
   * takes an encoding of axes permutation which maps +z,+y,+x to one of: *
   *         +z: 1, +y: 2, +x: 3 ; -z:-1,-y:-2,-x:-3                      *
   *   each.  enforces that the input be a legitimate permutation         *
   *   (there are 6*4*2 = 48 = (3! * 2^3) legitimate permutations)        *
   *   (( including the identity (1, 2, 3) ))                             *
   * example for semantics: if zto is -2, then old z+ will map to new y-  *
   */
  image3d<codomain>* fliprot(short int zto, short int yto, short int xto) {
    // note that this would not have worked too well in-place anyway
    // especially in the case of rotation for non-cubic image domain.
    // i find it easier to return a pointer to a new instance, caller must delete it.
    short int ztdir = (zto>=0) ? 1 : -1;
    short int ytdir = (yto>=0) ? 1 : -1;
    short int xtdir = (xto>=0) ? 1 : -1;
    std::set<short int> verifier;
    //    vtdir == signum vto
    //    vto * vtodir == absolute value of vto
    verifier.insert(zto*ztdir);
    verifier.insert(yto*ytdir);
    verifier.insert(xto*xtdir);
    if ( (verifier.find(1) == verifier.end()) ||
	 (verifier.find(2) == verifier.end()) ||
	 (verifier.find(3) == verifier.end()) ) {
      std::cerr << "Invalid input to fliprot, returning null pointer." << std::endl;
      return 0;
    } // --if erroneous input to fliprot
    dimIndexType z, y, x;         // indices for old image
    dimIndexType zi, yi, xi;      // incides for new image
    dimSizeType new_N_z = 0;      // new sizes
    dimSizeType new_N_y = 0;
    dimSizeType new_N_x = 0;
    // pointers for what indices map each axis
    dimIndexType *ztop;
    dimIndexType *ytop;
    dimIndexType *xtop;
    // decided to go with simple pointers
    switch(zto*ztdir) {
    case 1:  new_N_z = N_z;  ztop = &zi; break;
    case 2:  new_N_y = N_z;  ztop = &yi; break;
    default: new_N_x = N_z;  ztop = &xi; // case 3
    }
    
    switch(yto*ytdir) {
    case 1:  new_N_z = N_y;  ytop = &zi; break;
    case 2:  new_N_y = N_y;  ytop = &yi; break;
    default: new_N_x = N_y;  ytop = &xi; // case 3
    }

    switch(xto*xtdir) {
    case 1:  new_N_z = N_x;  xtop = &zi; break;
    case 2:  new_N_y = N_x;  xtop = &yi; break;
    default: new_N_x = N_x;  xtop = &xi; // case 3
    }

    // set references for what z, y, x correspond to
    dimIndexType &ztoir = *ztop; // "the z-to index reference" etc.
    dimIndexType &ytoir = *ytop; 
    dimIndexType &xtoir = *xtop;
    // the following are used in the initialization for clarity
    dimIndexType begin_ext_zto = (ztdir>0) ? min_z : max_z; // beginning extremum from z
    dimIndexType begin_ext_yto = (ytdir>0) ? min_y : max_y; // beginning extremum from y
    dimIndexType begin_ext_xto = (xtdir>0) ? min_x : max_x; // beginning extremum from x
    image3d<codomain> *newImgP;
    try {
      newImgP = new image3d<codomain>;
      newImgP->set_dims_no_init(new_N_z, new_N_y, new_N_x);
    } catch(...) {
      std::cerr<< "Exception occured in image3d::fliprot()!" << std::endl;
      if(newImgP) delete newImgP;
      throw;
    }
    codomain*** orImOr  = getOrigin();   // original image's origin
    codomain*** newImOr = newImgP->getOrigin(); // new image's origin
    for(z=min_z, ztoir = begin_ext_zto; z<=max_z; z++, ztoir+=ztdir) {
      for(y=min_y, ytoir = begin_ext_yto; y<=max_y; y++, ytoir+=ytdir) {
	for(x=min_x, xtoir = begin_ext_xto; x<=max_x; x++, xtoir+=xtdir) {
	  newImOr[zi][yi][xi] = orImOr[z][y][x];
 	}
      }
    }
    return newImgP;
  } // --image3d::fliprot()

  /** create a stretched version of this image by (integral) multiplier,   *
   * return pointer to stretched image                                     */
  /* (idea: doubles as copy constructor if you pass in 1 as multiplier)    */
  image3d<codomain>* operator*(short int multiplier) {
    if (multiplier == 0) {
      std::cout << "Zero multiplier passed to operator* -- returning null value"
		<< std::endl;
      return 0;
    }
    dimSizeType new_N_z = static_cast<dimSizeType> (N_z * multiplier);
    dimSizeType new_N_y = static_cast<dimSizeType> (N_y * multiplier);
    dimSizeType new_N_x = static_cast<dimSizeType> (N_x * multiplier);
    if(bv::debug) {
      std::cout << "Creating stretched image with dims ";
      std::cout << new_N_z << " x " << new_N_y << " x " << new_N_x << std::endl;  
    }
    if (0 == matrix_set) {
      std::cerr << "WARNING: trying to multiply image3d without a matrix3d!" << std::endl;
      std::cerr << "Returning null!" << std::endl;
      return 0;
    }
    image3d<codomain> *newImgP;
    newImgP = new image3d<codomain>;
    newImgP->set_dims_no_init(new_N_z, new_N_y, new_N_x);
    //    newImgP->set_dims(new_N_z, new_N_y, new_N_x);
    codomain*** orImOr  = getOrigin();   // original image's origin
    codomain*** newImOr = newImgP->getOrigin(); // new image's origin
    dimIndexType z, y, x, new_z, new_y, new_x;
    short int i, j, k;
    codomain val;
    for(z = min_z, new_z = newImgP->min_z; z<=max_z; z++, new_z += multiplier) {
      for(y = min_y, new_y = newImgP->min_y; y<=max_y; y++, new_y += multiplier) {
	for(x = min_x, new_x = newImgP->min_x; x<=max_x; x++, new_x += multiplier) {
	  val = orImOr[z][y][x];
  	  for(i = 0; i<multiplier; i++) {
  	    for(j = 0; j<multiplier; j++) {
  	      for(k = 0; k<multiplier; k++) {
		newImOr[new_z+i][new_y+j][new_x+k] = val;
  	      }
  	    }
  	  }
	}
      }
    }
    return newImgP;
  } // --image3d::operator*()

  /** fill the already-initialized bool image3d with 'random' values   *
   *  won't complain if codomain is not bool--will fill with 0s and 1s *
   */
  void randomize(double ranthresh = 0.8) {
    // higher ranthresh means more 1s.
    // ranthresh should be between 0 and 1 for of random image
    srand ( time(0) ); // so that we get different random numbers every 'time'
    dimIndexType z, y, x;
    codomain*** im = getOrigin();
    for(z=min_z; z<=max_z; z++) {
      for(y=min_y; y<=max_y; y++) {
	for(x=min_x; x<=max_x; x++) {
	  // generate a "random" probability
	  double prob = (1.0*rand())/(RAND_MAX);
	  //	  std::cout << prob << std::endl;
	  if (ranthresh > prob) {
	    im[z][y][x] = 1;
	  } else {
	    im[z][y][x] = 0;
	  }
	} // --for x
      } // --for y
    } // --for z
  } // --image3d::randomize()

  /** calculates minimum and maximum values of image *
   * precondition: image is already initialized!     *
   */
  void calcMinMax() {
    if (matrix_set)
      image_data->calcMinMax(&minOfImage, &maxOfImage);
    minmax_set = true;
  } // --image3d::calcMinMax()
  
  void reverse() { 
    if(matrix_set) {
      image_data->reverse();
      std::cout << "Image successfully 'reversed'... we think." << std::endl;
    }
    else std::cerr << "3dmatrix not even allocated to be reversed yet!";    
  } // --image3d::reverse()
  

  bettiType betti_zero() {
    // b0 is a data member reserved for this.
    if(bv::debug)
      std::cout << "Running image3d::betti_zero()..." << std::endl;
    // b0 = num_connected_components(6, false, false); //  DON'T prune it
    b0 = num_connected_components(6, false, true); //  prune it
    return b0;
  } // --image3d::betti_zero()

  bettiType betti_two() {
    // b2 is a data member reserved for this.
    if(bv::debug)
      std::cout << "Running image3d::betti_two()..." << std::endl;
    if((N_z>2) && (N_y>2) && (N_x>2)) { // are 0-components possible?
      multilabel_matrix->reverse_discrete(); // make 0s the new 1s
      b2 = num_connected_components(26, true, false); // don't prune it
      multilabel_matrix->reverse_discrete(); // 1s 1s again to prepare for euler
    } else { // all voxels are on some face, no 'reversing' necessary
      b2=0;
    }
    return b2;
  } // --image3d::betti_two()

  /* betti_one() MUST be called after betti_zero() and betti_two() */
  /* we could use safeguards to ensure that, among other things... */
  bettiType betti_one() {
    // b1 is a data member reserved for this.
    if(bv::debug)
      std::cout << "Running image3d::betti_one()..." << std::endl;
    // a bit inefficient: just so euler can have sensicle [sic]
    eulerType E = euler();
    // E gets the Euler number of the image (ones)
    // E and bi are guaranteed to be nonnegative
    if(bv::debug)
      std::cout << "in image3d::betti_one() E is " << E << std::endl;
    b1 = b0 + b2 - E;
    if(bv::debug)
      std::cout << "in image3d::betti_one() b1 is " << b1 << std::endl;
    return b1;
  } // --image3d::betti_zero()

  /** initialize, and if necessary, create, multilabel_matrix */
  void mm_init() {
    if(multilabel_matrix == 0) {
      if(bv::debug)
				std::cout << "About to allocate label matrix" << std::endl;
      try {
				multilabel_matrix = new matrix3d<labelType>(N_z, N_y, N_x);
      } catch(...) { // probably bad_alloc
				std::cerr << "in image3d::mm_init(): allocation of"
		  	<< " a new matrix3d failed." << std::endl;
				throw;
      }      
      if(bv::debug)
				std::cout << "Allocated label matrix successfully" << std::endl;
    }
		multilabel_matrix->init_to_zeros();
  } // --matrix3d::mm_init()

  /** matrix3d::thresholdize() creates 'binary' image from image_data *
   *  places the result into multilabel_matrix
   *  this choice is based on space economy: 
   *  we don't need to hold 3 matrices in memory.
   *  (one for the input, one for binary, one for labels)
   */
  void thresholdize() {
    // instance variable image3d::threshold is an "implicit parameter"
    if (multilabel_matrix == 0) {
      mm_init();
    }
    // it's fine if there's already stuff in multilabel_matrix 
    // from last iteration... soon it'll be wiped out
    codomain *x_ptr = image_data->getFirst();
    codomain *x_bound = image_data->getPast();
    labelType *mm_ptr = multilabel_matrix->getFirst();
    while(x_ptr != x_bound) {
      *mm_ptr++ = (threshold < *x_ptr++) ? 1 : 0;
    }
    /// don't need    
    // sanity check
    if(mm_ptr != multilabel_matrix->getPast()) {
      std::cerr << "Absurd exit condition in thresholdize()." << std::endl;
    }
  } // --image3d::thresholdize()



  /** matrix3d::saveOpenDxFile()  *
   *  This function create a array file to be read by openDx
   */
	void saveOpenDxFile(){
		ofstream myfile("exampleVol.data");
  	if (myfile.is_open()){
			float val;
				
			codomain*** orImOr  = getOrigin();   // original image's origin [a bit faster addressing]
			// a bit slow but linear up to division so it really doesn't matter much
			dimIndexType z,y,x;
			for(z=min_z; z<=max_z; z++) {
				for(y=min_y; y<=max_y; y++) {
					for(x=min_x; x<=max_x; x++) {
						val = orImOr[z][y][x];
						if (val < 46.0){
							 val = 0.0f;	
						}else{
							//if ( (x>-10) & (x<5)& (y>-0) & (y<5)& (z>-10) & (z<5) ){
							//	val = 20.0f;	
								//printf("%f ",val);
							//}else{
								val = 10.0f;
								//printf("%f ",val);	
							//}								
						}		
						//printf("%f ",val);	 						
  					myfile  << fixed << setprecision(5) << val << " ";
					}
						myfile << endl;
				}
			}				

			//myfile  << fixed << setprecision(5) << val << endl;
			myfile.close();
  	}
  	else 
	    cout << "Unable to open file";
	
	}// --image3d::saveOpenDxFile()


  /** matrix3d::saveOpenDxFileFHTS()  *
   *  This function create a array file to be read by openDx with the FHTS labels
   */
	void saveOpenDxFileFHTS(historytree * htd){
			
        //htd->myTravesa();
		int flag =0;
		ofstream myfile("exampleVolFHTS.data");
        if (myfile.is_open()){
			float val;
			codomain auxs;	
			codomain*** orImOr  = getOrigin();   // original image's origin [a bit faster addressing]
			// a bit slow but linear up to division so it really doesn't matter much
			dimIndexType z,y,x;
			int cont1,cont2,cont3;
			cont1=cont2=cont3=0;
			for(z=min_z; z<=max_z; z++) {
				for(y=min_y; y<=max_y; y++) {
					for(x=min_x; x<=max_x; x++) {
						val = orImOr[z][y][x];	
						orImOr[z][y][x] = 0.0f;	
						//~ if (val < 34.0){
							 //~ cont1++;	
							 //~ orImOr[z][y][x] = 0.0f;	
						//~ }else{	
							//~ cont2++;	
							//~ orImOr[z][y][x]= 10.00f;	
						//~ }		
					}
				}
			}				
	
			std::queue<historynode*> s;
		  s.push(htd->hroot);	 
			auxs=20.00000;
			while(!(s.empty())) {
				historynode* np_; // node pointer
				historynode* chp_; // child pointer
				np_ = s.front();
	
				s.pop();
				//long long lglg ;
				//if (np_->gray>lglg){	
				std::vector<coord>::iterator theIterator;
			  if (flag > 2 && flag <6 ){						
    			for (theIterator = np_->listCoord.begin(); theIterator != np_->listCoord.end();
	        	theIterator++){
						cont3++;		
						coord mycoord = (*theIterator);
						if ( (mycoord.zpos>50) || (mycoord.ypos >50) ||  (mycoord.xpos >50) || (mycoord.zpos<-50) || (mycoord.ypos<-50) ||  (mycoord.xpos <-50) ){
							printf("(%d %d %d)[%d] ",mycoord.zpos,mycoord.ypos,mycoord.xpos,np_->gray);
						}else{		
							//if (orImOr[mycoord.zpos][mycoord.ypos][mycoord.xpos]!=0.0){	
							printf("*",mycoord.zpos,mycoord.ypos,mycoord.xpos,np_->gray);		
							orImOr[mycoord.zpos][mycoord.ypos][mycoord.xpos] = auxs;			
							//}
						}
					}
    		}
					
            if (np_->hfirstchild!=0){
					chp_ = (np_->hfirstchild);			
					s.push(chp_);		
					while(chp_->hnext != 0) {				
						historynode* aux = chp_->hnext;				
						if ( aux != np_){	
							s.push(aux);	
						}
						chp_ = aux;	
					}									
				}else{
					printf("Acabou .... \n ");	
					flag++;
				}
			}	
		
			//printf(" cont1=[%d] -- cont2=[%d] -- cont3=[%d] \n",cont1,cont2,cont3);	
			
			for(z=min_z; z<=max_z; z++) {
				for(y=min_y; y<=max_y; y++) {
					for(x=min_x; x<=max_x; x++) {
						val = orImOr[z][y][x];					
  					myfile  << fixed << setprecision(5) << val << " ";	
					}
						myfile << endl;
				}
			}						

			myfile.close();
  	}
  	else 
	    cout << "Unable to open file";
	
	}// --image3d::saveOpenDxFileFHTS()

  /** returns foreground components history tree pointer to totally break (pseudo-)encapsulation */
  historytree* get_htp() { 
    if(HTP==0) {
      std::cout << "Warning!  null HTP being returned from image3d::get_htp()" << std::endl;
    }
    return HTP;
  }

  /** returns background components history tree pointer to totally break (pseudo-)encapsulation */
  historytree_bgnd* get_htbp() { 
    if(HTBP==0) {
      std::cout << "Warning!  null HTBP being returned from image3d::get_htbp()" << std::endl;
    }
    return HTBP;
  }

  /** 2nd version of history tree construction
      DO NOT CALL IN CONJUNCTION WITH ANY OF THE PERSISTENCE STUFF
      (will explode the memory!)
      constructs THEN prunes
  */
  void construct_history_tree(unsigned minsize=0, unsigned minheight=0) {
    // just for 'b0' (foreground)
    if(HTP) delete HTP;
    HTP = new historytree();
    historytree& HT = *HTP;
    if(bv::debug) std::cout << "Executing image3d::construct_history_tree()" << std::endl;
    // Yung had suggested to use a list of lists, and sort the gray values.
    // This (the map) way is not noticeably slower.
    dimIndexType z, y, x;
    if(all_voxels.empty()) {     // if calling the first time
      //      std::cout << "Calling the first time" << std::endl;
      codomain*** im = getOrigin();
      for(z = min_z; z <= max_z; z++) {
				for(y = min_y; y <= max_y; y++) { 
					for(x = min_x; x <= max_x; x++) {
						coord coo = {z, y, x};
						all_voxels[ im[z][y][x] ].push_back(coo);
						// the outside set of braces means find value in the map
					}
				}
      }
      unset_matrix(); // because we can! (to save memory)
      im = 0;
    }
    // go through each gray value, starting with largest
    if(bv::debug) std::cout << "Traversing all_voxels" << std::endl;
    mm_init(); // initialize multilabel matrix (to all zeros)
    labelType*** mm = multilabel_matrix->getOrigin();
    //    if(multilabel_matrix->is_all_zero()) {
    //      std::cout << "multilabel matrix really is all zeros" << std::endl;
    //    }
    std::cout << "bv::debug: About to make new FGDSF" << std::endl;
    if(FGDSF) delete FGDSF; // save space
    FGDSF = new dsf_for_tree_gen();
    dsf_for_tree_gen& myDSF = (*FGDSF);
    labelType* neighbors[26]; // an array of codomain pointers for the 6 (potential) neighbors
    // Will(?) use the same one for 26-connectivity, that's why it has size 26.
    labelType common, temp;
    unsigned i; // "dummy" loop var used below
    //    unsigned prnumnod=0, curnumnod=0;
    for(typename list_of_lists::reverse_iterator gray_it = (all_voxels.rbegin());
	gray_it != (all_voxels.rend()); gray_it++) {
      // the non-zero voxels of mm correspond to 'set'
      // cast because unsigned chars are by default printed like characters not numbers
      if(bv::debug) {
	std::cout
	  << "Processing the " <<  ((gray_it->second).size()) << " voxel(s) at "
	  << "grayvalue " << (static_cast<unsigned int>(gray_it->first)) << std::endl;
      }
      for(typename coordlist::const_iterator coord_it = (gray_it->second.begin());
	  coord_it != (gray_it->second.end()); coord_it++) {
	z = (*coord_it).zpos;
	y = (*coord_it).ypos;
	x = (*coord_it).xpos;
	//	std::cerr << "accessed coordit stuff" << std::endl;
	// check for 6-connectivity here
	neighbors[0] = (z > min_z) ? (&(mm[z-1][y][x])) : 0;
	neighbors[1] = (z < max_z) ? (&(mm[z+1][y][x])) : 0;
	neighbors[2] = (y > min_y) ? (&(mm[z][y-1][x])) : 0;
	neighbors[3] = (y < max_y) ? (&(mm[z][y+1][x])) : 0;
	neighbors[4] = (x > min_x) ? (&(mm[z][y][x-1])) : 0;
	neighbors[5] = (x < max_x) ? (&(mm[z][y][x+1])) : 0;
	for(i = 0; i<6; i++) {
	  if((neighbors[i]) && (common = (*(neighbors[i]))))
	    // we DO want the value of the assignment
	    break;
	}

	// i remember reading somewhere NOT to alter the for loop variable...
	// hence the necessity of using "break" as opposed to an inner loop
	if(i<6) { // this means that at least one neighbor already a component
	  i++;
	  while(i<6) {
	    if(neighbors[i] && (temp = (*(neighbors[i]))))
	      common = myDSF.unite(common, temp);
	    i++;
	  }
	  mm[z][y][x] = common;
	  myDSF.increment_size(common);
	} else {
	  mm[z][y][x] = myDSF.makeSet(gray_it->first);
	}
      } // added all at this grayvalue
      //      myDSF.history_update(HT, minsize, minheight); // passed by reference
      myDSF.history_update(HT, gray_it->first); // HT passed by reference
      // at every new update, a new level is created and at least one entry is guaranteed.
      // at least one voxel.  if new, new node.  if joins something, updates that, which now gets new node.
      //      std::cout << "Processed gray " << (gray_it->first) << std::endl;
      if(bv::debug) {
				std::cout << "number of labels in DSF: " << myDSF.get_num_labels() << std::endl;
				std::cout << "number of current components: " << myDSF.get_num_sets() << std::endl;
				std::cout << "number of nodes in history tree: " << HT.numnodes << std::endl;
      }
    } // added all voxels
    if(bv::debug) std::cout << "Done traversing all_voxels for foreground " <<
		"component tree construction (the map way)" << std::endl;    
    if(bv::debug) std::cout << "setting HT (fgnd) root now..." << std::endl;
    myDSF.set_HT_root(HT);
    if(bv::debug) std::cout << "HT (fgnd) root set." << std::endl;

    if(bv::debug) {
      std::cout << "Displaying entire non-pruned tree." << std::endl;
      HT.traverse_nonpruned();
    }
    if((minsize | minheight) != 0) {
      //prune_history_tree(minsize, minheight);
    }
  } // --image3d::construct_history_tree()

  /** prune EXISTING history tree of foreground components, based on params */
  void prune_history_tree(unsigned minsize, unsigned minheight) {
    if(bv::debug)
      std::cout << "++image3d::prune_history_tree()" << std::endl;
    if(HTP==0) {
      std::cerr << "Error: no history tree to be pruned!" << std::endl;
      return;
    }
    // DO call prune in subsequent calls with minsize=0, minheight=0.
    HTP->prune(minsize, minheight);
    if(bv::debug) {
      std::cout << "Displaying pruned tree with minsize=" << minsize 
		<< ", minheight=" << minheight << std::endl;
      HTP->traverse_nonpruned();
    }
    if(bv::debug) {
      std::cout << "--image3d::prune_history_tree()" << std::endl;      
    }
  } // --image3d::prune_history_tree()


  /** 3rd version of history tree construction
      "ideally" called with no params
      @return ptr to new historytree corresponding to this binned image3d with life of its own (this image3d will never know about it again)
      it is the (ultimate) callee's responsibility to free the returned object
  */
  historytree* construct_history_tree_indep(unsigned minsize=0, unsigned minheight=0) {
    // just for 'b0' (foreground)
    //    if(HTP) delete HTP;
    //    HTP = new historytree();
    historytree* myHTP = new historytree(); // local to this fnc
    historytree& HT = *myHTP; // reference set as before.  smoove!
    std::cout << "Executing image3d::construct_history_tree_indep()" << std::endl;
    // Yung had suggested to use a list of lists, and sort the gray values.
    // This (the map) way is not noticeably slower.
    dimIndexType z, y, x;
	if(all_voxels.empty()) {     // if calling the first time
		//std::cout << "Calling the first time" << std::endl;
		codomain*** im = getOrigin();
		for(z = min_z; z <= max_z; z++) {
			for(y = min_y; y <= max_y; y++) { 
				for(x = min_x; x <= max_x; x++) {
					coord coo = {z, y, x};
					all_voxels[ im[z][y][x] ].push_back(coo);
					// the outside set of braces means find value in the map
				}
			}
      	}
      	unset_matrix(); // because we can! (to save memory)
      	im = 0;
    }
    // go through each gray value, starting with largest
    if(bv::debug) std::cout << "Traversing all_voxels" << std::endl;
    mm_init(); // initialize multilabel matrix (to all zeros)
    labelType*** mm = multilabel_matrix->getOrigin();
	  
    dsf_for_tree_gen myDSF;  // why use new when don't have to?
    labelType* neighbors[26]; // an array of codomain pointers for the 6 (potential) neighbors
    // Will(?) use the same one for 26-connectivity, that's why it has size 26.
    labelType common, temp;
    unsigned i; // "dummy" loop var used below
    //    unsigned prnumnod=0, curnumnod=0;
    for(typename list_of_lists::reverse_iterator gray_it = (all_voxels.rbegin());
		gray_it != (all_voxels.rend()); gray_it++) {
      	// the non-zero voxels of mm correspond to 'set'
      	// cast because unsigned chars are by default printed like characters not numbers
      	if(bv::debug) {
			std::cout
	  		<< "#Processing the " <<  ((gray_it->second).size()) << " voxel(s) at "
	  		<< "#grayvalue " << (static_cast<unsigned int>(gray_it->first)) << std::endl;
      	}

      	for(typename coordlist::const_iterator coord_it = (gray_it->second.begin());
	  	  coord_it != (gray_it->second.end()); coord_it++) {
			z = (*coord_it).zpos;
			y = (*coord_it).ypos;
			x = (*coord_it).xpos;

			// check for 6-connectivity here
			neighbors[0] = (z > min_z) ? (&(mm[z-1][y][x])) : 0;
			neighbors[1] = (z < max_z) ? (&(mm[z+1][y][x])) : 0;
			neighbors[2] = (y > min_y) ? (&(mm[z][y-1][x])) : 0;
			neighbors[3] = (y < max_y) ? (&(mm[z][y+1][x])) : 0;
			neighbors[4] = (x > min_x) ? (&(mm[z][y][x-1])) : 0;
			neighbors[5] = (x < max_x) ? (&(mm[z][y][x+1])) : 0;
			for(i = 0; i<6; i++) {
			  if((neighbors[i]) && (common = (*(neighbors[i]))))
				// we DO want the value of the assignment
				break;
			}
			//	std::cerr << "accessed neighbors in mm" << std::endl;
			// i remember reading somewhere NOT to alter the for loop variable...
			// hence the necessity of using "break" as opposed to an inner loop
			if(i<6) { // this means that at least one neighbor already a component
	  			i++;
	  			while(i<6) {
	    			if(neighbors[i] && (temp = (*(neighbors[i]))))
	      				common = myDSF.unite(common, temp);
	    			i++;
	  			}
	  			mm[z][y][x] = common;
	  			myDSF.increment_size(common);
			} else {
	  			mm[z][y][x] = myDSF.makeSet(gray_it->first);
			}
		} // added all at this grayvalue
		//      myDSF.history_update(HT, minsize, minheight); // passed by reference
		myDSF.history_update(HT, gray_it->first); // HT passed by reference
	} // added all voxels

    if(bv::debug) std::cout << "Done traversing all_voxels for foreground " <<
		"component tree construction (the map way)" << std::endl;    
	
	if(bv::debug) std::cout << "setting HT (fgnd) root now..." << std::endl;
		myDSF.set_HT_root(HT);
	
	if(bv::debug) std::cout << "HT (fgnd) root set." << std::endl;

	if(bv::debug) {
		std::cout << "Displaying entire non-pruned tree." << std::endl;
		HT.traverse_nonpruned();
	}
	
	if((minsize | minheight) != 0) {
		prune_history_tree(minsize, minheight);
	}
	return myHTP;
  } // --image3d::construct_history_tree_indep()

		
		
 /** MATAV3
    4rd version of history tree construction
    This function used all gray value as threshold level.
    "ideally" called with no params
     @return ptr to new historytree corresponding to this binned image3d with life of its own (this image3d will never know about it again)
     it is the (ultimate) callee's responsibility to free the returned object
  */
  historytree* construct_history_tree_long_bins(unsigned minsize=0, unsigned minheight=0) {

    historytree* myHTP = new historytree(); // local to this fnc
    historytree& HT = *myHTP; // reference set as before.  smoove!
    std::cout << "Executing image3d::construct_history_tree_float_bins(...)" << std::endl;

    dimIndexType z, y, x;
    if(all_voxels.empty()) {     // if calling the first time
        //std::cout << "Calling the first time" << std::endl;
        codomain*** im = getOrigin();
        for(z = min_z; z <= max_z; z++) {
            for(y = min_y; y <= max_y; y++) {
                for(x = min_x; x <= max_x; x++) {
                    coord coo = {z, y, x};
                    all_voxels[ im[z][y][x] ].push_back(coo);
                    // the outside set of braces means find value in the map
                }
            }
        }
        unset_matrix(); // because we can! (to save memory)
        im = 0;
    }
			
			
    // go through each gray value, starting with largest
    if(bv::debug) std::cout << "Traversing all_voxels" << std::endl;
			
    printf(":: Number of beans (different levels):  %d  \n",all_voxels.size());
    myHTP->nBins = 	all_voxels.size();
    mm_init(); // initialize multilabel matrix (to all zeros)
    labelType*** mm = multilabel_matrix->getOrigin();
	  
    dsf_for_tree_gen myDSF;  // why use new when don't have to?
    labelType* neighbors[26]; // an array of codomain pointers for the 6 (potential) neighbors
    // Will(?) use the same one for 26-connectivity, that's why it has size 26.
    labelType common, temp;
    unsigned i; // "dummy" loop var used below

    int total_steps = all_voxels.size();
    float actual_steps= 0.0;
    int per =5;
    int contadorGERAL = 0;
    int at =0;
    for(typename list_of_lists::reverse_iterator gray_it = (all_voxels.rbegin());
					gray_it != (all_voxels.rend()); gray_it++) {
      // the non-zero voxels of mm correspond to 'set'
      // cast because unsigned chars are by default printed like characters not numbers
        if(bv::debug) {
                printf("#Processing the %d  voxel(s) at %lld \n ",gray_it->second.size(),gray_it->first);
        }
        long long  au= 			gray_it->first;
        actual_steps+=1.0;
        if ( per <  (int) (ceil((actual_steps/total_steps)*100.0))){
                printf(" #[%d %] \n",per);
                per = per + 5;
        }
        contadorGERAL+=gray_it->second.size();
        std::map<labelType,std::vector<coord> > labelCoodMap;
        labelCoodMap.clear();
        int contador =0;
        int tt=0;
        for(typename coordlist::const_iterator coord_it = (gray_it->second.begin());
	  	  		coord_it != (gray_it->second.end()); coord_it++) {
				
            z = (*coord_it).zpos;
            y = (*coord_it).ypos;
            x = (*coord_it).xpos;

            neighbors[0] = (z > min_z) ? (&(mm[z-1][y][x])) : 0;
            neighbors[1] = (z < max_z) ? (&(mm[z+1][y][x])) : 0;
            neighbors[2] = (y > min_y) ? (&(mm[z][y-1][x])) : 0;
            neighbors[3] = (y < max_y) ? (&(mm[z][y+1][x])) : 0;
            neighbors[4] = (x > min_x) ? (&(mm[z][y][x-1])) : 0;
            neighbors[5] = (x < max_x) ? (&(mm[z][y][x+1])) : 0;
            for(i = 0; i<6; i++) {
                if((neighbors[i]) && (common = (*(neighbors[i]))))
                // we DO want the value of the assignment
                break;
            }
            // i remember reading somewhere NOT to alter the for loop variable...
            // hence the necessity of using "break" as opposed to an inner loop
            if(i<6) { // this means that at least one neighbor already a component
                i++;
                while(i<6) {
                    if(neighbors[i] && (temp = (*(neighbors[i])))){
                            common = myDSF.unite(common, temp);
                    }
                    i++;
                }
                mm[z][y][x] = common;
                unsigned int labelInc = myDSF.increment_size(common);
                labelCoodMap[labelInc].push_back((*coord_it));
            } else {
                mm[z][y][x] = myDSF.makeSet(gray_it->first);
                labelCoodMap[mm[z][y][x]].push_back((*coord_it));
            }
        } // added all at this grayvalue
        tt=0;
        std::map<labelType,std::vector<coord> >::iterator itt;
        for(itt=labelCoodMap.begin();itt != labelCoodMap.end(); itt++){
            tt+=(*itt).second.size();
        }
        if (tt!=gray_it->second.size())
            printf("************ Error!! tt!=gray_it->second.size() [%d != %d] => %lld  \n",tt,gray_it->second.size(),gray_it->first);
							
								
        myDSF.history_update(HT, gray_it->first,labelCoodMap,tt); // HT passed by reference
        std::map<labelType,std::vector<coord> >::iterator it;
        for(it=labelCoodMap.begin(); it != labelCoodMap.end(); it++){
            at+=(*it).second.size();
        }

    } // added all voxels
    if(bv::debug) std::cout << "Done traversing all_voxels for foreground " <<
        "component tree construction (the map way)" << std::endl;

			
    if(bv::debug) std::cout << "setting HT (fgnd) root now..." << std::endl;

    myDSF.set_HT_root(HT);

    if(bv::debug) std::cout << "HT (fgnd) root set." << std::endl;


    myHTP->nBins = all_voxels.size();
    myHTP->numnodes =HT.numnodes;
    historynode *t;
    std::cout << "number of labels in DSF: " << myDSF.get_num_labels() << std::endl;
    std::cout << "number of current components: " << myDSF.get_num_sets() << std::endl;
    std::cout << "number of nodes in history tree: " << HT.numnodes << std::endl;
    std::cout << "contadorGERAL (number of elements in all_voxels-coord): " << contadorGERAL << std::endl;
    std::cout << "memory used: " << contadorGERAL * (sizeof(historynode)) << " Bytes  ()" <<sizeof(historynode)<< std::endl;
    std::cout << "memory used: " << contadorGERAL * (sizeof(float)+(9*sizeof(t))+(3*sizeof(long long))+(6*sizeof(int))+sizeof(bool)) << " Bytes " << std::endl;
    printf("$ AT: numeber of coord that is suposed to add in M = %d\n",at);
    if(bv::debug) {
        std::cout << "Displaying entire non-pruned tree." << std::endl;
        HT.traverse_nonpruned();
    }
    //myHTP->myTravesa();

    return myHTP;
} // --image3d::construct_history_tree_indep()

int calculeFloatBins(){		
		float*** orImOr  = getOrigin(); 	
		dimIndexType z, y, x;		
		std::map<float, coordlist> listFloat;
		for(z = min_z; z <= max_z; z++) {
			for(y = min_y; y <= max_y; y++) { 
				for(x = min_x; x <= max_x; x++) {
					coord coo = {z, y, x};
					listFloat[ orImOr[z][y][x] ].push_back(coo);
					// the outside set of braces means find value in the map
				}
			}
		}
		int n_bins=listFloat.size();		
		printf("::calculeBins => Number of beans:  %d  \n",n_bins);		
		return 	n_bins;
	}

	int calculeDoubleBins(){		
		double*** orImOr  = getOrigin(); 	
		dimIndexType z, y, x;	
		std::map<double, coordlist> listDouble;	
		for(z = min_z; z <= max_z; z++) {
			for(y = min_y; y <= max_y; y++) { 
				for(x = min_x; x <= max_x; x++) {
					coord coo = {z, y, x};
					listDouble[ orImOr[z][y][x] ].push_back(coo);
					// the outside set of braces means find value in the map
				}
			}
		}
	
		int n_bins=listDouble.size();		
		printf("::calculeBins => Number of beans:  %d  \n ",n_bins);		
		return 	n_bins;
	}
	// go through each gray value, starting with largest


 
			
  /** version of history tree construction for background voxels
      might parametrize it: (mincomponentsize and pruneby)
      or have these be postprocessing things.
      DO NOT CALL IN CONJUNCTION WITH ANY OF THE PERSISTENCE STUFF
      (will explode the memory!)
      constructs THEN prunes
  */
  void construct_history_tree_bgnd(unsigned minsize=0, unsigned minheight=0) {
    // just for 'b2' (background)
    if(HTBP) delete HTBP;
    HTBP = new historytree_bgnd();
    historytree_bgnd& HT = *HTBP;
    if(bv::debug) std::cout << "Executing image3d::construct_history_tree_bgnd()" << std::endl;
    // Yung had suggested to use a list of lists, and sort the gray values.
    // This (the map) way is not noticeably slower.
    dimIndexType z, y, x;
    if(all_voxels.empty()) {     // if calling the first time
      std::cout << "in image3d::construct_history_tree_bgnd() should NEVER get here!" << std::endl;
      //      std::cout << "Calling the first time" << std::endl;
      codomain*** im = getOrigin();
      for(z = min_z; z <= max_z; z++) {
	for(y = min_y; y <= max_y; y++) { 
	  for(x = min_x; x <= max_x; x++) {
	    coord coo = {z, y, x};
	    all_voxels[ im[z][y][x] ].push_back(coo);
	    // the outside set of braces means find value in the map
	  }
	}
      }
      unset_matrix(); // because we can! (to save memory)
      im = 0;
    }
    // go through each gray value, starting with largest
    if(bv::debug) std::cout << "Traversing all_voxels" << std::endl;
    mm_init(); // initialize multilabel matrix (to all zeros)
    float *** mm = multilabel_matrix->getOrigin();
    //    if(multilabel_matrix->is_all_zero()) {
    //      std::cout << "multilabel matrix really is all zeros" << std::endl;
    //    }
    std::cout << "debug: About to make new BGDSF" << std::endl;
    if(BGDSF) delete BGDSF; // save space
    BGDSF = new dsf_for_tree_bgnd_gen();
    dsf_for_tree_bgnd_gen& myDSF = (*BGDSF);  // declaring local reference
    float * neighbors[26]; // an array of codomain pointers for the 6 (potential) neighbors
    // Will(?) use the same one for 26-connectivity, that's why it has size 26.
    float  common, temp;
    unsigned i; // "dummy" loop var used below
    // this time forward-iterate
    for(typename list_of_lists::iterator gray_it = (all_voxels.begin());
	gray_it != (all_voxels.end()); gray_it++) {
      // the non-zero voxels of mm correspond to 'set'
      // cast because unsigned chars are by default printed like characters not numbers
      if(bv::debug) {
	std::cout
	  << "Processing the " <<  ((gray_it->second).size()) << " voxel(s) at "
	  << "grayvalue " << (static_cast<unsigned int>(gray_it->first)) << std::endl;
      }
      for(typename coordlist::const_iterator coord_it = (gray_it->second.begin());
	  coord_it != (gray_it->second.end()); coord_it++) {
	z = (*coord_it).zpos;
	y = (*coord_it).ypos;
	x = (*coord_it).xpos;
	//	std::cerr << "accessed coordit stuff" << std::endl;
	// check for 26-connectivity here
	bool boundary = false;
	if(z > min_z) neighbors[0] = (&(mm[z-1][ y ][ x ]));
	else { neighbors[0] = 0; boundary = true; }
	if(z < max_z) neighbors[1] = (&(mm[z+1][ y ][ x ]));
	else { neighbors[1] = 0; boundary = true; }
	if(y > min_y) neighbors[2] = (&(mm[ z ][y-1][ x ]));
	else { neighbors[2] = 0; boundary = true; }
	if(y < max_y) neighbors[3] = (&(mm[ z ][y+1][ x ]));
	else { neighbors[3] = 0; boundary = true; }
	if(x > min_x) neighbors[4] = (&(mm[ z ][ y ][x-1]));
	else { neighbors[4] = 0; boundary = true; }
	if(x < max_x) neighbors[5] = (&(mm[ z ][ y ][x+1]));
	else { neighbors[5] = 0; boundary = true; }
	neighbors[ 6] = ((z > min_z) && (y > min_y))                 ?  (&(mm[z-1][y-1][ x ])) : 0;
	neighbors[ 7] = ((z > min_z) && (y < max_y))                 ?  (&(mm[z-1][y+1][ x ])) : 0;
	neighbors[ 8] = ((z < max_z) && (y > min_y))                 ?  (&(mm[z+1][y-1][ x ])) : 0;
	neighbors[ 9] = ((z < max_z) && (y < max_y))                 ?  (&(mm[z+1][y+1][ x ])) : 0;
	neighbors[10] = ((z > min_z) &&                (x > min_x))  ?  (&(mm[z-1][ y ][x-1])) : 0;
	neighbors[11] = ((z > min_z) &&                (x < max_x))  ?  (&(mm[z-1][ y ][x+1])) : 0;
	neighbors[12] = ((z < max_z) &&                (x > min_x))  ?  (&(mm[z+1][ y ][x-1])) : 0;
	neighbors[13] = ((z < max_z) &&                (x < max_x))  ?  (&(mm[z+1][ y ][x+1])) : 0;
	neighbors[14] =                ((y > min_y) && (x > min_x))  ?  (&(mm[ z ][y-1][x-1])) : 0;
	neighbors[15] =                ((y > min_y) && (x < max_x))  ?  (&(mm[ z ][y-1][x+1])) : 0;
	neighbors[16] =                ((y < max_y) && (x > min_x))  ?  (&(mm[ z ][y+1][x-1])) : 0;
	neighbors[17] =                ((y < max_y) && (x < max_x))  ?  (&(mm[ z ][y+1][x+1])) : 0;
	neighbors[18] = ((z > min_z) && (y > min_y) && (x > min_x))  ?  (&(mm[z-1][y-1][x-1])) : 0;
	neighbors[19] = ((z > min_z) && (y > min_y) && (x < max_x))  ?  (&(mm[z-1][y-1][x+1])) : 0;
	neighbors[20] = ((z > min_z) && (y < max_y) && (x > min_x))  ?  (&(mm[z-1][y+1][x-1])) : 0;
	neighbors[21] = ((z > min_z) && (y < max_y) && (x < max_x))  ?  (&(mm[z-1][y+1][x+1])) : 0;
	neighbors[22] = ((z < max_z) && (y > min_y) && (x > min_x))  ?  (&(mm[z+1][y-1][x-1])) : 0;
	neighbors[23] = ((z < max_z) && (y > min_y) && (x < max_x))  ?  (&(mm[z+1][y-1][x+1])) : 0;
	neighbors[24] = ((z < max_z) && (y < max_y) && (x > min_x))  ?  (&(mm[z+1][y+1][x-1])) : 0;
	neighbors[25] = ((z < max_z) && (y < max_y) && (x < max_x))  ?  (&(mm[z+1][y+1][x+1])) : 0;
	
  	// essentially same thing as we did for 6-connectivity
	for(i = 0; i<26; i++) {
	  if((neighbors[i]) && (common = (*(neighbors[i]))))
	    // we DO want the value of the assignment
	    break;
	}
	//	std::cerr << "accessed neighbors in mm" << std::endl;
	// i remember reading somewhere NOT to alter the for loop variable...
	// hence the necessity of using "break" as opposed to an inner loop
	if(i<26) { // this means that at least one neighbor already a component
	  i++;
	  while(i<26) {
	    if(neighbors[i] && (temp = (*(neighbors[i])))) {
	      common = myDSF.unite(common, temp);
	    }
	    i++;
	  }
  	  if(boundary) { // if we're at a boundary, can only destroy a component
  	    mm[z][y][x] = 1; // to make a point
	    myDSF.increment_size(1); // what the hell, just increment "size" of "boundary".
	    myDSF.unite(common, 1);
  	  } else {
  	    mm[z][y][x] = common;
	    myDSF.increment_size(common);
  	  }
  	} else {
  	  //  neighborless unbounded component of 0s: DO treat it like a component (already exists)
  	  if(boundary) {
	    mm[z][y][x] = 1;
	    myDSF.increment_size(1);
	  } else {
	    mm[z][y][x] = myDSF.makeSet(gray_it->first);
	  }
  	}
      } // -- done with this set of voxels with the same grayvalue
      myDSF.history_update(HT, gray_it->first); // HT passed by reference
      // at every new update, a new level is created and at least one entry is guaranteed.
      // at least one voxel.  if new, new node.  if joins something, updates that, which now gets new node.
      //      std::cout << "Processed gray " << (gray_it->first) << std::endl;
      if(bv::debug) {
	std::cout << "number of labels in DSF: " << myDSF.get_num_labels() << std::endl;
	std::cout << "number of current components: " << myDSF.get_num_sets() << std::endl;
	std::cout << "number of nodes in history tree: " << HT.numnodes << std::endl;
      }
    } // added all voxels
    if(bv::debug) std::cout << "Done traversing all_voxels for background " <<
		"component tree construction (the map way)" << std::endl;    
    if(bv::debug) std::cout << "setting HT (bgnd) root now..." << std::endl;
    myDSF.set_HT_root(HT);
    if(bv::debug) std::cout << "HT (bgnd) root set." << std::endl;
    //    std::cout << "number of labels in DSF: " << myDSF.get_num_labels() << std::endl;
    //    std::cout << "number of current components: " << myDSF.get_num_sets() << std::endl;
    //    std::cout << "number of nodes in history tree: " << HT.numnodes << std::endl;
    if(bv::debug) {
      std::cout << "Displaying entire non-pruned tree." << std::endl;
      HT.traverse_nonpruned();
    }
    if((minsize | minheight) != 0) {
      prune_history_tree_bgnd(minsize, minheight);
    }
  } // --image3d::construct_history_tree_bgnd()



  /** prune EXISTING history tree of background components, based on params */
  void prune_history_tree_bgnd(unsigned minsize, unsigned minheight) {
    if(bv::debug)
      std::cout << "++image3d::prune_history_tree_bgnd()" << std::endl;
    if(HTBP==0) {
      std::cerr << "Error: no historytree_bgnd to be pruned!" << std::endl;
      return;
    }
    // DO call prune in subsequent calls with minsize=0, minheight=0.
    HTBP->prune(minsize, minheight);
    if(bv::debug) {
      std::cout << "Displaying pruned tree_bgnd with minsize=" << minsize 
		<< ", minheight=" << minheight << std::endl;
      HTBP->traverse_nonpruned();
    }
    if(bv::debug) {
      std::cout << "--image3d::prune_history_tree_bgnd()" << std::endl;      
    }
  } // --image3d::prune_history_tree_bgnd()


  /** 3rd? version of bgnd history tree construction
      "ideally" called with no params
      @return ptr to a new historytree_bgnd corresponding to this binned image3d with life of its own (this image3d will never know about it again)
      it is the (ultimate) callee's responsibility to free the returned object
  */
  historytree_bgnd* construct_history_tree_bgnd_indep(unsigned minsize=0, unsigned minheight=0) {
    // just for 'b2' (background)
    //    if(HTBP) delete HTBP;
    //    HTBP = new historytree_bgnd();
    historytree_bgnd* myHTBP = new historytree_bgnd(); // local to this fnc!
    historytree_bgnd& HT = *myHTBP; // reference set as before.  smoove!
    if(bv::debug) std::cout << "Executing image3d::construct_history_tree_bgnd()" << std::endl;
    // Yung had suggested to use a list of lists, and sort the gray values.
    // This (the map) way is not noticeably slower.
    dimIndexType z, y, x;
    if(all_voxels.empty()) {     // if calling the first time
      std::cerr << "in image3d::construct_history_tree_bgnd() should NEVER get here!" << std::endl;
      //      std::cout << "Calling the first time" << std::endl;
      codomain*** im = getOrigin();
      for(z = min_z; z <= max_z; z++) {
	for(y = min_y; y <= max_y; y++) { 
	  for(x = min_x; x <= max_x; x++) {
	    coord coo = {z, y, x};
	    all_voxels[ im[z][y][x] ].push_back(coo);
	    // the outside set of braces means find value in the map
	  }
	}
      }
      unset_matrix(); // because we can! (to save memory)
      im = 0;
    }
    // go through each gray value, starting with largest
    if(bv::debug) std::cout << "Traversing all_voxels" << std::endl;
    mm_init(); // initialize multilabel matrix (to all zeros)
    float *** mm = multilabel_matrix->getOrigin();
    //    if(multilabel_matrix->is_all_zero()) {
    //      std::cout << "multilabel matrix really is all zeros" << std::endl;
    //    }
    //    std::cout << "debug: About to make a new BGDSF" << std::endl;
    //    if(BGDSF) delete BGDSF; // save space
    //    BGDSF = new dsf_for_tree_bgnd_gen();
    //    dsf_for_tree_bgnd_gen& myDSF = (*BGDSF);  // declaring local reference
    dsf_for_tree_bgnd_gen myDSF; // why use reference?  (won't use it again.)
    float * neighbors[26]; // an array of codomain pointers for the 6 (potential) neighbors
    // Will(?) use the same one for 26-connectivity, that's why it has size 26.
    float  common, temp;
    unsigned i; // "dummy" loop var used below
    // this time forward-iterate
    for(typename list_of_lists::iterator gray_it = (all_voxels.begin());
	gray_it != (all_voxels.end()); gray_it++) {
      // the non-zero voxels of mm correspond to 'set'
      // cast because unsigned chars are by default printed like characters not numbers
      if(bv::debug) {
	std::cout
	  << "Processing the " <<  ((gray_it->second).size()) << " voxel(s) at "
	  << "grayvalue " << (static_cast<unsigned int>(gray_it->first)) << std::endl;
      }
      for(typename coordlist::const_iterator coord_it = (gray_it->second.begin());
	  coord_it != (gray_it->second.end()); coord_it++) {
	z = (*coord_it).zpos;
	y = (*coord_it).ypos;
	x = (*coord_it).xpos;
	//	std::cerr << "accessed coordit stuff" << std::endl;
	// check for 26-connectivity here
	bool boundary = false;
	if(z > min_z) neighbors[0] = (&(mm[z-1][ y ][ x ]));
	else { neighbors[0] = 0; boundary = true; }
	if(z < max_z) neighbors[1] = (&(mm[z+1][ y ][ x ]));
	else { neighbors[1] = 0; boundary = true; }
	if(y > min_y) neighbors[2] = (&(mm[ z ][y-1][ x ]));
	else { neighbors[2] = 0; boundary = true; }
	if(y < max_y) neighbors[3] = (&(mm[ z ][y+1][ x ]));
	else { neighbors[3] = 0; boundary = true; }
	if(x > min_x) neighbors[4] = (&(mm[ z ][ y ][x-1]));
	else { neighbors[4] = 0; boundary = true; }
	if(x < max_x) neighbors[5] = (&(mm[ z ][ y ][x+1]));
	else { neighbors[5] = 0; boundary = true; }
	neighbors[ 6] = ((z > min_z) && (y > min_y))                 ?  (&(mm[z-1][y-1][ x ])) : 0;
	neighbors[ 7] = ((z > min_z) && (y < max_y))                 ?  (&(mm[z-1][y+1][ x ])) : 0;
	neighbors[ 8] = ((z < max_z) && (y > min_y))                 ?  (&(mm[z+1][y-1][ x ])) : 0;
	neighbors[ 9] = ((z < max_z) && (y < max_y))                 ?  (&(mm[z+1][y+1][ x ])) : 0;
	neighbors[10] = ((z > min_z) &&                (x > min_x))  ?  (&(mm[z-1][ y ][x-1])) : 0;
	neighbors[11] = ((z > min_z) &&                (x < max_x))  ?  (&(mm[z-1][ y ][x+1])) : 0;
	neighbors[12] = ((z < max_z) &&                (x > min_x))  ?  (&(mm[z+1][ y ][x-1])) : 0;
	neighbors[13] = ((z < max_z) &&                (x < max_x))  ?  (&(mm[z+1][ y ][x+1])) : 0;
	neighbors[14] =                ((y > min_y) && (x > min_x))  ?  (&(mm[ z ][y-1][x-1])) : 0;
	neighbors[15] =                ((y > min_y) && (x < max_x))  ?  (&(mm[ z ][y-1][x+1])) : 0;
	neighbors[16] =                ((y < max_y) && (x > min_x))  ?  (&(mm[ z ][y+1][x-1])) : 0;
	neighbors[17] =                ((y < max_y) && (x < max_x))  ?  (&(mm[ z ][y+1][x+1])) : 0;
	neighbors[18] = ((z > min_z) && (y > min_y) && (x > min_x))  ?  (&(mm[z-1][y-1][x-1])) : 0;
	neighbors[19] = ((z > min_z) && (y > min_y) && (x < max_x))  ?  (&(mm[z-1][y-1][x+1])) : 0;
	neighbors[20] = ((z > min_z) && (y < max_y) && (x > min_x))  ?  (&(mm[z-1][y+1][x-1])) : 0;
	neighbors[21] = ((z > min_z) && (y < max_y) && (x < max_x))  ?  (&(mm[z-1][y+1][x+1])) : 0;
	neighbors[22] = ((z < max_z) && (y > min_y) && (x > min_x))  ?  (&(mm[z+1][y-1][x-1])) : 0;
	neighbors[23] = ((z < max_z) && (y > min_y) && (x < max_x))  ?  (&(mm[z+1][y-1][x+1])) : 0;
	neighbors[24] = ((z < max_z) && (y < max_y) && (x > min_x))  ?  (&(mm[z+1][y+1][x-1])) : 0;
	neighbors[25] = ((z < max_z) && (y < max_y) && (x < max_x))  ?  (&(mm[z+1][y+1][x+1])) : 0;
	
  	// essentially same thing as we did for 6-connectivity
	for(i = 0; i<26; i++) {
	  if((neighbors[i]) && (common = (*(neighbors[i]))))
	    // we DO want the value of the assignment
	    break;
	}
	//	std::cerr << "accessed neighbors in mm" << std::endl;
	// i remember reading somewhere NOT to alter the for loop variable...
	// hence the necessity of using "break" as opposed to an inner loop
	if(i<26) { // this means that at least one neighbor already a component
	  i++;
	  while(i<26) {
	    if(neighbors[i] && (temp = (*(neighbors[i])))) {
	      common = myDSF.unite(common, temp);
	    }
	    i++;
	  }
  	  if(boundary) { // if we're at a boundary, can only destroy a component
  	    mm[z][y][x] = 1; // to make a point
	    myDSF.increment_size(1); // what the hell, just increment "size" of "boundary".
	    myDSF.unite(common, 1);
  	  } else {
  	    mm[z][y][x] = common;
	    myDSF.increment_size(common);
  	  }
  	} else {
  	  //  neighborless unbounded component of 0s: DO treat it like a component (already exists)
  	  if(boundary) {
	    mm[z][y][x] = 1;
	    myDSF.increment_size(1);
	  } else {
	    mm[z][y][x] = myDSF.makeSet(gray_it->first);
	  }
  	}
      } // -- done with this set of voxels with the same grayvalue
      myDSF.history_update(HT, gray_it->first); // HT passed by reference
      // at every new update, a new level is created and at least one entry is guaranteed.
      // at least one voxel.  if new, new node.  if joins something, updates that, which now gets new node.
      //      std::cout << "Processed gray " << (gray_it->first) << std::endl;
      if(bv::debug) {
	std::cout << "number of labels in DSF: " << myDSF.get_num_labels() << std::endl;
	std::cout << "number of current components: " << myDSF.get_num_sets() << std::endl;
	std::cout << "number of nodes in history tree: " << HT.numnodes << std::endl;
      }
    } // added all voxels
    if(bv::debug) std::cout << "Done traversing all_voxels for background " <<
		"component tree construction (the map way)" << std::endl;    
    if(bv::debug) std::cout << "setting HT (bgnd) root now..." << std::endl;
    myDSF.set_HT_root(HT);
    if(bv::debug) std::cout << "HT (bgnd) root set." << std::endl;
    //    std::cout << "number of labels in DSF: " << myDSF.get_num_labels() << std::endl;
    //    std::cout << "number of current components: " << myDSF.get_num_sets() << std::endl;
    //    std::cout << "number of nodes in history tree: " << HT.numnodes << std::endl;
    if(bv::debug) {
      std::cout << "Displaying entire non-pruned tree." << std::endl;
      HT.traverse_nonpruned();
    }
    if((minsize | minheight) != 0) {
      prune_history_tree_bgnd(minsize, minheight);
    }
    return myHTBP;
  } // --image3d::construct_history_tree_bgnd_indep()



  static inline bool first_compare(myPairType a, myPairType b) {
    return (a.first < b.first);
  } // --first_compare() [static]


  void bettis() {
    // image3d::threshold and image3d::multilabel_matrix are "implicit parameters"
    // to do:
    // 1) threshold the image, creating multilabel_matrix as necessary (thresholdize)
    thresholdize();
    // 2) call betti_zero, which is responsible for pruning the resultant image
    betti_zero();
    // 3) must call betti_two next.
    betti_two();
    // 4) now for betti_one, which calls euler number and uses previous results
    betti_one();
  } // --image3d::bettis(0-ary)
    
  void bettis(bettiType &b0, bettiType &b1, bettiType &b2) {
    // image3d::threshold and image3d::multilabel_matrix are "implicit parameters"
    // to do:
    // 1) threshold the image, creating multilabel_matrix as necessary (thresholdize)
    thresholdize();
    // 2) call betti_zero, which is responsible for pruning the resultant image
    b0 = betti_zero();
    // 3) must call betti_two next.
    b2 = betti_two();
    // 4) now for betti_one, which calls euler number and uses previous results
    b1 = betti_one();
  } // --image3d::bettis(3-ary)

  /** image3d::persistent(steps) driver, calls the 4-arg version */
  void persistent(unsigned steps = 1, unsigned smallest = 1) {
    bettiType dum1, dum2, dum3;
    this->smallest = smallest;
    persistent(steps, dum1, dum2, dum3);
  } // --image3d::persistent(0 or 1 arg)
  
  /** image3d::persistent(steps and 3 unsigned int refs)
   * Computes a sequence of betti numbers at different thresholds based on steps arg
   * Rationale:
   * We have zeros wherever the image is less than or equal to a threshold,
   *  ones everywhere else.  So to obtain any kind of 'signature', 
   *  it is pointless to look at threshold>=maxOfImage (will always give us all zeros)
   *  since it will not be a discriminating property of our image.
   * However, thresholding at minOfImage *is* (generally) meaningful.
   *  if minOfImage==maxOfImage, we report NO Betti numbers.
   *  it can be thought of as the null sequence, appropriate for constant, featureless images.
   * We COULD run our iterations so that we have one huge connected component (1, 0, 0)
   *  at first (for threshold in the open interval (-inf, minOfImage)) and finally 
   *  (for threshold in (maxOfImage, +inf)) one huge void (0, 0, 0) 
   *  but since every (bounded voronoi) image will have these triplets in its "signature", they
   *  are not distinguishing properties so there is no benefit to thinking of them as such.
   * (Nov 6--Now that I think about it, the b_i references are quite meaningless)
   */
  void persistent(unsigned steps, bettiType &b0, bettiType &b1, bettiType &b2) {
    if(!minmax_set)
      calcMinMax(); // set minOfImage and maxOfImage to correct values
    if (steps==0)
      throw std::invalid_argument("image3d::persistent requires positive, integral 'steps' argument");
    if (maxOfImage == minOfImage) {
      std::cout << 
	"image3d::persistent() warning: min and max values of image are the same!"
		<< std::endl;
    }
    codomain increment = (maxOfImage - minOfImage) / steps;
    if (bv::debug) {
      std::cout << "min: " << minOfImage << ", max: " << maxOfImage << std::endl;
      std::cout << "Increment is " << increment << std::endl;
    }
    if (increment==0) // could happen if our codomain is of 'integral' type (e.g., bool)
      increment = 1;
    bettiType prev_b0 = 1, prev_b1 = 0, prev_b2 = 0;
    // not very idiomatic 'C', but a bit more userfri3ndly: step 1 happens at minOfImage
    // October: 'fixed the loop' so always starts with "1 0 0" and ends with "0 0 0"
    threshold = minOfImage;
    // the following is based on our "domain knowledge: our entire box will be 
    // one big component after thresholding below minimum value of image
    std::cout << "1 0 0"
	      << " (thresh < " << minOfImage << ")" << std::endl;
    for(unsigned step = 1; step <= steps; step++) {
      // threshold is a data member that our algorithms in this class use
      bettis(b0, b1, b2);
      if ((b0 == prev_b0) && (b1 == prev_b1) && (b2 == prev_b2)) {
	  // don't print anything
      } else {
	if(0 != b0) {
	  std::cout << b0 << " " << b1 << " " << b2 
		    << " (step " << step << ", thresh@" << threshold << " )" << std::endl;
	} else {
	  std::cout << b0 << " " << b1 << " " << b2 
		    << " (step " << step << ", thresh >= " << threshold << " )" << std::endl;
	  return;
	}
      }
      prev_b0 = b0;
      prev_b1 = b1;
      prev_b2 = b2;
      threshold+=increment;
    } // --for
    if (b0 != 0) {
      // b0 *could* have been zero for thresh < maxval
      // if all components after thresholding turned out to be smaller than allowed
      // we also have the 'domain knowledge' that b0 == 0 forces b1 = b2 = 0
      // (very buddhic: no component ==> no hole, no cavity)
      std::cout << "0 0 0"
	      << " (thresh >= " << maxOfImage << ")" << std::endl;
    }
    // so, even if steps was 1, we print at least 2 sets of persistent Betti numbers
  } // --image3d::persistent(4args)

  /** calculate and display what we think are also true sequences
      first get a list of all extrema, i.e., values at voxels that are 
      max or min of 6-neighbors... then, interpolate so that 
      thresholding based on > vs >= does not make a difference
      Caveat: We don't know what this does to tunnels, we THINK it doesn't affect them
      ... because they appear only at non-degenerate critical points (according to morse theory)
      This whole thing fell through (so far) because critical points aren't just extrema.
      DON'T use this!
  */

  
  /** calculate and display the *true* sequence of Betti numbers *
   *  threshold at all grayvalues... determine that set first.   *
   *  otherwise much like "persistent".                          *
   */
  void true_sequence(unsigned smallest = 1) {
    this->smallest = smallest;
    if (!(matrix_set && image_data)) {
      std::cerr << "Matrix not set, and trying to find true_sequence() !" << std::endl;
      std::cerr << "image3d::true_sequence() returning control to caller." << std::endl;
      return;
    }
    // debug
    std::cout << "in image3d::true_sequence()" << std::endl;
    codomain *x_ptr = image_data->getFirst();
    codomain *x_bound = image_data->getPast();
    std::set<codomain> grays; // all gray values go here
    while(x_ptr != x_bound) {
      grays.insert(*x_ptr);
      x_ptr++;
    }
    // STL sets (containing numeric type) are always sorted.

    // minOfImage and maxOfImage are instance variables
    maxOfImage = (*(grays.rbegin()));
    // underlying assumption throughout: image is at least 1x1x1, 
    // which means that grays must contain at least 1 value
    typename std::set<codomain>::const_iterator it = grays.begin();
    minOfImage = (*it);
    bettiType b0=1, b1=0, b2=0; // locals overriding the instance variables
    bettiType prev_b0 = b0, prev_b1 = b1, prev_b2 = b2;
    std::cout << "1 0 0" << " (thresh < " << minOfImage << ")" << std::endl;
    if(minOfImage != maxOfImage) {
      //    codomain grayval;
      for( ; it != grays.end(); it++) {
	// threshold is an instance variable
	threshold = (*it);
	//  thresholdize() is called by bettis()
	if (threshold < maxOfImage) {
	  bettis(b0, b1, b2); // passed by reference
	  if ((b0 != prev_b0) || (b1 != prev_b1) || (b2 != prev_b2)) {
	    std::cout << b0 << " " << b1 << " " << b2
		      << " (thresh @ " << threshold << ")" << std::endl;
	    prev_b0 = b0; prev_b1 = b1; prev_b2 = b2;
	  }
	}
      }
    }
    // must have at least one component whenever threshold is less than maxOfImage
    // so we are guaranteed that the tuple differs from "0 0 0"
    std::cout << "0 0 0"
	      << " (thresh >= " << maxOfImage << " )" << std::endl;
  }  // --image3d::true_sequence()
  
private:

  /* The INCREMENTAL version of the euler number calculation (Yung's Dec 9 way)   *
   * Rather than processing for all points of the multilabel matrix,              *
   * returns the change to euler number from coordinate of a voxel that has       *
   * become nonzero in multilabel_matrix since the last call to euler_incremental *
   *   (i.e., this function ASSUMES that the passed coord WAS ZERO last time!)    *
   * pre: instance variable E initialized to 0 BEFORE the first call ever,        *
   *  !!! NOT to be reinitilized here !!!                                         *
   * big side effect: UPDATE E.  E is thusly available to the caller, as well     *
   * as the change to E, which is the return value.                               */
  eulerType euler_incremental(const coord coo) {
    if(bv::debug)
      std::cout << "Running image3d::euler_incremental()..." << std::endl;
    dimIndexType z, y, x;
    labelType*** mm = multilabel_matrix->getOrigin();
    // YASC:
    if(0 == mm[coo.zpos][coo.ypos][coo.xpos]) {
      std::cerr << "Freakout!" << std::endl;
    }
    static const euIndexType v0=0x01, v1=0x02, v2=0x04, v3=0x08;
    static const euIndexType v4=0x10, v5=0x20, v6=0x40, v7=0x80;
    // convention: vk = 2 to the power of k.  each voxel is a bit of the index.

    /// here is how we think of our 2x2x2 subpatterns:
    ///
    ///         increasing z (behind)
    ///            ^
    ///           /
    ///          v4-------v5---> increasing x (to the right of)
    ///         /        / |
    ///        /        /  |
    ///       v0------v1   |
    ///       |       |   v7
    ///       |       |  /    (v6 is not shown so that clutter may be reduced,
    ///       |       | /          you can picture it below v4, behind v2.)
    ///       v2------v3
    ///       |
    ///       |
    ///       V
    ///      increasing y (below)
    ///
    /// sort of in line with our notion of "z varies fastest, then y, then x"

    // the Euler differential for each possible 2x2x2 (sub)pattern qij
    //  ['q' preferred to 'Q' since it looks less like a zero]
    //   as in Toriwaki & Yonekura (p.190) initialized based on 6-connectivity
    // i is the number of 1-voxels at vertex
    // j is consistent with anti-symmetry mentioned below
    //    ==> i has to match the number of "on" (1) bits in the index
    // some symmetries (for sanity checks):
    //   swap high 4 bits & low 4 bits:
    //   [ reflection of 2x2x2 graph to flip z coord bits | all vk=v((k+4) mod 8) ]
    //      ==> symmetry in "main (major) diagonal" of 16x16 table
    //   qij = "-" q(8-i)j :  ["-" means reversing the 1- and 0-voxels]
    //                    [ note that this is valid for i=4 as well, so, q4j = "-" q4j ]
    //      ==> (anti-)symmetry in "origin" (indicated by an 'o' in a comment) of 16x16 table
    //   (combined) ==> (anti-)symmetry in *minor* diagonal of 16x16 table

    static const euLocalType q01 =  0;
    static const euLocalType q11 =  1;
    static const euLocalType q21 =  0,  q22 =  2,  q23 =  2;
    static const euLocalType q31 = -1,  q32 =  1,  q33 =  3; // q33 was (wrongly) 0
    static const euLocalType q41 =  0,  q42 = -2,  q43 = -2, q44 = 0, q45 = 0, q46 = 4;
    static const euLocalType q51 = -1,  q52 = -3,  q53 = -1;
    static const euLocalType q61 =  0,  q62 = -2,  q63 = -6;
    static const euLocalType q71 =  1;
    static const euLocalType q81 =  0;

    static const euLocalType euTab[256] =
     // number of "on" (1) bits among the four low bits of index:
     // 0   1   1   2   1   2   2   3 | 1   2   2   3   2   3   3   4
     // -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
     // low 4 bits of index (v3+...+v0):
     // 0   1   2   3   4   5   6   7 | 8   9   A   B   C   D   E   F
     // -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
    /*                                |                                */
      {q01,q11,q11,q21,q11,q21,q22,q31,q11,q22,q21,q31,q21,q31,q31,q41,  // 0 | high      | 0

       q11,q21,q22,q31,q22,q31,q33,q42,q23,q32,q32,q43,q32,q43,q44,q51,  // 1 | 4 bits    | 1

       q11,q22,q21,q31,q23,q32,q32,q43,q22,q33,q31,q42,q32,q44,q43,q51,  // 2 |  of       | 1

       q21,q31,q31,q41,q32,q43,q44,q51,q32,q44,q43,q51,q45,q52,q52,q61,  // 3 | index     | 2

       q11,q22,q23,q32,q21,q31,q32,q43,q22,q33,q32,q44,q31,q42,q43,q51,  // 4 |(v7+...+v4)| 1

       q21,q31,q32,q43,q31,q41,q44,q51,q32,q44,q45,q52,q43,q51,q52,q61,  // 5 |  <---<    | 2

       q22,q33,q32,q44,q32,q44,q45,q52,q33,q46,q44,q53,q44,q53,q52,q62,  // 6 |           | 2

       q31,q42,q43,q51,q43,q51,q52,q61,q44,q53,q52,q62,q52,q62,q63,q71,  // 7 |           | 3
/*  ---                               o                              --- */
       q11,q23,q22,q32,q22,q32,q33,q44,q21,q32,q31,q43,q31,q43,q42,q51,  // 8 |    >--->  | 1

       q22,q32,q33,q44,q33,q44,q46,q53,q32,q45,q44,q52,q44,q52,q53,q62,  // 9 | number    | 2

       q21,q32,q31,q43,q32,q45,q44,q52,q31,q44,q41,q51,q43,q52,q51,q61,  // A | of "on"   | 2

       q31,q43,q42,q51,q44,q52,q53,q62,q43,q52,q51,q61,q52,q63,q62,q71,  // B | bits      | 3

       q21,q32,q32,q45,q31,q43,q44,q52,q31,q44,q43,q52,q41,q51,q51,q61,  // C | among     | 2

       q31,q43,q44,q52,q42,q51,q53,q62,q43,q52,q52,q63,q51,q61,q62,q71,  // D | high      | 3

       q31,q44,q43,q52,q43,q52,q52,q63,q42,q53,q51,q62,q51,q62,q61,q71,  // E | 4 bits    | 3

       q41,q51,q51,q61,q51,q61,q62,q71,q51,q62,q61,q71,q61,q71,q71,q81}; // F |           | 4
     /*                               |                                 */
    // euTab tells us the local euler number for an 8-bit index (treated as a bitfield)
    // (value at vertex vi corresponds to the bit that represents 2^i as an unsigned)

    // we simulate a coat of 0s around our entire image
    // without enhancing / changing it or allocating a new one.
    bool ok0, ok1, ok2, ok3, ok4, ok5, ok6, ok7;
    euIndexType curin;
    eulerType delta_E = 0; // the change in E for the addition of this new voxel
    // calculate change due to changes at the 8 2x2x2 "subcubes" of the 
    // 3x3x3 cube centered around coo
    for(z = -1 + coo.zpos; z <= coo.zpos; z++) {
      for(y = -1 + coo.ypos; y <= coo.ypos; y++) {
	for(x = -1 + coo.xpos; x <= coo.xpos; x++) {
	  //	  std::cout << z << " " << y << " " << x << ";";
	  ok0 = ok1 = ok2 = ok3 = ok4 = ok5 = ok6 = ok7 = true;
	  if(z <  min_z) ok0 = ok1 = ok2 = ok3 = false;
	  if(z == max_z) ok4 = ok5 = ok6 = ok7 = false;
	  if(y <  min_y) ok0 = ok1 = ok4 = ok5 = false;
	  if(y == max_y) ok2 = ok3 = ok6 = ok7 = false;
	  if(x <  min_x) ok0 = ok2 = ok4 = ok6 = false;
	  if(x == max_x) ok1 = ok3 = ok5 = ok7 = false;

	  curin = 0;

	  if(ok0 && (0!= mm[ z ][ y ][ x ])) curin+=v0;
	  if(ok1 && (0!= mm[ z ][ y ][x+1])) curin+=v1;
	  if(ok2 && (0!= mm[ z ][y+1][ x ])) curin+=v2;
	  if(ok3 && (0!= mm[ z ][y+1][x+1])) curin+=v3;
	  if(ok4 && (0!= mm[z+1][ y ][ x ])) curin+=v4;
	  if(ok5 && (0!= mm[z+1][ y ][x+1])) curin+=v5;
	  if(ok6 && (0!= mm[z+1][y+1][ x ])) curin+=v6;
	  if(ok7 && (0!= mm[z+1][y+1][x+1])) curin+=v7;

	  delta_E += (euTab[curin]);
	  //	  std::cout << "curin: " << static_cast<unsigned int>(curin);
	  // we KNOW that the vi correponding to coo WAS in curin
	  // (because the multilabel_matrix WAS UPDATED before call to this function)
	  // and now shall be taken away
	  euIndexType vi = (0x01 << (4*(coo.zpos - z) + 2*(coo.ypos - y) + (coo.xpos - x)));
	  //	  std::cout << static_cast<unsigned int>(vi) << " ";
	  curin -= vi;
	  //	  std::cout << " and then: " << static_cast<unsigned int>(curin) << std::endl;
	  delta_E -= (euTab[curin]);
	  // subtract from delta_E the "previous" (before coo) euler number for this 2x2x2 cube
	} // --for x
      } // --for y
    } // --for x
    if (0!= (delta_E%8)) {
      std::cerr << "in image3d::euler_incremental(): delta_E is " << delta_E << std::endl;
      std::cerr << "NOT a multiple of 8." << std::endl;
    }
    delta_E /= 8; // refrain from ">>3" since E can be negative
    //    std::cout << delta_E << " " << E;
    (this->E) += delta_E; // maintain the instance variable E as always the ACTUAL Euler number (not times 8)
    //    std::cout << std::endl;
    return delta_E;
  } // --image3d::euler_incremental()

  /** Calculate the euler number of 1s based on 6/26 connectivity
   * for the thresholded image in the instance variable multilabel_matrix   
   */
  eulerType euler() {
    if(bv::debug)
      std::cout << "Running image3d::euler()..." << std::endl;
    dimIndexType z, y, x;
    //    codomain*** im = image_data->getOrigin();
    labelType*** mm = multilabel_matrix->getOrigin();
    static const euIndexType v0=0x01, v1=0x02, v2=0x04, v3=0x08;
    static const euIndexType v4=0x10, v5=0x20, v6=0x40, v7=0x80;
    // convention: vk = 2 to the power of k.  each voxel is a bit of the index.

    /// here is how we think of our 2x2x2 subpatterns:
    ///
    ///         increasing z (behind)
    ///            ^
    ///           /
    ///          v4-------v5---> increasing x (to the right of)
    ///         /        / |
    ///        /        /  |
    ///       v0------v1   |
    ///       |       |   v7
    ///       |       |  /    (v6 is not shown so that clutter may be reduced,
    ///       |       | /          you can picture it below v4, behind v2.)
    ///       v2------v3
    ///       |
    ///       |
    ///       V
    ///      increasing y (below)
    ///
    /// sort of in line with our notion of "z varies fastest, then y, then x"

    ///
    // the Euler differential for each possible 2x2x2 (sub)pattern qij
    //  ['q' preferred to 'Q' since it looks less like a zero]
    //   as in Toriwaki & Yonekura (p.190) initialized based on 6-connectivity
    // i is the number of 1-voxels at vertex
    // j is consistent with anti-symmetry mentioned below
    //    ==> i has to match the number of "on" (1) bits in the index
    // some symmetries (for sanity checks):
    //   swap high 4 bits & low 4 bits:
    //   [ reflection of 2x2x2 graph to flip z coord bits | all vk=v((k+4) mod 8) ]
    //      ==> symmetry in "main (major) diagonal" of 16x16 table
    //   qij = "-" q(8-i)j :  ["-" means reversing the 1- and 0-voxels]
    //                    [ note that this is valid for i=4 as well, so, q4j = "-" q4j ]
    //      ==> (anti-)symmetry in "origin" (indicated by an 'o' in a comment) of 16x16 table
    //   (combined) ==> (anti-)symmetry in *minor* diagonal of 16x16 table

    static const euLocalType q01 =  0;
    static const euLocalType q11 =  1;
    static const euLocalType q21 =  0,  q22 =  2,  q23 =  2;
    static const euLocalType q31 = -1,  q32 =  1,  q33 =  3; // q33 was (wrongly) 0
    static const euLocalType q41 =  0,  q42 = -2,  q43 = -2, q44 = 0, q45 = 0, q46 = 4;
    static const euLocalType q51 = -1,  q52 = -3,  q53 = -1;
    static const euLocalType q61 =  0,  q62 = -2,  q63 = -6;
    static const euLocalType q71 =  1;
    static const euLocalType q81 =  0;

    //    int i= q01+q11+q21+q22+q23+q31+q32+q33+q41+q42+q43+q44+q45+q46;
    // i+= q51+q52+q53+q61+q62+q63+q71+q81;
    static const euLocalType euTab[256] =
     // number of "on" (1) bits among the four low bits of index:
     // 0   1   1   2   1   2   2   3 | 1   2   2   3   2   3   3   4
     // -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
     // low 4 bits of index (v3+...+v0):
     // 0   1   2   3   4   5   6   7 | 8   9   A   B   C   D   E   F
     // -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
    /*                                |                                */
      {q01,q11,q11,q21,q11,q21,q22,q31,q11,q22,q21,q31,q21,q31,q31,q41,  // 0 | high      | 0

       q11,q21,q22,q31,q22,q31,q33,q42,q23,q32,q32,q43,q32,q43,q44,q51,  // 1 | 4 bits    | 1

       q11,q22,q21,q31,q23,q32,q32,q43,q22,q33,q31,q42,q32,q44,q43,q51,  // 2 |  of       | 1

       q21,q31,q31,q41,q32,q43,q44,q51,q32,q44,q43,q51,q45,q52,q52,q61,  // 3 | index     | 2

       q11,q22,q23,q32,q21,q31,q32,q43,q22,q33,q32,q44,q31,q42,q43,q51,  // 4 |(v7+...+v4)| 1

       q21,q31,q32,q43,q31,q41,q44,q51,q32,q44,q45,q52,q43,q51,q52,q61,  // 5 |  <---<    | 2

       q22,q33,q32,q44,q32,q44,q45,q52,q33,q46,q44,q53,q44,q53,q52,q62,  // 6 |           | 2

       q31,q42,q43,q51,q43,q51,q52,q61,q44,q53,q52,q62,q52,q62,q63,q71,  // 7 |           | 3
/*  ---                               o                              --- */
       q11,q23,q22,q32,q22,q32,q33,q44,q21,q32,q31,q43,q31,q43,q42,q51,  // 8 |    >--->  | 1

       q22,q32,q33,q44,q33,q44,q46,q53,q32,q45,q44,q52,q44,q52,q53,q62,  // 9 | number    | 2

       q21,q32,q31,q43,q32,q45,q44,q52,q31,q44,q41,q51,q43,q52,q51,q61,  // A | of "on"   | 2

       q31,q43,q42,q51,q44,q52,q53,q62,q43,q52,q51,q61,q52,q63,q62,q71,  // B | bits      | 3

       q21,q32,q32,q45,q31,q43,q44,q52,q31,q44,q43,q52,q41,q51,q51,q61,  // C | among     | 2

       q31,q43,q44,q52,q42,q51,q53,q62,q43,q52,q52,q63,q51,q61,q62,q71,  // D | high      | 3

       q31,q44,q43,q52,q43,q52,q52,q63,q42,q53,q51,q62,q51,q62,q61,q71,  // E | 4 bits    | 3

       q41,q51,q51,q61,q51,q61,q62,q71,q51,q62,q61,q71,q61,q71,q71,q81}; // F |           | 4
     /*                               |                                 */
    // euTab tells us the local euler number for an 8-bit index (treated as a bitfield)
    // (value at vertex vi corresponds to the bit that represents 2^i as an unsigned)

    // we simulate a coat of 0s around our image
    // without enhancing / changing it or allocating a new one.
    bool ok0, ok1, ok2, ok3, ok4, ok5, ok6, ok7;
    euIndexType curin;
    // eulerType E (unsigned) is a data member reserved for this.
    E = 0; // needs to happen every time euler() is called.
    for(z=min_z-1; z<=max_z; z++) {
      for(y=min_y-1; y<=max_y; y++) {
	for(x=min_x-1; x<=max_x; x++) {
	  ok0 = ok1 = ok2 = ok3 = ok4 = ok5 = ok6 = ok7 = true;
	  if(z <  min_z) ok0 = ok1 = ok2 = ok3 = false;
	  if(z == max_z) ok4 = ok5 = ok6 = ok7 = false;
	  if(y <  min_y) ok0 = ok1 = ok4 = ok5 = false;
	  if(y == max_y) ok2 = ok3 = ok6 = ok7 = false;
	  if(x <  min_x) ok0 = ok2 = ok4 = ok6 = false;
	  if(x == max_x) ok1 = ok3 = ok5 = ok7 = false;

	  curin = 0;
	  if(ok0 && (0!= mm[ z ][ y ][ x ])) curin+=v0;
	  if(ok1 && (0!= mm[ z ][ y ][x+1])) curin+=v1;
	  if(ok2 && (0!= mm[ z ][y+1][ x ])) curin+=v2;
	  if(ok3 && (0!= mm[ z ][y+1][x+1])) curin+=v3;
	  if(ok4 && (0!= mm[z+1][ y ][ x ])) curin+=v4;
	  if(ok5 && (0!= mm[z+1][ y ][x+1])) curin+=v5;
	  if(ok6 && (0!= mm[z+1][y+1][ x ])) curin+=v6;
	  if(ok7 && (0!= mm[z+1][y+1][x+1])) curin+=v7;
	  if(bv::debug) {
	    if((-1==z) && (-1==y) && (-1==x)) {
	      std::cout << "2x2x2 at " << z << ", " << y << ", " << x;
	      std::cout << "; curin: " << std::hex << curin;
	      std::cout << " (" << std::dec << curin << ")" 
			<< "; local euler: " << euTab[curin] << std::endl;
	    }
	  }
	  E += euTab[curin];
	} // --for x
      } // --for y
    } // --for z
    if (0!= (E%8)) {
      std::cerr << "in image3d::euler(): Euler number 'times 8' is " << E << std::endl;
      std::cerr << "NOT a multiple of 8." << std::endl;
    }
    return E/8; // refrain from ">>3" since E can be negative
  } // --image3d::euler()
  
  int unset_matrix() { 
    if(matrix_set) {
      if(image_data) { 
				delete image_data; // delete calls matrix3d destructor on *image_data
				image_data = 0;
      } else {
	if(bv::debug)
	  std::cerr << "We Have a Problem! [ image3d::unset_matrix() ]" 
		    << std::endl;
	return 1;
      }
      matrix_set = false;
    }
    else {
      std::cout << "Matrix not set to begin with!" << std::endl;
    }
    return 0;
  } // --image3d::unset_matrix()
  
public:
  /** dump image to stdout (can take a few days) */
  void dump() {    
    dimIndexType z, y, x;
    codomain*** im = getOrigin();
    for(z=min_z; z<=max_z; z++) {
      std::cout << "For  z = " << z << ":" << std::endl;
      for(y=min_y; y<=max_y; y++) {
	for(x=min_x; x<=max_x; x++) {
	  std::cout << "im@["<<z<<"]["<<y<<"]["<<x<<"] = "
		    << im[z][y][x] << std::endl;
	} // --for x
      } // --for y
    } // --for z
    std::cout << std::endl;
  } // --image3d::dump()

  /** Dumps a central "cube" based on the "radius" parameter. *
   *  No range checking--use at your own risk!                */
  void dump_core(int rad) {
    dimIndexType z, y, x;
    codomain*** im = getOrigin();
    for(z= -rad; z<=rad; z++) {
      std::cout << "For  z = " << z << ":" << std::endl;
      for(y= -rad; y<=rad; y++) {
	for(x= -rad; x<=rad; x++) {
	  std::cout << "im@["<<z<<"]["<<y<<"]["<<x<<"] = "
		    << im[z][y][x] << std::endl;
	} // --for x
      } // --for y
    } // --for z
    std::cout << std::endl;
  } // --image3d::dump()

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
  } // --image3d::is_all_zero()

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
  } // --image3d::is_all_nonzero()

  ~image3d() {
    if(bv::debug)
      std::cout << "Destroying image3d" << std::endl;
    if(matrix_set) 
      unset_matrix();
    if(multilabel_matrix)
      delete multilabel_matrix;
    // may want to destruct DSF(s) and historytree(s) somewhere else after assigning
    // well, let these stay until that's decided.
//      if(myDSFpb0) delete myDSFpb0;
//      if(myDSFpb2) delete myDSFpb2;
    if(FGDSF) delete FGDSF;
    if(BGDSF) delete BGDSF;
    if(HTP) delete HTP;
    if(HTBP) delete HTBP;
  }
}; // -- image3d<codomain> class template definition

