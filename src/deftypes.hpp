#ifndef _DEFTYPES_DEFINED
#define _DEFTYPES_DEFINED
// (relatively) simple types, used by image3d and others
typedef short int dimIndexType;
typedef unsigned short int dimSizeType; // has to be "unsigned dimIndexType"
typedef unsigned long long int totSizeType; // must handle max_dimSize cubed
typedef unsigned int labelType; // should be 32 bits, up it if not so
// labelType has to be comparable using '<'

// labelType has to be comparable using '<'
typedef unsigned int bettiType; // for all our Betti numbers
typedef unsigned char euIndexType; // indices to euler table need only be 8 bits
typedef signed long int eulerType;     // for the Euler number, which can be negative
typedef eulerType euLocalType; // same so that we can trust the signed arithmetic

// holds the coordinates of a voxel, for Yung's (Dec 2003) way of calculating Betti numbers
struct coord {
  dimIndexType zpos;
  dimIndexType ypos;
  dimIndexType xpos;
  // need to declare REAL coordinate type.
/*    static coord mean(coord a, unsigned weight_a,  */
/*  		    coord b, unsigned weight_b) {  */
/*      // by value: caller should make sure to declare own var */
/*      unsigned totalweight = weight_a + weight_b; */
/*      coord av; */
/*      av.zpos = (weight_a * a.zpos + weight_b * b.zpos)*1.0/totalweight; */
/*      av.ypos = (weight_a * a.ypos + weight_b * b.ypos)*1.0/totalweight; */
/*      av.xpos = (weight_a * a.xpos + weight_b * b.xpos)*1.0/totalweight; */
/*      return av; */
/*    } */
};

struct fcoord {
  float zpos;
  float ypos;
  float xpos;
  // define operator=
  static fcoord mean(fcoord a, unsigned weight_a, 
		    fcoord b, unsigned weight_b) { 
    // by value: caller should make sure to declare own var
    unsigned totalweight = weight_a + weight_b;
    fcoord av;
    //    av.zpos = (weight_a * a.zpos + weight_b * b.zpos)*1.0/totalweight;
    //    av.ypos = (weight_a * a.ypos + weight_b * b.ypos)*1.0/totalweight;
    //    av.xpos = (weight_a * a.xpos + weight_b * b.xpos)*1.0/totalweight;
    av.zpos = (weight_a * a.zpos + weight_b * b.zpos)/totalweight;
    av.ypos = (weight_a * a.ypos + weight_b * b.ypos)/totalweight;
    av.xpos = (weight_a * a.xpos + weight_b * b.xpos)/totalweight;
    return av;
  }
};

/* dimIndexType is currently unsigned short, so a coord should take 6 bytes
   a pointer is 4 bytes on these machines, so not too much is being lost by this approach.
   has the advantage of being "transparent": 
   can be used "directly" both for image_data as well as associated label images
*/


#endif
