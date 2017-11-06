#ifndef PERBIN_H
#define PERBIN_H
/** @package botanist
    @date 2005.06.22
    @author Deniz Sarioz
    @description types to internally organize/represent many trees at once
    all of the "our 2004 pruning" variety (for now)
*/
#include<map>
#include <QString>
/** a more extensible parameter list structure */
struct params_t {
  int minsize;
  int delta;
  QString str; // string rep of all parameters

  /** define less than in lexicographic order */  
  bool operator<(params_t other) const {	
    if (minsize < other.minsize) return true;
    if (minsize > other.minsize) return false;
    return (delta < other.delta);
  }

  params_t(int minsize_in,int delta_in) {
    minsize = minsize_in;
    delta = delta_in;
    str = QString( "(%1,%2)" ) . arg(minsize). arg(delta);
  }

  params_t() {
    params_t(0,0);
  }

};

typedef std::map<params_t, edgytree*> edgies_t;
typedef std::map<params_t, historytree*> histories_t;

struct perbinsize_t { 
  historytree*      htfgp;
  historytree_bgnd* htbgp;
  histories_t fg_histories;
  edgies_t fg_edgies;
  edgies_t bg_edgies;

  /** perbinsize_t constructor performs some inits to be safe 
   */
  perbinsize_t() {
    htfgp = 0;
    htbgp = 0;
    fg_edgies.clear();
    bg_edgies.clear();
  } // --struct perbinsize_t ctor
  
  /** perbinsize_t destructor to guard against memory leaks
      pre: all of these things were initialized with new, and exist where expected
  */
  ~perbinsize_t() {
    if(htfgp) {
      delete htfgp;
      // std::cout << "delete'ing htfgp of bin" << std::endl;
    }
    if(htbgp) {
      delete htbgp;
      // std::cout << "delete'ing htbgp of bin" << std::endl;
    }
    // recursively delete edgy trees
    for(edgies_t::iterator it = fg_edgies.begin(); 
			it != fg_edgies.end(); it++) {
      //std::cout << "deleting associated FG edgytree pruned by " + it->first.str << std::endl;
      if(it->second) delete (it->second);
    }
    for(edgies_t::iterator it = bg_edgies.begin(); 
			it != bg_edgies.end(); it++) {
      //std::cout << "deleting associated BG edgytree pruned by " << it->first.str << std::endl;
      if(it->second) delete (it->second);
    }
    // recursively delete histories
    for(histories_t::iterator it = fg_histories.begin();
            it != fg_histories.end(); it++) {
      //std::cout << "deleting associated FG edgytree pruned by " + it->first.str << std::endl;
      if(it->second) delete (it->second);
    }
  } // --struct perbinsize_t dtor

}; // --struct perbinsize_t def



/** "smart array" allmybins_t */
struct allmybins_t {
  //perbinsize_t* _bins[256]; /** < array of pointers to perbinsize_t */
    std::map<int, perbinsize_t*> binsMap;
  //~ /** constructor makes sure to set all ptrs to null */
  //~ allmybins_t() {
    //~ for(int i = 0; i < 256; i++) {
      //~ _bins[i] = 0;
    //~ }
  //~ }  // --allmybins_t ctor

  /** accessor */
  perbinsize_t* operator[](int i) {	
    std::map<int, perbinsize_t*>::iterator it = binsMap.find(i);
         if(bv::debug) printf("perbinsize_t* operator[]  %d \n",binsMap.size());
	
		if (it!= binsMap.end() ){
             if(bv::debug) printf(" it->first %d\n",it->first);
			return it->second;
		}else{
             if(bv::debug) printf(" NTEM %d \n",i);
			return 0;	
		}
  }  // --allmybins_t()[] accessor

  void makeNewBin(int i){
        if(bv::debug) printf("perbinsize_t* operator[]  %d \n",binsMap.size());
		std::map<int, perbinsize_t*>::iterator it = binsMap.find(i);
		if (it!= binsMap.end() ){	
             if(bv::debug) printf(" CANOT added bin %d. This bin already exist \n");
		}else{
           if(bv::debug) printf(" will added bin %d \n",i);
            binsMap[i]= new perbinsize_t();
		}
	  //~ if(!_bins[i]){
      //~ _bins[i] = new perbinsize_t;
      //~ std::cout << "added bin " << i << std::endl;
    //~ }
  } // --allmybins_t::makeNewBin()


  /** maybe create and then return (that) new edgy tree
      @precondition bin binno already initialized
      @precondition htfgp and htbgp of bin binno already initialized to MEANINGFUL things
  */
edgytree* get_edgy(int binno, int minsize, int delta) {
	perbinsize_t* curbinp = binsMap[binno];

	if(! curbinp) {
		std::cout << "Error, no bin in allmybins_t::get_edgy()" << binno << std::endl;
		return 0;
	}
    std::cout << "get_edgy():: parameter (" <<	binno << ", minsize_i= " << minsize<< " delta="<< (delta) << std::endl;
	
    params_t mytuple (minsize,delta);
	if(true) { // based on (foreground) historytree
		if(! curbinp->htfgp) {
			std::cout << "Error, no htfgp set for binno " << binno << std::endl;
			return 0;
		} // else...		
		edgies_t& fgeds = curbinp->fg_edgies;
		printf(" size => %d \n",fgeds.size());	
		if( fgeds.find(mytuple) != fgeds.end()) { // found
			return fgeds[mytuple];
		} else { // create new
            std::cout << "Ta errado. Nao achou o parameter for bins (" <<	binno << ", minsize_i= " << minsize<< "  delta="<< static_cast<unsigned int>  (delta) << std::endl;
			printf("Ta errado. Nao achou o parameter for bins :-) \n");	
            exit(0); // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			curbinp->htfgp->prune(minsize, 0);
			edgytree* etp = curbinp->htfgp->generate_pruned_edgy_tree();
			curbinp->htfgp->prune(0,0); // is this hack (still) necessary to be safe?
			etp->flatten();
            etp->pruneByHeight(delta);
			etp->flatten(); // maybe overkill?  but playing it safe.
			fgeds[mytuple] = etp;
			return etp;
		}
	}else{ // based on (background) historytree_bgnd
		if(! curbinp->htbgp) {
			std::cout << "Error, no htbgp set for binno " << binno << std::endl;
			return 0;
		} // else...
		edgies_t& bgeds = curbinp->bg_edgies;
		if( bgeds.find(mytuple) != bgeds.end()) { // found
			return bgeds[mytuple];
		} else { // create new
			curbinp->htbgp->prune(minsize, 0);
			edgytree* etp = curbinp->htbgp->generate_pruned_edgy_tree();
			curbinp->htbgp->prune(0,0); // is this hack (still) necessary to be safe?
			etp->flatten();
            etp->pruneByHeight(delta);
			etp->flatten(); // maybe overkill?  but playing it safe.
			bgeds[mytuple] = etp;
			return etp;
		}
	}
} // --allmybins_t::get_edgy()


/**


*/
historytree* get_historyTree(int binno, int minsize, int delta) {
  perbinsize_t* curbinp = binsMap[binno];

  if(! curbinp) {
      std::cout << "Error, no bin in allmybins_t::get_historyTree()" << binno << std::endl;
      return 0;
  }
  std::cout << "get_historyTree():: parameter (" <<	binno << ", minsize_i= " << minsize<< " delta="<< (delta) << std::endl;

  params_t mytuple (minsize,delta);
  if(true) { // based on (foreground) historytree
      if(! curbinp->htfgp) {
          std::cout << "Error, no htfgp set for binno " << binno << std::endl;
          return 0;
      } // else...
      histories_t& fght = curbinp->fg_histories;
      printf(" size => %d \n",fght.size());
      if( fght.find(mytuple) != fght.end()) { // found
          return fght[mytuple];
      } else { // create new
          std::cout << "Ta errado in get_historyTree(). Nao achou o parameter for bins (" <<	binno << ", minsize_i= " << minsize<< "  delta="<< static_cast<unsigned int>  (delta) << std::endl;
          printf("Ta errado. Nao achou o parameter for bins :-) \n");
          exit(0); // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      }
  }
} // --allmybins_t::get_historyTree()

  /** reset everything! (formerly in destructor, this is more elegant)
   */
  void reset() {
 			std::cout << "ATTENTION destructing perbinsize_t object for #bins=: ";
			exit(0);
	 } // --allmybins_t::reset()


  ~allmybins_t() {
    reset();
  } // --allmybins_t dtor  

};

struct allmydesciptors_t {

    std::map<QString, allmybins_t*> structuresMap;

    /** accessor */
    allmybins_t* operator[](QString i) {
        std::map<QString, allmybins_t*>::iterator it = structuresMap.find(i);
        if(bv::debug) printf("allmydesciptors_t* operator[]  %d \n",structuresMap.size());
        if (it!= structuresMap.end() ){
            cout << " it->first "<< it->first.toAscii().constData()<<"\n";
            return it->second;
        }else{
            printf(" NTEM\n");
            return 0;
        }
    }// --allmydesciptors_t()[] accessor

    void makeNewStructures(QString i){
        if(bv::debug) cout<< " makeNewStructures for "<<i.toAscii().constData()<<"\n";
        std::map<QString, allmybins_t*>::iterator it = structuresMap.find(i);
        if (it!= structuresMap.end() ){
            if(bv::debug) printf(" CANOT added bin %d. This bin already exist \n");
        }else{
            if(bv::debug) cout<< "will added structures "<<i.toAscii().constData()<<"\n";
            structuresMap[i]= new allmybins_t();
        }
    }

    void reset() {
       std::cout << "ATTENTION destructing allmydesciptors_t object for #bins=: ";
       exit(0);
    }// --allmybins_t::reset()

    ~allmydesciptors_t() {
      reset();
    } // --allmydesciptors_t dtor
};

#endif // PERRBIN_H
