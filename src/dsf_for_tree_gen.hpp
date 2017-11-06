/**
    @file dsf_for_tree_gen.cpp
    @date 2004.06.20--28
    back again to DSF WITH union-by-rank AND path-compression
    this is to be used as "helper" to tree (of history of components etc) generation
    @author Deniz for both this and preds.
    @implements class dsf_for_tree_gen
    @description
    Disjoint-set forest ADT (with union-by-rank and path compression)
    Keep track of births __based on gray values__.  (i.e., "birth" is grayval... another idea: ordinal <-> gray)
    assumes this is for the foreground components (first ones born with high gray value)
    @notes
    Rather than taking in a type parameter,
    Makes each "set element" be a labelType ( defined in BettiVoxel.h ).
    (which is an unsigned-compatible type)
    See the Cormen-Leiserson-Rivest(-Stein) Algorithms textbook
    for rationale and merits of DSF w/ union-by-rank and path-compression
    Slightly modified it such that makeSet is no-ary: implicitly makes 
       a new set with a new (next) label.
    Take care of "number of components" (sets) here to decrease the burden of caller
    
    PROBABLY HAVE TO MODIFY IT FOR b2, the background components (param vs 2 versions?)
*/
// pre: labelType is an integral type (can safely be used as index to array/std::vector)

/** class dsf_for_tree_gen */
class dsf_for_tree_gen {
  std::vector<unsigned int> parent;
  std::vector<unsigned int> rank;
  std::vector<long long> born; // index is actual label, value is time in modificated gray-value
  // (make it a type parameter rather than hard-coding float? ###)
  std::vector<bool> updated; // "just updated" // leave it for now
  std::vector<unsigned> size;
  ///  std::vector<fcoord> centroid; // important only for active component
  // is there any way we can get away with storing 'size' and 'centroid' only in the historytree?
  std::vector<historynode*> HNodeP;
  // if representative, born does NOT pertain to label, but to entire component.
  // the two can be different since we are back to union-by-rank instead of union-by-age
  //  std::vector<labelType> atTime; // index is tick, value is affected component's label.  position 0 is bogus scratch.
  //  std::vector<bool> whatHappen; // true for births, false for deaths.  position 0 is bogus scratch.
  unsigned int nextLabel;
  unsigned int num_sets;
  unsigned int final_representative;
  //  labelType tick; // not really for labels, but why not make it compatible?
public:
  /** dsf_for_tree_gen constructor */
  dsf_for_tree_gen() {
    nextLabel = 2;
    final_representative = 0;
    // in the label image, '0' is reserved for background, non-set
    //  '1' came in for the bounded components of the complement...
    //  it does not affect what we do with 'foreground'.
    parent.push_back(0); // make 0 own parent: useless
    parent.push_back(1); // make 1 own parent: possibly useful?
    rank.push_back(0);
    rank.push_back(0);
    born.push_back(0); // whatever: birthtimes for the "null" component and the "boundary" component
    born.push_back(0); //   shall never be used, even if the components are.
    updated.push_back(false);
    updated.push_back(false); // does not matter what these two are, for now
    size.push_back(0);
    size.push_back(0);
    HNodeP.push_back(0);
    HNodeP.push_back(0); // well, in fact, we should have a node for the bgnds...
    num_sets = 0; // might want to make this 1 for 'b2'... hence, might want ctor to have param
  } // --dsf_for_tree_gen::dsf_for_tree_gen() [constructor]
  
  /** The DSF makeSet: implicitly "take" the next label */
  unsigned int makeSet(long long atgray) {
    parent.push_back(nextLabel);
    rank.push_back(0);
    born.push_back(atgray);
    size.push_back(1); // when just born, component size is 1
    updated.push_back(true);
    HNodeP.push_back(0); // null ptr. /might/ declare in history-related fnc below
    num_sets++;
    return nextLabel++; // return value then increment
  } // --dsf_for_tree_gen::makeSet()

  /** The DSF findSet (which does the path-compression) */
  unsigned int findSet(unsigned int x) {
    if(x != parent[x])
      parent[x] = findSet(parent[x]);
    return parent[x];
  } // --dsf_for_tree_gen::findSet()

  unsigned int increment_size(unsigned int x) {
    x = findSet(x);
    size[x]++;
    updated[x] = true;
		return x;
  }

  /** The dsf_for_tree_gen union BY RANK (again) operation
   *  union is a reserved word in C++, so we rename it to 'unite'
   */
  unsigned int unite(unsigned int x, unsigned int y) {
    x = findSet(x);
    y = findSet(y);
    if (x==y) return x; // DON'T decrement num_sets: do nothing at all
    // (in CLR's version of DSF, unite is called with representatives only)
    // now we have established that x and y are sitting on different trees
    unsigned int newRoot;
    if (rank[x] > rank[y]) {
      newRoot = parent[y] = x; // make x y's parent (and representative)
      // ### THIS HARDCODES b0-SPECIFIC BEHAVIOR! BAD!
      //      HNodeP[y] = HNodeP[x];
      size[x] += size[y]; // update component size
      if(born[x] < born[y]) {
				born[x] = born[y]; // new rep gets birthstamp of the older merging component
      }
    } else {
      newRoot = parent[x] = y; // make y x's parent (and representative)
      size[y] += size[x]; // update component size
      //      HNodeP[x] = HNodeP[y];
      if(born[y] < born[x]) {
				born[y] = born[x]; // new rep gets birthstamp of the older merging component
				// don't need "old" birthdates here for "persistence": that will be in the history tree.
      }
      if (rank[x] == rank[y]) {
        rank[y]++;
      }
    }
    num_sets--;
    updated[x] = updated [y] = true; // is there really a reason to specify dead component?...
    return newRoot;
  } // --dsf_for_tree_gen::unite()

  /** update history tree.  don't try to construct-as-pruned, just update.         *
   *  @params historytree reference and the current gray value                     *
   *  size changes are committed to ht as components change or die                 *
   *  at the end, all will be dead, except one--which will change at the last step *
   */
  void history_update(historytree& ht, long long curgray, std::map< labelType,std::vector <coord> > labelCoodMap,int tt) {
    int conta =0;
		int labelT =-1;
		int labelT2 =-1;
		int labelT3 =-1;		
		int labelT4 =-1;
		int labelT5 =-1;			
		for(unsigned int label = 1; label < nextLabel; label++) {       
			if(updated[label]) {
				if(bv::debug) printf("						Label:%d \n",(int)label);	
				if(label == parent[label]) { // if representative
					// N.B.: parent[label] is NOT NECESSARILY label's representative (==> HTree parent),
					//       but if it is equal to label, then it definitely is.
					// create new node in historytree (only for representatives!)
					historynode* nnp = new historynode(); // height 0 by default
					nnp->gray = curgray; // gray of this component in historytree--not necessarily birthgray
					nnp->size = size[label];
					ht.numnodes++;
					//nnp->listCoord = labelCoodMap[label];
					for(int i=0; i< labelCoodMap[label].size();i++){
                        coord mycoord = labelCoodMap[label][i];
						nnp->listCoord.push_back(labelCoodMap[label][i]);
						conta++;
					}							
					labelT =label;
					labelT2 =labelCoodMap[label].size();		
					if(HNodeP[label]) { // if already has current hnode, set hparent to the new hnode
						nnp -> addChild(HNodeP[label]);	
					} else {
					/// if(bv::debug) std::cout << "setting label " << label << 
					/// " to hnode (size: " << size[label] << ")" << std::endl;
					}
					HNodeP[label] = nnp; // regardless, set label's hnode ptr to the new hnode
					updated[label] = false; // set unupdated ONLY representatives!!!
				}
			}
		}
	// now set correct parents for nonrepresentatives that just changed
		for(unsigned int label = 1; label < nextLabel; label++) {
			// pre: parent will ALWAYS have greater size
			if(updated[label]) {
				if(HNodeP[label]) {
					unsigned int rep = findSet(label); // rep, not parent!
					if(0 == HNodeP[rep]) { // sanity check 
						std::cerr << "in dsf_for_tree_gen::history_update() rep has no HNode!!!" << std::endl;
						throw std::logic_error("in dsf_for_tree_gen::history_update() rep has no HNode!!!");
					}
					for(int i=0; i< labelCoodMap[label].size();i++){
                        coord mycoord = labelCoodMap[label][i];
						HNodeP[rep]->listCoord.push_back(labelCoodMap[label][i]);
						conta++;	
					}	
					labelT3 =label;	
					labelT4 =rep;
					labelT5 =labelCoodMap[label].size();			
					HNodeP[rep] -> addChild(HNodeP[label]);
					//	  std::cerr << "setting " << label << " as a child of " << parent[label] << std::endl;
				} else {
					for(int i=0; i< labelCoodMap[label].size();i++){
						HNodeP[findSet(label)]->listCoord.push_back(labelCoodMap[label][i]);
						conta++;
					}		
				// IT'S OK for HNodeP to be zero for updated non-representative labels
				// that means nothing was created for it in the tree.  naturally, if more than one label
				// with the same gray turns out to be neighbors, there is no need for an HNode for each, just the rep.
				//	  std::cerr << "For some reason HNodeP uninitialized for label " << label 
				//		    << " in dsf_for_tree_gen::history_update()" << std::endl;
				}	
				updated[label] = false;
			}
		
		}	
		if (conta!=tt){
				printf("ERROR the number of voxels add to VERTEX is wrong contado=%d previsto=%d curgray=%lld [%d, %d, %d, %d, %d]\n",conta,tt,curgray,labelT,labelT2,labelT3,labelT4,labelT5);		 
		}
    ///   std::cout << "HNodeP[2] is: " << HNodeP[2] << std::endl;
	//   if(bv::debug) std::cout << "--dsf_for_tree_gen::history_update()" << std::endl;
} // --dsf_for_tree_gen::history_update()

  /** get_num_sets returns ACTUAL number of components, none of the utilitarian/fake ones (?) */
  unsigned int get_num_sets() {
    return num_sets;
  } // --dsf_for_tree_gen::get_num_sets()

  unsigned int get_num_labels() { // 0 is a fake label, 1 is also fake for foreground components
    return nextLabel-2;
  }

  /** set root of historytree at the very end of everything! */
  void set_HT_root(historytree& ht) {
    // label 1 is problematic for now
    final_representative = 0;
    for(unsigned int label = 2; label < nextLabel; label++) {
      if(label == parent[label]) {
				if(final_representative != 0) {
	  			std::cerr << "Weirdness in set_HT_root(): More than one root!" << std::endl;
				}
				final_representative = label;
				if(bv::debug) {
	  			std::cout << "final representative label is " << label << std::endl;
	  			if(0 == HNodeP[final_representative]) {
	    			std::cout << "strangeness in set_HT_root(): HNodeP[final_representative] is null!" << std::endl;
	  			}
				}
				ht.hroot = HNodeP[final_representative];
				// no return for errorcheck
      }
    }
    if((bv::debug) && (final_representative == 0))
      std::cerr << "Error in set_HT_root(): no final representative!" << std::endl;
  }
}; // --class dsf_for_tree_gen
