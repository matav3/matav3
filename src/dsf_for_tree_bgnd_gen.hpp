/**
    @file dsf_for_tree_bgnd_gen.cpp
    @date 2004.08.16 by modifying dsf_for_tree_gen.cpp very slightly
    DSF for creation of historytree_bgnd of background voxels' components
    @author Deniz for both this and preds.
    @implements class dsf_for_tree_bgnd_gen
    @notes
    Rather than taking in a type parameter,
    Makes each "set element" be a labelType ( defined in BettiVoxel.h ).
    (which is an unsigned-compatible type)
    Slightly modified it such that makeSet is no-ary: implicitly makes 
       a new set with a new (next) label.
*/
// pre: labelType is an integral type (can safely be used as index to array/std::vector)
class dsf_for_tree_bgnd_gen {
  std::vector<labelType> parent;
  std::vector<unsigned char> rank;
  std::vector<unsigned char> born; // index is actual label, value is time in modificated gray-value
  // (make it a type parameter rather than hard-coding float? ###)
  std::vector<bool> updated; // "just updated" // leave it for now
  std::vector<unsigned> size;
  ///  std::vector<fcoord> centroid; // important only for active component
  // is there any way we can get away with storing 'size' and 'centroid' only in the historytree_bgnd?
  std::vector<historynode_bgnd*> HNodeP;
  // if representative, born does NOT pertain to label, but to entire component.
  // the two can be different since we are back to union-by-rank instead of union-by-age
  //  std::vector<labelType> atTime; // index is tick, value is affected component's label.  position 0 is bogus scratch.
  //  std::vector<bool> whatHappen; // true for births, false for deaths.  position 0 is bogus scratch.
  labelType nextLabel;
  labelType num_sets;
  labelType final_representative;
  //  labelType tick; // not really for labels, but why not make it compatible?
public:
  /** dsf_for_tree_bgnd_gen constructor */
  dsf_for_tree_bgnd_gen() {
    nextLabel = 2;
    final_representative = 0;
    // in the label image, '0' is reserved for background, non-set
    //  '1' came in for the bounded components of the complement...
    //  it does not affect what we do with 'foreground'.
    parent.push_back(0); // make 0 own parent: useless
    parent.push_back(1); // make 1 own parent: possibly useful?
    rank.push_back(0);
    rank.push_back(0);
    born.push_back(0); // birthtimes for the "null" component and the "boundary" component
    born.push_back(0); //   bgnd component IS born at gray 0.
    updated.push_back(false);
    updated.push_back(false);
    size.push_back(0);
    size.push_back(0); // make size of ONLY this 0: to be interpreted as infinite
    HNodeP.push_back(0);
    historynode_bgnd* nnp = new historynode_bgnd(true); // param true means "this is a boundary node"
    nnp->gray = 0;
    nnp->size = 0;
    HNodeP.push_back(nnp);
    num_sets = 1; // for boundary (only for bgnd)
  } // --dsf_for_tree_bgnd_gen::dsf_for_tree_bgnd_gen() [constructor]
  
  /** The DSF makeSet: implicitly "take" the next label */
  labelType makeSet(unsigned char atgray) {
    parent.push_back(nextLabel);
    rank.push_back(0);
    born.push_back(atgray);
    size.push_back(1); // when just born, component size is 1
    updated.push_back(true);
    HNodeP.push_back(0); // null ptr. /might/ declare in history-related fnc below
    num_sets++;
    return nextLabel++; // return value then increment
  } // --dsf_for_tree_bgnd_gen::makeSet()

  /** The DSF findSet (which does the path-compression) */
  labelType findSet(labelType x) {
    if(x != parent[x])
      parent[x] = findSet(parent[x]);
    return parent[x];
  } // --dsf_for_tree_bgnd_gen::findSet()

  void increment_size(labelType x) {
    x = findSet(x);
    size[x]++;
    updated[x] = true;
  }

  /** The dsf_for_tree_bgnd_gen union BY RANK (again) operation
   *  union is a reserved word in C++, so we rename it to 'unite'
   */
  labelType unite(labelType x, labelType y) {
    x = findSet(x);
    y = findSet(y);
    if (x==y) return x; // DON'T decrement num_sets: do nothing at all
    // (in CLR's version of DSF, unite is called with representatives only)
    // now we have established that x and y are sitting on different trees
    labelType newRoot;
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
  } // --dsf_for_tree_bgnd_gen::unite()

  /** update history tree.  don't try to construct-as-pruned, just update.         *
   *  @params historytree_bgnd reference and the current gray value                     *
   *  size changes are committed to ht as components change or die                 *
   *  at the end, all will be dead, except one--which will change at the last step *
   */
  void history_update(historytree_bgnd& ht, long long curgray) {
    if(bv::debug)
      std::cout << "in history_update() for gray " << static_cast<unsigned int> (curgray) << std::endl;
    //    for(each label in DSF)
    for(labelType label = 1; label < nextLabel; label++) {      
      if(updated[label]) {
				if(label == parent[label]) { // if representative
					// N.B.: parent[label] is NOT NECESSARILY label's representative (==> HTree parent),
					//       but if it is equal to label, then it definitely is.
					// create new node in historytree_bgnd (only for representatives!)
					historynode_bgnd* nnp;
					if(label == findSet(1)) nnp = new historynode_bgnd(true);
					else nnp = new historynode_bgnd(false);
					nnp->gray = curgray; // gray of this component in historytree_bgnd--not necessarily birthgray
					nnp->size = size[label];
					ht.numnodes++;
					if(HNodeP[label]) { // if already has current hnode, set hparent to the new hnode
						nnp -> addChild(HNodeP[label]);
					} else {
						
					}
					HNodeP[label] = nnp; // regardless, set label's hnode ptr to the new hnode
					updated[label] = false; // set unupdated ONLY representatives!!!
				}
      }
    }
    // now set correct parents for nonrepresentatives that just changed
    for(labelType label = 1; label < nextLabel; label++) {
      // pre: parent will ALWAYS have greater size
      if(updated[label]) {
				if(HNodeP[label]) {
					labelType rep = findSet(label); // rep, not parent!
					if(0 == HNodeP[rep]) { // sanity check
						std::cerr << "in dsf_for_tree_bgnd_gen::history_update() rep has no HNode!!!" << std::endl;
						throw std::logic_error("in dsf_for_tree_bgnd_gen::history_update() rep has no HNode!!!");
					}
					HNodeP[rep] -> addChild(HNodeP[label]);
					//	  std::cerr << "setting " << label << " as a child of " << parent[label] << std::endl;
				} else {
					// IT'S OK for HNodeP to be zero for updated non-representative labels
					// that means nothing was created for it in the tree.  naturally, if more than one label
					// with the same gray turns out to be neighbors, there is no need for an HNode for each, just the rep.
					//	  std::cerr << "For some reason HNodeP uninitialized for label " << label 
					//		    << " in dsf_for_tree_bgnd_gen::history_update()" << std::endl;
				}
				updated[label] = false;
      }
    }
  } // --dsf_for_tree_bgnd_gen::history_update()

  /** get_num_sets returns ACTUAL number of components, none of the utilitarian/fake ones (?) */
  labelType get_num_sets() {
    return num_sets;
  } // --dsf_for_tree_bgnd_gen::get_num_sets()

  labelType get_num_labels() { // 0 is a fake label.  for bgnd, 1 is not fake.
    return nextLabel-1;
  }

  /** set root of historytree_bgnd at the very end of everything! */
  void set_HT_root(historytree_bgnd& ht) {
    // label 1 is problematic for now
    std::cout << "in dsf_for_tree_bgnd_gen::set_HT_root()" << std::endl;
    final_representative = 0;
    for(labelType label = 1; label < nextLabel; label++) {
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
    if(bv::debug && (final_representative == 0))
      std::cerr << "Error in set_HT_root(): no final representative!" << std::endl;
  }
}; // --class dsf_for_tree_bgnd_gen
