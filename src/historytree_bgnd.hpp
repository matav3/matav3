#ifndef HISTORYNODEANDTREE_BGND
#define HISTORYNODEANDTREE_BGND
/** historytree_bgnd.cpp
    @date plagiarized from historytree.cpp on 2004.08.16
    @author deniz
    @implements
      struct historynode_bgnd
      struct historytree_bgnd
    decided to make everything public for easy access to fields, 
    yet provide some abstraction for ease of use.
    EVERY historynode_bgnd MUST BE ALLOCATED using 'new'.
    once root is properly set in historytree_bgnd, everything taken care of by historytree_bgnd's destructor.
*/

// don't think we need to write new versions of edgynode/edgytree for implementing Approach 2'

namespace bv {
  extern const bool debug;
}

struct historynode_bgnd {
  historynode_bgnd* hnext; // means: next SIBLING (in the order that was added)
  historynode_bgnd* hparent; // so that we never confuse it with DSF parent
  historynode_bgnd* hfirstchild; // ad-hoc queue (of children) implementation
  historynode_bgnd* hlastchild;
  unsigned size; // of components at particular point in time, in number of voxels... 
  // size can be interpreted differently (e.g., ignored) for unbounded components.
  unsigned char gray;  // the gray value ("time"): INCREASES towards the root for BGND.
  unsigned char height; // height of subtree rooted at this node
                        //  (absolute value of gray difference to the leaves)
  bool alive; // for pruning: if not alive, the subtree rooted at this node is dead
  bool unbounded; // node corresponds to the unbounded component at certain threshold

  historynode_bgnd(bool unbdd = false) {
    hnext = hparent = hfirstchild = hlastchild = 0;
    size = 0;
    height = 0; // "initialization saved my life."
    gray = 0; // why not...
    alive = true;
    unbounded = unbdd;
  }

  /** make (*chp) a child of (*this) [pre: it is not yet] */
  void addChild(historynode_bgnd* chp) {
    if(chp==0) { 
      std::cerr << "error: attempt to pass null pointer as child!" << std::endl;
      throw std::invalid_argument("tried to add null child");
      return;
    }
    if(chp->hparent) {
      std::cerr << "error in historytree_bgnd::addChild() : ch already has a parent!!!" << std::endl;
      throw std::logic_error("Parent already exists for this hnode!");
    }
    if(hlastchild) { // add new sibling
      hlastchild->hnext = chp;
      hlastchild = chp;
    } else {
      hfirstchild = hlastchild = chp;
    }
    chp->hparent = this; // set new child's parent
    // update height of this node based on new child
    unsigned char tmp = (chp->height) + (ABSDIF(gray, (chp->gray)));
    if (height < tmp) height = tmp;
  } // --historynode_bgnd::addChild()

  /** make (*par) the parent of (*this) [pre: (*this) has no parent yet] */
  void setParent(historynode_bgnd* par) {
    if(par==0) {
      std::cerr << "err in setParent(): null argument!" << std::endl;
      throw std::invalid_argument("in setParent(): null argument!");
    }
    if(hparent!=0) {
      std::cerr << "err in setParent(): Parent already exists for this hnode!" << std::endl;
      throw std::logic_error("Parent already exists for this hnode!");
    }
    par->addChild(this);
  } // --historynode_bgnd::setParent()

  /** deep destructor for (subtree rooted at) this node */
  ~historynode_bgnd() {
    ////    if(bv::debug) std::cerr << "historynode_bgnd destructor called!!!" << std::endl;
    historynode_bgnd* ch = hfirstchild;
    while(ch != 0) { // while still have children
      historynode_bgnd* nextch = ch->hnext; // remember next child
      delete ch; // get rid of first child (recursion)
      ch = nextch; // make next child first child
    }
  }
}; // --struct historynode_bgnd


/** struct historytree_bgnd: tree of historynodelist    *
 *  every node should be allocated with new (why?) */
struct historytree_bgnd {
  historynode_bgnd* hroot;
  unsigned numnodes;
  unsigned setalive; // number of nodes that are set alive: calculated by prune()
  unsigned char prunedbyheight;
  unsigned prunedbysize;
  /** historytree_bgnd "harmless" no-arg constructor */
  historytree_bgnd() {
    hroot = 0; // null.  HAVE to do this.
    numnodes = 1; // because we always have the first one...
    setalive = 0;
    prunedbyheight = 0;
    prunedbysize = 0;
  }
  
  /** precondition:  this historytree_bgnd is root-ed, with correct size & height @ nodes
   *  postcondition: sets alive field in top-nodes-to-be-pruned to false
   *                 i.e.: if i am dead, don't even bother with my children.
   *                 NB: does not result exactly in logical pruned tree--
   *                     if i am dead and my parent is alive, there might 
   *                     be logical nodes between us that are alive.
   */
  void prune(unsigned minsize=0, unsigned minheight=0) {
    std::stack<historynode_bgnd*> s;
    unsigned char mh = static_cast<unsigned char> (minheight);
    if(!s.empty()) {      // should NEVER happen
      std::cerr << "huge error: stack s non-empty at beginning"
		<< " of historytree_bgnd::prune()" << std::endl;
      throw std::runtime_error("stack nonempty at beginning of historytree_bgnd::prune()");
    }
    if(hroot == 0) {
      std::cerr << "in historytree_bgnd::prune() hroot is null!!!" << std::endl;
      throw std::logic_error("in historytree_bgnd::prune() hroot is null");
    }
    s.push(hroot);
    historynode_bgnd* np; // node pointer
    historynode_bgnd* chp; // child pointer
    setalive = 0;
    while(!(s.empty())) {
      np = s.top();
      s.pop();
      // Yung Kong wrote on 2004 August 15th for approach 2': "Do not prune unbounded branch by size but prune it by height"
      if(((np->unbounded) || (minsize <= (np->size))) 
	 && (mh <= (np->height))) { // not to be pruned
	(np->alive) = true; // set node alive
	setalive++;
	chp = (np->hfirstchild);
	while(chp != 0) { // and add ptrs to all children to the stack
	  s.push(chp);
	  chp = (chp->hnext);
	}
      } else { // otherwise, 'prune' and forget about children
	(np->alive) = false;
      }
    }
    if(setalive!=numnodes) // another hack so that unpruning the tree has no effect
      { 
	std::cout << setalive << " out of " << numnodes 
		  << " nodes of Bgnd tree alive after pruning (minsize=" << minsize 
		  << ", minheight=" << static_cast<unsigned int>(mh) << ")." << std::endl;
      }
    prunedbysize = minsize;
    //    std::cout << "prunedbysize set to " << prunedbysize 
    //	      << " in historytree_bgnd::prune()" << std::endl; 
    prunedbyheight = mh;
    //    std::cout << "prunedbyheight is " << static_cast<unsigned int>(prunedbyheight) 
    //	      << " in historytree_bgnd::prune()" << std::endl;
  } // --historytree_bgnd::prune()

  /** generate and return pointer to edgytree based on this (pruned) historytree_bgnd            *
   *  corresponds to "all_leaves == true" version of generate_pruned_simple_tree() sans bug *
   *  for bgnd: implement this last. ####  */
  edgytree* generate_pruned_edgy_tree() {
    // if all_leaves is set to true, then "extra" nodes will be 
    // added to leaves that are in the logical tree but were never in the "physical" tree
    // e.g., node has height 7 and 2 children at heights 3 and 4: if pruned at height 5, 
    // the edgytree will show this node to have 2 children at height 5.
    if(bv::debug) {
      std::cout << "prunedbyheight is " 
		<< static_cast<unsigned int> (prunedbyheight) 
		<< " in generate_pruned_edgy_tree()" << std::endl;
    }
    edgytree* etp = new edgytree();
    if(!(hroot->alive)) return etp;
    historynode_bgnd* hnp;
    historynode_bgnd* hchp;
    edgynode* enp;
    edgynode* echp;
    std::stack<historynode_bgnd*> hns; // history node stack
    std::stack<edgynode*> ens;  // edgy node stack
    hns.push(hroot);
    etp->eroot = new edgynode;
    etp->eroot->edgetoparent = 0; // root (and only root) has no parent.
    ens.push(etp->eroot);
    // deep copy ONLY 'alive' nodes
    while(!hns.empty()) {
      // sanity check:
      if(ens.empty()) { 
	std::cerr << "Non-empty hns but empty ens.  Returning null." << std::endl;
	return 0;
      }
      hnp = hns.top();
      hns.pop();
      enp = ens.top();
      ens.pop();
      if(hnp->hparent) { // only root has no parent, so 
// maybe doing this check every time is not the most efficient way
	(enp->edgetoparent) = ABSDIF(hnp->gray, hnp->hparent->gray);
      }
      if(hnp->alive) {
	hchp = (hnp->hfirstchild);
	edgynode* epp = 0; // ptr to edgy "predecessor" (left-sibling) node
	while(hchp != 0) { // and add ptrs to selected children to stacks
	  if(hchp->alive) {
	    hns.push(hchp);
	    enp->elastchild = echp = new edgynode;
	    ens.push(echp);
	    if(epp==0) {
	      enp->efirstchild = echp;
	    } else {
	      epp->enext = echp;
	    }
	    epp = echp;
	    echp = (echp->enext);
	  } else if((prunedbysize <= (hchp->size)) 
		    && (prunedbyheight < (ABSDIF((hchp->gray), (hnp->gray))))) {
	    std::cout << "ABSDIF(hchp->gray, hnp->gray) is " << 
	      static_cast<unsigned int>(ABSDIF((hchp->gray), (hnp->gray))) << std::endl;
// the latter MUST be a strict inequality: otherwise ENTIRE branch is pruned.
// if this has height that was prunedby, then it OUGHT TO be a logical leaf too.
// otherwise,
// children set dead must have height smaller than prunedbyheight,
// insert into physical tree the minimum-height acceptable ancestor in logical tree 
	    // (as long as child had enough size)
// all logical children of this new node will be really pruned (don't show at all)
// hence we do not use the stack
// (in fact, in this paradigm there is no more need to keep track of / check 'alive'...)
	    // construct new edgynode at height prunebyheight
	    enp->elastchild = echp = new edgynode;
	    /// assuming foreground component tree: nodes have strictly greater grays than their parents
	    ///	 echp->gray =  (hchp->gray) - prunedbyheight;
	    // not assuming anything--
	    echp->edgetoparent = (ABSDIF((hchp->gray) , (hnp->gray))) - prunedbyheight;
	    // perform sanity check ?
	    if(epp==0) {
	      enp->efirstchild = echp;
	    } else {
	      epp->enext = echp;
	    }
	  }
	  hchp = (hchp->hnext);
	}
	enp->elastchild = epp;
      }
    }
    // another sanity check:
    if(!ens.empty()) {
      std::cerr << "There are still nodes in the edgy node stack!  Returning null." << std::endl;
      return 0;
    }
    return etp;
  } // --generate_pruned_edgy_tree()
  

  /** recursively traverse subtree */
  void traverse_nonpruned_subtree(historynode_bgnd* hnp) {
    if(hnp==0) {
      std::cerr << "Error in traverse_nonpruned_subtree!  null pointer passed!" << std::endl;
      return;
    }
    if(hnp->alive) {
      std::cout << "Size: " << (hnp->size) << ", Gray: " 
		<< static_cast<unsigned>(hnp->gray) << ", Height: " 
		<< static_cast<unsigned>(hnp->height) << std::endl;
      historynode_bgnd *chp = hnp->hfirstchild;
      if(chp) {
	std::cout << "Children: [ " << std::endl;
	while(chp != 0) {
	  traverse_nonpruned_subtree(chp);
	  chp = chp->hnext;
	}
	std::cout << "]" << std::endl;
      }
    }
  } // --traverse_nonpruned_subtree()

  void traverse_nonpruned() {
    if(hroot == 0) {
      std::cerr << 
	"hroot is null in historytree_bgnd::traverse_nonpruned()"
		<< std::endl;
      return;
    }
    traverse_nonpruned_subtree(hroot);
  } // --traverse_nonpruned()

  /** is the tree empty?  yes if non-rooted AND numnodes is 0 */
  bool empty() {
    // watch it!    /// return (hroot==0);
    return ((numnodes==0) && (hroot==0));
  }

  /** deep destructor for historytree_bgnd */
  ~historytree_bgnd() {
    delete hroot;
  }
}; // --struct historytree_bgnd
#endif
