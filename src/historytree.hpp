#ifndef HISTORYNODEANDTREE
#define HISTORYNODEANDTREE
/** historytree.cpp
    2004.06.21--2004.07.01
    @author deniz
    @implements
      struct simplenode   # no-frill-whatsoever version of below
      struct simpletree   # no-frill-whatever version of below
      struct historynode
      struct historytree
    decided to make everything public for easy access to fields, 
    yet provide some abstraction for ease of use

    EVERY historynode MUST BE ALLOCATED using 'new'.
    once root is properly set in historytree, everything taken care of by historytree's destructor.0

    this file is ONLY FOREGROUND components.
    started historytree_bgnd.cpp for the tree of  background components of the complement
	
		starting in 09.09.2007
		By: Lucas M. Oliveira.
		Add new function for new Tree-Step simplification.


*/



namespace bv {
  extern const bool debug;
}


struct historynode {
    historynode* hnext; // means: next SIBLING (in the order that was added)
    historynode* hparent; // so that we never confuse it with DSF parent
    historynode* hfirstchild; // ad-hoc queue (of children) implementation
    historynode* hlastchild;
    historynode* unsimplified; // link for the corresopndent node in the unsimplified tree. I just wanna save memory and dont copy all the voxels in the simplified tree
    unsigned size; // of components at particular point in time, in number of voxels
    long long  gray;  // the gray value ("time"): DECREASES towards the root for fgnd.
    float gray_f;
    unsigned int height; // height of subtree rooted at this node
                //  (absolute value of gray difference to the leaves)
    bool alive; // for pruning: if not alive, the subtree rooted at this node is dead
    std::vector<coord> listCoord;


    // New atribute to version IV
    historynode** IPCD; //list of all the immediate proper critical descendants of c in T.

    // set by the procedure computeOutTree(). At the end of each iteration of the main loop Fast
    // Algorithm, if c is any vertex of the FHT that is to be output, then c.newIPCD will have been
    // set to a Slist of all the immediate proper critical descendants of c in the FHT that is to be output.
    historynode** newIPCD;

    //farthestDescendant is a  descendant v of T for which v.graylevel is maximal; in particular,
    //c.farthestDescendant is c itself if c is a leaf.
    historynode* farDescendant;

    //For any critical vertex c that is not a leaf, c.successor isthe immediate proper critical
    //descendant of c that lies on the path from c to c.farthestDescendant.
    historynode* successor;

    //The critical vertices c of T for which c.deltaLimit > delta are just the critical vertices of
    // an FHT that is produced by pruning (T,L,S) by branch length delta
    long long deltaLimit;

    //c.label is set by computeOutTree(), but only for critical vertices c that satisfy
    // c.deltaLimit > delta. For any such vertex c after execution of the "FOR d IN r.IPCD DO" loop
    // within the main loop of Algorithm 2:
    //      If c is a vertex of T that would have been deleted by the
    //      "WHILE e[j] <= delta1 DO" loop of the Inefficient Algorithm,
    //      then c.label is the value of e[j] at the iteration that
    //      would have deleted c.
    //
    //      If c would NOT have been deleted by the
    //      "WHILE e[j] <= delta1 DO" loop, then c.label > delta1.
    long long label;

    int sizeIPCD; // hold the next possition (or last after the simplificatio) to a IPCD be add in IPCD list
    int sizeNewIPCD	;
    int maximumK;
    int maxSizeNewIPCD;
    int maxSizeIPCD;
    int numberVoxel;
  historynode() {
    hnext = hparent = hfirstchild = hlastchild = farDescendant = successor =0;
    size = 0;
    height = 0; // "initialization saved my life."
    gray = 0; // why not...
		gray_f = 0.0;	
    alive = true;
        listCoord.clear();
		deltaLimit =0;
		label = 0;	
        maximumK =0;
        numberVoxel=0;
        sizeIPCD=sizeNewIPCD=maxSizeNewIPCD=maxSizeIPCD=0;
  }



/** MATAv3
make (*chp) a child of (*this) [pre: it is not yet]
*/
  void addChild(historynode* chp) {
    if(chp==0) { 
      std::cerr << "error: attempt to pass null pointer as child!" << std::endl;
      throw std::invalid_argument("tried to add null child");
      return;
    }
    if(chp->hparent) {
      std::cerr << "error in historytree::addChild() : ch already has a parent!!!" << std::endl;
      throw std::logic_error("Parent already exists for this hnode!");
    }
    if(hlastchild) { // add new sibling
      hlastchild->hnext = chp;
      hlastchild = chp;
    } else {
      hfirstchild = hlastchild = chp;
    }
    chp->hparent = this; // set new child's parent
    long long tmp = (chp->height) + ABSDIF(chp->gray,gray);
    if (height < tmp) height = tmp;  // update height of this node based on new child
		if(bv::debug) printf("Add child: pai height %d , filho height %d \n",height,chp->height);			
  } // --addChild()

/** MATAv3
make (*par) the parent of (*this) [pre: (*this) has no parent yet]
*/
  void setParent(historynode* par) {
    if(par==0) {	
      std::cerr << "err in setParent(): null argument!" << std::endl;
      throw std::invalid_argument("in setParent(): null argument!");
    }
    if(hparent!=0) {
      std::cerr << "err in setParent(): Parent already exists for this hnode!" << std::endl;
      throw std::logic_error("Parent already exists for this hnode!");
    }
    par->addChild(this);
  } // --setParent()

/** MATAv3
 deep destructor for (subtree rooted at) this node
*/
  ~historynode() {
    //if(bv::debug) std::cerr << "historynode destructor called!!!" << std::endl;
    historynode* ch = hfirstchild;
    while(ch != 0) { // while still have children
      historynode* nextch = ch->hnext; // remember next child
      delete ch; // get rid of first child (recursion)
      ch = nextch; // make next child first child
    }
  }
}; // --struct historynode


/** MATAv3
 struct historytree: tree of historynodelist    *
 every node should be allocated with new (why?) */
struct historytree {
  historynode* hroot;
  unsigned numnodes;
  unsigned setalive; // number of nodes that are set alive: calculated by prune()
  long long prunedbyheight;
  long long prunedbysize;
  unsigned properAncestor;
  long long omega;
	long long delta1;	
	int nBins;
	long long myPosInf;
	long long myNegInf;
    historynode** validAncestors;
  /** historytree "harmless" no-arg constructor */
  historytree() {
    hroot = 0; // null.  HAVE to do this.
    numnodes = 0;
    setalive = 0;
    prunedbyheight = 0;
    prunedbysize = 0;
    nBins =0;
    delta1=0;
    myNegInf=myPosInf=0;
  }
  

/** MATAv3
The two following function create a copy of the unsimplified tree (clone) to be used in the simplification.
*/
void cloneTree(historytree* newTree){
    if(bv::debug) printf( " \n cloneTree() \n");
    if (hroot==0){
        printf("YOU ARE SURE THAT YOU WANT CLONE AN EMPTY TREE?");
        return;
    }
    newTree->hroot = new historynode();
    newTree->nBins =nBins;
    newTree->myNegInf= myNegInf;
    newTree->myPosInf = myPosInf;
    cloneTreeRecursive(hroot,newTree->hroot);


}

void cloneTreeRecursive( historynode* original, historynode* clone ){
    if(bv::debug) printf( " \n cloneTreeRecursive %d \n",original->gray);
    //COPY NODE INFORMATION
    clone->size = original->size; // of components at particular point in time, in number of voxels
    clone->gray = original->gray;  // the gray value ("time"): DECREASES towards the root for fgnd.
    clone->gray_f= original->gray;
    clone->height=original->height; // height of subtree rooted at this node
    clone->alive =original->alive; // for pruning: if not alive, the subtree rooted at this node is dead
    clone->numberVoxel= original->numberVoxel;
    clone->unsimplified = original->unsimplified;

    //This is provisory.
    for(int i=0;i< original->listCoord.size();i++){
        clone->listCoord.push_back(original->listCoord[i]);
    }

    if(original->hfirstchild==0) {
      //end of the run. Only copy information.
      return;
    }

    historynode *chp = original->hfirstchild;
    historynode *fc_clone = new historynode();
    clone->addChild(fc_clone);
    cloneTreeRecursive(chp,fc_clone);//recursive first Child

    if (chp->hnext){
        historynode *nc = chp->hnext;
        historynode *nc_clone = new historynode();
        clone->addChild(nc_clone);
        cloneTreeRecursive(nc,nc_clone); //recursive nex Child
        historynode* pnc = nc->hnext;
        while (pnc!=0){
            nc = pnc->hnext;
            historynode *ot_clone = new historynode();
            clone->addChild(ot_clone);
            cloneTreeRecursive(pnc,ot_clone); //recursive nex Child
            pnc = nc;
        }
    }
}

  /** MATAV3
   * precondition:  this historytree is root-ed, with correct size & height @ nodes
   *  postcondition: sets alive field in top-nodes-to-be-pruned to false
   *                 i.e.: if i am dead, don't even bother with my children.
   *                 NB: does not result exactly in logical pruned tree--
   *                     if i am dead and my parent is alive, there might 
   *                     be logical nodes between us that are alive.
   */
  void prune(int minsize=0, int minheight=0) {
    std::stack<historynode*> s;
    long long mh = static_cast<long long> (minheight);
    if(!s.empty()) {      // should NEVER happen
      std::cerr << "huge error: stack s non-empty at beginning"
		<< " of historytree::prune()" << std::endl;
      throw std::runtime_error("stack nonempty at beginning of historytree::prune()");
    }
    if(hroot == 0) {
      std::cerr << "in historytree::prune() hroot is null!!!" << std::endl;
      throw std::logic_error("in historytree::prune() hroot is null");
    }
			
    std::cout << ":: Simplification by  Componente Size k="<< minsize <<" :: " << std::endl;
    s.push(hroot);
    historynode* np; // node pointer
    historynode* chp; // child pointer
    setalive = 0;
    while(!(s.empty())) {
			np = s.top();
			s.pop();
			if((minsize <= (np->size)) && (mh <= (np->height))) { // not to be pruned
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
			<< " nodes in Fgnd tree alive after pruning (minsize=" << minsize 
			<< ", minheight=" << static_cast<unsigned int>(mh) << ")." << std::endl;
		}
		prunedbysize = minsize;
		prunedbyheight = mh;
		//QString name = "pruneByCompSize_"+QString(QString::number(this->omega));
	} // --historytree::prune()

		

/** MATAV3
    generate and return pointer to edgytree based on this (pruned) historytree            *
*   corresponds to "all_leaves == true" version of generate_pruned_simple_tree() sans bug */
  edgytree* generate_pruned_edgy_tree() {
    // if all_leaves is set to true, then "extra" nodes will be 
    // added to leaves that are in the logical tree but were never in the "physical" tree
    // e.g., node has height 7 and 2 children at heights 3 and 4: if pruned at height 5, 
    // the edgytree will show this node to have 2 children at height 5.
    if(bv::debug) {
      std::cout << "generate_pruned_edgy_tree() h is  " 
		<< static_cast<unsigned int> (prunedbyheight) 
		<< " in generate_pruned_edgy_tree()" << std::endl;
    }
    edgytree* etp = new edgytree();
    if(!(hroot->alive)) return etp;
    historynode* hnp;
    historynode* hchp;
    edgynode* enp;
    edgynode* echp;
    std::stack<historynode*> hns; // history node stack
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
			//printf("{a} top = %lld & hparen= %lld \n",hnp->gray, hnp->hparent->gray);						
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
					} else 
						if((prunedbysize <= (hchp->size)) 
								&& (prunedbyheight < ABSDIF(hchp->gray, hnp->gray))) {
                                //std::cout << "ABSDIF(hchp->gray, hnp->gray) is " <<
                                //static_cast<unsigned int>(ABSDIF(hchp->height, hnp->height)) << std::endl;
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
                                // assuming foreground component tree: nodes have strictly greater grays than their parents
                                ///	 echp->gray =  (hchp->gray) - prunedbyheight;
                                echp->edgetoparent = (hchp->gray) - (hnp->gray) - prunedbyheight;
                                //printf("{b} top = %lld & hparen= %lld \n",hchp->gray, hnp->hparent->gray);
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
  

/** MATAv3
recursively traverse subtree
*/
  void traverse_nonpruned_subtree(historynode* hnp) {
    if(hnp==0) {
      std::cerr << "Error in traverse_nonpruned_subtree!  null pointer passed!" << std::endl;
      return;
    }
    if(hnp->alive) {
		historynode *chp = hnp->hfirstchild;
		if (chp){
			if (chp->hnext){
				historynode *aaa = chp->hnext;
				int cccc=2;
				historynode* pnc = aaa->hnext;
				while (pnc!=0){
					aaa = pnc->hnext;
					pnc = aaa;
					cccc++;
				}				
				//std::cout << "--->>> I'm critical with "<<aaa;
				//printf(" Size %d  , Gray %lld Height %d \n ",hnp->size,hnp->gray,hnp->height);				
			}
			
		}else{
			//std::cout << "--->>> I'm node.";
			//printf("Size %d  , Gray %lld Height %d \n ",hnp->size,hnp->gray,hnp->height);	
		}
		if(chp) {
			//std::cout << "Children: [ " << std::endl;
			while(chp != 0) {
				traverse_nonpruned_subtree(chp);
				chp = chp->hnext;
			}
			//std::cout << "]" << std::endl;
      	}
    }
  } // --traverse_nonpruned_subtree()


void traverse_nonpruned() {
    if(hroot == 0) {
      std::cerr << 	"hroot is null in historytree::traverse_nonpruned()"<< std::endl;
      return;
    }
    traverse_nonpruned_subtree(hroot);
  } // --traverse_nonpruned()

 /** is the tree empty?  yes if non-rooted AND numnodes is 0 */
  bool empty() {
    // watch it!    /// return (hroot==0);
    return ((numnodes==0) && (hroot==0));
  }

  /** deep destructor for historytree */
  ~historytree() {
    delete hroot;
  }


 /** MATAv3
 Copy all the voxel for of node to their parent by recursively traverse subtree.
 After call this function the nodes has all the voxel for their descendant..
 */
void updateAllVoxesTree(historynode* hnp) {
    //printf("updateAllVoxesTree(%d)  \n",hnp->gray);
    if (hnp->hfirstchild==0){
        copyVoxelToParent(hnp);
        //printf("[no]  [%lld,%d] \n",hnp->gray,hnp->listCoord.size());
    }else{
        historynode *chp = hnp->hfirstchild;
        //printf("[P]-[%lld,%d]  [FC]-[%lld,%d] \n",hnp->gray,hnp->listCoord.size(),chp->gray,chp->listCoord.size());
        updateAllVoxesTree(chp);
        while(chp->hnext!=0) {
            historynode *aux = chp->hnext;
            //printf("[P]-[%lld,%d]  [OC]-[%lld,%d] \n",hnp->gray,hnp->listCoord.size(),aux->gray,aux->listCoord.size());
            updateAllVoxesTree(aux);
            //copyVoxelToParent(aux);
            chp = aux;
        }
        //printf("[P]-[%lld,%d]  [AVO]-[%lld,%d] \n",hnp->gray,hnp->listCoord.size(),hnp->hparent->gray,hnp->hparent->listCoord.size());
        copyVoxelToParent(hnp);
        // copy  children
    }
} //

/** MATAv3
This function find the clone of a node in the unsimplified tree. The idea is have only one copy of the voxels in
the unsiplified tree. If you need the voxels of a node, get a copy in the  unsiplified tree.
*/

 historynode*  cloneInUnsTree(historynode *node, historytree *unsimpliedT,int k){

     std::queue<historynode*> s;
     s.push(hroot);
     std::queue<historynode*> sUns;
     sUns.push(unsimpliedT->hroot);
     while(!(s.empty())) {
         historynode* np_; // node pointer
         historynode* chp_; // child pointer
         np_ = s.front();
         s.pop();
         historynode* np_Uns; // node pointer
         historynode* chp_Uns; // child pointer
         np_Uns = sUns.front();
         sUns.pop();
         if(bv::debug) printf("\n<%d>((%d)) ",np_->gray,np_->numberVoxel);
         if (node==np_){
             if(bv::debug) printf("<%d>() np_->listCoord.size()=%d,  <%d> np_->listCoord.size()=%d \n",np_->gray, np_->listCoord.size(),np_Uns->gray, np_Uns->listCoord.size());
             return np_Uns;
         }

         if (np_->hfirstchild!=0){
             chp_ = (np_->hfirstchild);
             if(bv::debug) printf(" FC(%d)<%d> ",chp_->gray,chp_->numberVoxel);
             s.push(chp_);
             if(bv::debug) printf(" + ");
             historynode* pnc_ = chp_->hnext;

             //verify if the node in the unsimplified is smaller than k
             chp_Uns = (np_Uns->hfirstchild);
             if (chp_Uns==0){ if(bv::debug) printf(" AQUIII A \n");}
             historynode* pnc_Uns;
             while (chp_Uns->numberVoxel<k ){
                if (chp_Uns==0){ if(bv::debug) printf(" AQUIII B \n");}
                if(bv::debug) printf(" smallll* (%d)<%d>-------- \n",chp_Uns->gray,chp_Uns->numberVoxel);
                pnc_Uns = chp_Uns->hnext;
                chp_Uns = pnc_Uns;
             }
             sUns.push(chp_Uns);
             if(bv::debug) printf(" - ");
             pnc_Uns = chp_Uns->hnext;

             while(pnc_!= 0) {
                 s.push(pnc_);
                 if(bv::debug) printf(" NC(%d)<%d> ",pnc_->gray,pnc_->numberVoxel);
                 chp_ = pnc_->hnext;
                 pnc_ = chp_;
                 //verify if the node in the unsimplified is smaller than k
                 if (pnc_Uns==0){ if(bv::debug) printf(" AQUIII C \n");}
                 if (pnc_Uns->numberVoxel<k){
                    while (pnc_Uns->numberVoxel<k){
                        if(bv::debug) printf(" smallll (%d)<%d> ",pnc_Uns->gray,pnc_Uns->numberVoxel);
                        chp_Uns = pnc_Uns->hnext;
                        if (chp_Uns==0){ if(bv::debug) printf(" AQUIII D \n ");}
                        pnc_Uns = chp_Uns;
                    }
                    sUns.push(pnc_Uns);
                 }else{
                    sUns.push(pnc_Uns);
                    chp_Uns = pnc_Uns->hnext;
                    pnc_Uns = chp_Uns;
                 }

             }
         }
     }
     return node;
}

/** MATAv3  **/
void copyVoxelToParent(historynode* n){
    if (n->hparent!=0){
        historynode *pai= n->hparent;
        for(int i=0;i< n->listCoord.size();i++){
            pai->listCoord.push_back(n->listCoord[i]);
        }
    }
}

/** MATAv3
Copy all the # voxel for of node to their parent by recursively traverse subtree
*/
int updateNumberVoxelTree(historynode* hnp) {
     if (hnp->hfirstchild==0){
       hnp->numberVoxel=hnp->listCoord.size();
       return hnp->numberVoxel;
   }else{
       historynode *chp = hnp->hfirstchild;
       hnp->numberVoxel=hnp->numberVoxel + updateNumberVoxelTree(chp);
       while(chp->hnext!=0) {
           historynode *aux = chp->hnext;
           hnp->numberVoxel=hnp->numberVoxel+updateNumberVoxelTree(aux);
           chp = aux;
       }
       hnp->numberVoxel=hnp->numberVoxel+ hnp->listCoord.size();
       return hnp->numberVoxel;
   }
   return 0;
} //

/** MATAv3
Compute the ocuppance of the protein in the Tree. This function was building to be used in Paul Gottilieb LAB
*/

float computeOccuppanceTree(historynode* hnp) {
    if (hnp->hfirstchild==0){
       return hnp->listCoord.size() * hnp->gray;
   }else{
       historynode *chp = hnp->hfirstchild;
       float val =  computeOccuppanceTree(chp);
       while(chp->hnext!=0) {
           historynode *aux = chp->hnext;
           val=  val + computeOccuppanceTree(aux);
           chp = aux;
       }
       val = val + hnp->listCoord.size()*hnp->gray;
       return val;
   }
   return 0.0;
} //
/** MATAv3
*****************************************************************************
* This function find a h-th_proper_ancestor of v. Return an edgynode
* Ex: The 2-proper_ancestor is the grandparent of v. The 3-proper_ancestor
* is the grand_grandparent of v
******************************************************************************/

historynode* findProperAncestor(historynode* v){
	historynode* pa = v->hparent;
	if (v->hparent==0) return v;	
	while ( (v->gray-pa->gray) <= this->omega ) {		
		//chek if the ancentor have only one children and is a good-omega
		if ( (pa->gray >= v->gray - this->omega)  &&  (pa->hlastchild != pa->hfirstchild) ){			
			return pa;
		}else{	
			if (pa->hparent!=0) {
				historynode* chd = pa->hparent;					
				pa = chd;
			}else{
				return v;
			}
		}
	}
	return v;
}



/** MATAV3
****************************************************************************
* This function check is a node is critical vertex: has more than two children
* or is leave.
******************************************************************************/
bool isCriticalVertex(historynode* v){
	if (v->hfirstchild == 0 || v->hfirstchild != v->hlastchild){
		return true;
	}else{
		return false;
	}
}



/** MATAv3
****************************************************************************
* Funtion that traverse a tree and show information about the tree
* 
******************************************************************************/
void myTravesa() {
   std::stack<historynode*> s;
    printf("(***************************)-[myTravesa] <no> [pnc: [%lld,%d] \n",hroot->gray,hroot->size);	
	s.push(hroot);	 
    setalive = 0;
	int cont =0; 
	int contCoor =0;
	int contCv =0 ;
	int errorc=0;

    while(!(s.empty())) {
		historynode* np_; // node pointer
		historynode* chp_; // child pointer
		np_ = s.top();
		//s.pop();
		cont++;	
		contCoor+= np_->listCoord.size();	

		if (isCriticalVertex(np_)){
			contCv++;
			printf("-->> gray=%lld  [ ",np_->gray); 
		 }
		std::vector<coord>::iterator theIterator;
		for (theIterator = np_->listCoord.begin(); theIterator != np_->listCoord.end();
			theIterator++){	
			coord mycoord = (*theIterator);
		}					
				
		s.pop();

		if (np_->hfirstchild!=0){
			chp_ = (np_->hfirstchild);
			if (isCriticalVertex(np_)) printf(" (%lld) ",immediateCriticalVertice(chp_)->gray);	
			s.push(chp_);
			historynode* pnc_ = chp_->hnext;
			while(pnc_!= 0) {				
				s.push(pnc_);				
				if (isCriticalVertex(np_)) printf(" (%lld) ",immediateCriticalVertice(pnc_)->gray);
				chp_ = pnc_->hnext;
				pnc_ = chp_;				
			}
			
			if (isCriticalVertex(np_)) if (pnc_!= np_->hlastchild ) printf(" [%lld] ",immediateCriticalVertice(np_->hlastchild)->gray);
			
		}else{
			printf(" node "); 
		}
		if (isCriticalVertex(np_)) printf(" ]\n"); 
	}	
	std::cout << "\n(*)-[ this tree has " << cont <<
							" elements] and " << contCoor << "pixels= " << "  contCv="<< contCv<< " errorc=" << errorc << std::endl;
	
} 




/** MATAV3
* This function is used in the simplification STEP 3. Remove all node that
* are marked to be pruned (alive is false).
********************************************************************* *********/
 void ClearPrunedNode() {
    printf("::ClearPrunedNode: numnodes %d \n",numnodes);
    std::queue<historynode*> q;
    if(!q.empty()) {      // should NEVER happen
      std::cerr << "huge error: stack s non-empty at beginning"
				<< " of historytree::prune()" << std::endl;
      throw std::runtime_error("stack nonempty at beginning of historytree::prune()");
    }
		 
    if(hroot == 0) {
       std::cerr << "in historytree::prune() hroot is null!!!" << std::endl;
       throw std::logic_error("in historytree::prune() hroot is null");
    }
		int sonRem = 0; 
    q.push(hroot);
    historynode* np; // node pointer
    historynode* chp; // child pointer
   	int contGeral_v =0;
   	int contGeral_n =0;
    int nodesRem =0;
    while(!q.empty()) {
            np = q.front();
            q.pop();
            std::vector<coord>::iterator theIterator;
            for (theIterator = np->listCoord.begin(); theIterator != np->listCoord.end();theIterator++){
                coord mycoord = (*theIterator);
            }
            //The node alive is false. Need be removed
            if (!np->alive){
                int paiA = np->hparent->listCoord.size();
                //printf("{%lld,%lld} (p=%d,f=%d) -- [p=%d,f=%d]",np->hparent->gray,np->gray, np->hparent->listCoord.size(),np->listCoord.size(),np->hparent->numberVoxel,np->numberVoxel);
                sonRem = changeVoxelToParent_down(np);
                //printf("{%lld,%lld}    [p=%d,f=%d]  -- [p=%d,f=%d] \n",np->hparent->gray,np->gray, np->hparent->listCoord.size(),np->listCoord.size(),np->hparent->numberVoxel,np->numberVoxel);
                if (sonRem!=(np->hparent->listCoord.size()-paiA)){
                    printf(" ERRO paiA=%d paiD=%d rem=%d \n",paiA,np->hparent->listCoord.size(),sonRem);
                }
                removeFromTree(np);
                numnodes--;
                contGeral_v= contGeral_v+ sonRem;
            }else{
                contGeral_n= contGeral_n+1;
                contGeral_v= contGeral_v+ np->listCoord.size();
                if (np->hfirstchild!=0){
                    chp = np->hfirstchild;
                    q.push(chp);
                    while (chp->hnext!=0){
                        historynode* nex = chp->hnext;
                        q.push(nex);
                        chp = nex;
                    }
                }
            }
        }
		 
        int nn = 0;
        int cont =0;
        q.push(hroot);
        while(!q.empty()) {
            np = q.front();
            q.pop();
            nn++;
            cont+=np->listCoord.size();
            if (np->hfirstchild!=0){
                chp = np->hfirstchild;
                q.push(chp);
                while (chp->hnext!=0){
                    q.push(chp->hnext);
                    chp = chp->hnext;
                }
            }
        }
        numnodes = nn;
		//printf("numnodes %d   #voxels=%d contGeral_v=%d contGeral_n=%d\n",numnodes,cont,contGeral_v,contGeral_n);
} // --historytree::ClearPrunedNode()

/******************************************************************************
* This function remove a node. In this case only the element is removed. The
* children will be child of the parent of a node removed.
******************************************************************************/
void removeNode(historynode* v){
	if (v==hroot)
		return;		
	historynode* pai = v->hparent;
	numnodes--;	
	//case 1: is firschild
	if (pai->hfirstchild == v){
		//case 1.1: is leaf			
		if (v->hfirstchild ==0){
			if(v->hnext!=0){
				pai->hfirstchild = v->hnext;
			}else{
				pai->hfirstchild = 0;
				pai->hlastchild = 0;	
			}
		//case 1.2: have one or more children
		}else{
			historynode* frt = v->hfirstchild;
			pai->hfirstchild = frt;
			frt->hparent = pai;	
			// link with another children
			historynode* aux = frt;
			while (aux->hnext!=0){
				aux->hparent = pai;	
				aux = aux->hnext;	
			}
			if (v->hnext != 0){ 
				aux->hnext = v->hnext;	
			}else{
				aux->hnext = 0;
				pai->hlastchild= aux;
			}
			v->hfirstchild=0;
		}
			
	//case 2: is a simple children
	}else{
		//get  before				
		historynode* bef = pai->hfirstchild;
		while(bef->hnext !=v){
			bef = bef->hnext;
			if (bef==0)	
						printf("---- ----> [error ] \n");
		}			

		//case 2.1: is leaf				
		if (v->hfirstchild ==0){
			bef->hnext = v->hnext;
			if (pai->hlastchild==v){
				pai->hlastchild = bef;
			}
		//case 2.2: have one or more children
		}else{
			historynode* frt = v->hfirstchild;
			bef->hnext = frt;
			frt->hparent = pai;	
			// link with another children
			historynode* aux = frt;
			while (aux->hnext!=0){
				aux->hparent = pai;	
				aux = aux->hnext;	
			}
			aux->hnext = v->hnext;	
			if (pai->hlastchild==v){
				pai->hlastchild = aux;
			}
			v->hfirstchild=0;	
		}	
	}
	delete(v); // calls deep deep destructor	
}



/** MATAVA3
* This function is used to remove the node of tree
******************************************************************************/
void removeFromTree(historynode *v) {
    if (v==hroot) {
      std::cerr << "in removeFromTree(): ATTEMPT TO REMOVE ROOT FROM TREE!" << std::endl;
      return;
    }
    historynode* p = (v->hparent);	
    if(p != 0) {	
      historynode* c = (p->hfirstchild);
      if(c==0) {
				return;
      }
      if(v == c) { // v is the first child
				if(v == (p->hlastchild)) { // v is the only child. Rem it.
                //if(bv::debug) std::cerr << "removing some only child from tree" << std::endl;
	  			p->hfirstchild = 0;
	  			p->hlastchild = 0;
				} else { // set the v's next to firschild 
	  			p->hfirstchild = v->hnext;
				}
      } else { // there is at least one child that is not v
				while((c->hnext) != v) {
	  			c = (c->hnext);
	  			if(c==0) {
	    			if(bv::debug) std::cerr << "fatal error: parent does not know about v!" << std::endl;
	    			break;
	    			//return;
	  			}
				}
				if(c != 0) {
			 		if((p->hlastchild) == v) {
						(p->hlastchild) = c;
						(c->hnext) = 0;
			  	} else {
						(c->hnext) = (v->hnext);
			  	}
				}
      }
			v->hparent = 0;
			v->hnext = 0;
			v->hfirstchild = 0; // these are NECESSARY
			v->hlastchild = 0;
			delete(v); // calls deep deep destructor
    }else{
      std::cerr << "parent is null in removeFromTree() ??" << std::endl;
    }
} 



/*
The function change all voxels associated with one vextex to your parent
*/
  int changeVoxelToParent_down(historynode* n) {
    if(n==0) {
      std::cerr << "Error in traverse_nonpruned_subtree!  null pointer passed!" << std::endl;
      return;
    }    
    int ret =0;
    if(n->hfirstchild!=0){
            historynode* ch = n->hfirstchild;
            ret =	n->listCoord.size();
            ret += changeVoxelToParent_down(ch);
            while(ch->hnext!=0){
                historynode* aux = ch->hnext;
                ret+= changeVoxelToParent_down(aux);
                ch =aux;
            }
            historynode *pai= n->hparent;
            for(int i=0;i< n->listCoord.size();i++){
                pai->listCoord.push_back(n->listCoord[i]);
            }
            n->listCoord.clear();
            return ret;
    }else{
            historynode *pai= n->hparent;
            for(int i=0;i< n->listCoord.size();i++){
                pai->listCoord.push_back(n->listCoord[i]);
                ret++;
            }
            n->listCoord.clear();
            return ret;
    }
 }

/** MATAV3
This function retorn the immediate critical vertice of node
*/
historynode* immediateCriticalVertice(historynode *node){
	if (isCriticalVertex(node)){
			return node;
	}else{
		historynode* chp = node->hfirstchild;	
		while(!isCriticalVertex(chp)){	
			chp = chp->hfirstchild;		
		}
		return chp;
	}			
}

/******************************************************************************
This function find all absolute difference between two critical vertex u and v,
such that u is ancestor of v, u is not leaf, and L(u) - L(v) <= delta
******************************************************************************/
std::set<long long> differentLevelCVset(long long delta){
    std::queue<historynode*> qGeral;
		std::queue<historynode*> qDiff;

    qGeral.push(hroot);
    historynode* np; // node pointer
    historynode* chp; // child pointer
   	set<long long> myset;
		int count =0;
    while(!qGeral.empty()) {
			np = qGeral.front();
			qGeral.pop();					
			//add all non leaf children of np	
			if (np->hfirstchild!=0){	
				chp = immediateCriticalVertice(np->hfirstchild);
				if (chp->hfirstchild !=0){
						qGeral.push(chp);
						if (chp->gray-np->gray<=delta) qDiff.push(chp);
				}
					
				chp = np->hfirstchild;					
				while(chp->hnext!=0){
					chp = chp->hnext;
					historynode*  aux = immediateCriticalVertice(chp);	
					if (aux->hfirstchild !=0){
						qGeral.push(aux);
						if (aux->gray-np->gray<=delta) qDiff.push(aux);
					}
				}
			}
			historynode* ndif_p;
			historynode* ndif_c;	
			while(!qDiff.empty()) {
				ndif_p =qDiff.front();
				qDiff.pop();					
				myset.insert(ndif_p->gray - np->gray);
				
				if (ndif_p->hfirstchild!=0){	
					ndif_c = immediateCriticalVertice(ndif_p->hfirstchild);
					if ( (ndif_c->hfirstchild !=0) && (ndif_c->gray-np->gray<=delta) ){
						qDiff.push(ndif_c);
					}
				}
					
				ndif_c = ndif_p->hfirstchild;					
				while(ndif_c->hnext!=0){
					ndif_c = ndif_c->hnext;
					historynode*  aux = immediateCriticalVertice(ndif_c);	
					if ( (aux->hfirstchild !=0)&&(aux->gray-np->gray<=delta) ){
						qDiff.push(aux);
					}
				}
			}
		}
		//printf("\n # of cv :%d \n",myset.size());
	
		

		//printf("\n -------------- \n",myset.size());
		set<long long>::iterator it;
		int ct =0;
		for (it=myset.begin() ; it != myset.end(); it++){		
			long long a = static_cast<long long> (*it);
			//printf("[%d] %lld \n",ct++,a);	
		}
		return myset;
}


/******************************************************************************
This function find all absolute difference between two critical vertex u and v,
such that u is ancestor of v and u is not leaf.
******************************************************************************/
std::set<long long> differentLevelCVset(){
     std::queue<historynode*> qGeral;
		std::queue<historynode*> qDiff;

    qGeral.push(hroot);
    historynode* np; // node pointer
    historynode* chp; // child pointer
   	set<long long> myset;
		int count =0;
    while(!qGeral.empty()) {
			np = qGeral.front();
			qGeral.pop();					
			//add all non leaf children of np	
			if (np->hfirstchild!=0){	
				chp = immediateCriticalVertice(np->hfirstchild);
				if (chp->hfirstchild !=0){
					qGeral.push(chp);
					qDiff.push(chp);
				}
					
				chp = np->hfirstchild;					
				while(chp->hnext!=0){
					chp = chp->hnext;
					historynode*  aux = immediateCriticalVertice(chp);	
					if (aux->hfirstchild !=0){
						qGeral.push(aux);
						qDiff.push(aux);
					}
				}
			}
			historynode* ndif_p;
			historynode* ndif_c;	
			while(!qDiff.empty()) {
				ndif_p =qDiff.front();
				qDiff.pop();					
				myset.insert(ndif_p->gray - np->gray);
				
				if (ndif_p->hfirstchild!=0){	
					ndif_c = immediateCriticalVertice(ndif_p->hfirstchild);
					if ( (ndif_c->hfirstchild !=0) ){
						qDiff.push(ndif_c);
					}
				}
					
				ndif_c = ndif_p->hfirstchild;					
				while(ndif_c->hnext!=0){
					ndif_c = ndif_c->hnext;
					historynode*  aux = immediateCriticalVertice(ndif_c);	
					if ( (aux->hfirstchild !=0) ){
						qDiff.push(aux);
					}
				}
			}
		}
		//printf("\n # of cv :%d \n",myset.size());
	
		

		//printf("\n -------------- \n",myset.size());
		set<long long>::iterator it;
		int ct =0;
		for (it=myset.begin() ; it != myset.end(); it++){		
			long long a = static_cast<long long> (*it);
			//printf("[%d] %lld \n",ct++,a);	
		}
		return myset;
}

void saveImageTree(QString name){
        edgytree* etp = this->generate_pruned_edgy_tree();
		etp->flatten(); // once again, does not remove 'edge' from root down
		etp->setparents();
		QString NewickFileName = name+".new";
		QString XBMFileName =  name+".xbm";
		QString jpgName = name+".jpg";
		if(true) {
			std::cout  << "obtaining new Newick representation from edgytree and writing it to " << NewickFileName.toAscii().constData() << std::endl;
			std::ofstream NewickFile (NewickFileName.toAscii().constData(), std::ios::out); // replace if exists!
			if (NewickFile.is_open()) {
				etp->NEXUS_out(NewickFile); // edgytree::NEXUS_out() writes in the simpler Newick format rather than NEXUS
				NewickFile.close();
			} else {
				cout << "Error opening " << NewickFileName.toAscii().constData() << " for writing" << endl;
				return;
			}
		}
		
		QString command = "drawgram " +  // botanist warns user as first thing if no drawgram in path
		NewickFileName + " " + // input file
		XBMFileName;   // output file
		system(command.toAscii().constData());
		command ="";
		command = "convert -flip " +  // botanist warns user as first thing if no drawgram in path (-rotate -180)
		XBMFileName + " " + // input file
		jpgName;   // output file
		system(command.toAscii().constData());		
}			














/* #############################################################################

THIS FOLLING FUNCTION ARE NOT BEEN USED

#############################################################################*/






/******************************************************************************
* The rescursive funtion propose by prof YUNG to collpse class using a delta.
The function receive the root of the tree and retorn the pruned tree
*
******************************************************************************/
void collapseClasses(historynode* r, historynode* v, long long delta){
    std::vector<historynode*> immediateCV;
    immediateCV = getImmediateCV(v);
    for (int i=0;i<	immediateCV.size();i++){
        historynode* icv = immediateCV[i];
        //remove from T all proper descendant of v that are proper ancestors of icv
        historynode* aux = icv->hparent;
        while(aux!=v){
            for(int i=0;i<aux->listCoord.size();i++){
                aux->hparent->listCoord.push_back(aux->listCoord[i]);
            }
            historynode* paiaux = aux->hparent;
            aux->listCoord.clear();
            removeNode(aux);
            aux = paiaux;
        }
        if ( icv->gray <= v->gray+delta){
            collapseClasses(r,icv,delta);
            for(int i=0;i<icv->listCoord.size();i++){
                icv->hparent->listCoord.push_back(icv->listCoord[i]);
            }
            removeNode(icv);
        }else{
            icv->hparent = r;
            collapseClasses(icv,icv,delta);
        }
    }
}

/******************************************************************************
* The function receive a node e retorn the immediate critical vertex
*
******************************************************************************/
std::vector<historynode*> getImmediateCV(historynode* u){
  std::vector<historynode*> listCV;
    listCV.clear();

    if (u->hfirstchild==0)
        return listCV;

    // firschild immediate CV
    historynode* fc = u->hfirstchild;
    while( !isCriticalVertex(fc) ){
        fc = fc->hfirstchild;
    }
    listCV.push_back(fc);

    // other child immediate CV
    historynode* aux = u->hfirstchild;
    while(aux->hnext!=0){
        historynode* im = aux->hnext;
        while( !isCriticalVertex(im) ){
            im  = im ->hfirstchild;
        }
        listCV.push_back(im);
        aux = aux->hnext;
    }

    return listCV;
}


/******************************************************************************
* Simplification by  Collapsing Delta Equivalent Class
The step 3 of prune simplification. This is the non-incremental version 2.
The give the parameter delta to prune the tree.
*
******************************************************************************/
void pruneByDeltaStepC_recursive(int delta) {
    std::cout << ":: Simplification by  Collapsing Delta Equivalent Class - Recursive:: " << std::endl;
    collapseClasses(hroot,hroot,static_cast<long long> (delta));
}



/*******************************************************************************
* FUNCTION FOR THE TREE STEP SIMPLIFICATION  --------VERSION 2.--------
parameter:
omega: 	is used to collapsed class
allDiff: is used to say if will be used all different graylevel (true) or not (false)
printFht: control if the intermediate tree will be pruned.
increment: if printFht is true interval that the tree will be pruned. All tree with
                        less than 100 nodel will be pruned.
namePrint: where the tree will be salved.
*******************************************************************************/
void pruneByDeltaStepC_recursive(int omega, bool allDiff,bool printFht, int incremento, QString namePrint, int minNode) {
    std::cout << ":: Simplification by Interative Collapsing Delta Equivalent Class - Recursive:: " << std::endl;

    set<long long> myset;

    if (allDiff){
        myset = differentLevelCVset();
    }else{
        myset = differentLevelCVset((long long)omega);
    }
    set<long long>::iterator it;
    int tot = myset.size();
    int parcial = 0;
    int acumulado =0;
    //int numnodes_ant = numnodes;
    //printf("STEP 3 = 0.00 %  \n");
    for (it=myset.begin() ; it != myset.end(); it++){
        this->omega = static_cast<long long> (*it);
        parcial++;
        float percent = (float)(((float)parcial/tot)*100);
    //printf("STEP 3 = %.2f % [%d/%d] -(%lld) \n",percent,parcial,tot,this->omega);
        collapseClasses(hroot,hroot,this->omega);

        if(  ((numnodes<100) || (this->omega > acumulado)) ){// && (numnodes_ant>numnodes) ){
                //numnodes_ant = numnodes;
                acumulado = acumulado+incremento;
                QString name = namePrint+"_n"+QString(QString::number(numnodes))+"_e"+QString(QString::number(this->omega));
                if (printFht) saveImageTree(name);
        }
        if (numnodes <= minNode)
            break;
    }
}


/*******************************************************************************
* FUNCTION FOR THE TREE STEP SIMPLIFICATION  -----  VERSION 3 --------
parameter:
omega: 	is used to collapsed class
allDiff: is used to say if will be used all different graylevel (true) or not (false)
printFht: control if the intermediate tree will be pruned.
increment: if printFht is true interval that the tree will be pruned. All tree with
                        less than 100 nodel will be pruned.
namePrint: where the tree will be salved.
*******************************************************************************/

void pruneByDeltaStepBC_recursive(int omega, bool allDiff,bool printFht, int incremento, QString namePrint, int minNode) {
    std::cout << ":: Simplification by Length Branch/ Interative Collapsing Delta Equivalent Class  - Recursive:: " << std::endl;

    set<long long> myset;

    if (true){
        myset = differentLevelCVset();
    }else{
        myset = differentLevelCVset((long long)omega);
    }
    set<long long>::iterator it;
    int tot = myset.size();
    int parcial = 0;
    int step_p = tot/50;
    int step = tot/50;
    int acumulado =0;
    int numnodes_ant = numnodes;
    for (it=myset.begin() ; it != myset.end(); it++){
        this->omega = static_cast<long long> (*it);
        parcial++;
        float percent = (float)(((float)parcial/tot)*100);
    if (parcial>step){
            //printf("STEP 3 = %.2f  [%d/%d] -(%lld)\n",percent,parcial,tot,this->omega);
            step= step +step_p;
        }
        pruneByDeltaStepB(this->omega);
        //std::cout << ":: collapseClasses:: " << this->omega<< std::endl;
        collapseClasses(hroot,hroot,this->omega);
        if(  ((numnodes<100) || (this->omega>acumulado)) && (numnodes_ant>numnodes) ){
                numnodes_ant = numnodes;
                acumulado = acumulado+incremento;
                QString name = namePrint+"_n"+QString(QString::number(numnodes))+"_e"+QString(QString::number(this->omega));
                if (printFht) saveImageTree(name);
        }
        if (numnodes <= minNode)
            break;
    }
}



/******************************************************************************
* This function prune the tree by Collapsing Delta Equivalent Class. The input
*   is the tree building by pruneByDeltaStepB (step 2).
*		This funciont will get each critical vertex v and mark alive=false
*		to all equivalent class that have root in v.
******************************************************************************/
 void pruneByDeltaStepC(int omega) {
    std::queue<historynode*> q;
		//myTravesa();  
		flattenClass();
		//myTravesa();  
		this->omega = static_cast<long long> (omega);  
		std::cout << ":: Simplification by  Collapsing Delta Equivalent Class :: " << std::endl; 
    int leaff =1;
    if(!q.empty()) {      // should NEVER happen
      std::cerr << "huge error: stack s non-empty at beginning"
				<< " of historytree::prune()" << std::endl;
      throw std::runtime_error("stack nonempty at beginning of historytree::prune()");
    }
    if(hroot == 0) {
       std::cerr << "in historytree::prune() hroot is null!!!" << std::endl;
      throw std::logic_error("in historytree::prune() hroot is null");
    }

		if(bv::debug)	printf("(*)-[root]  [ %lld, %d, %d] \n",  hroot->gray, hroot->size,hroot->height);		 
		
        //PARTE 1: Mark all (true) node that will be removed of tree
    q.push(hroot);
    historynode* np; // node pointer
    //historynode* chp; // child pointer
   	
    while(!(q.empty())) { 
			np = q.front();
			q.pop();
			if (np->hfirstchild==0){
				if(bv::debug)	printf("(*)-[leaf]  [ %lld, %d, %d] \n",np->gray, np->size,np->height);		
				leaff++;
				np->alive = false;
			}else{
				//the vertex is not critical. It have only one children.
				if (!isCriticalVertex(np) ){
					if(bv::debug)	printf("(*)-[not is critical] [ %lld, %d, %d] \n",np->gray, np->size,np->height);	
					historynode* aux = np->hfirstchild;
					q.push(aux);			
				}else{
					if(bv::debug)	printf("(*)-[critical] [ %lld, %d, %d] \n",np->gray, np->size,np->height);	
					if (np->alive){					
						find_R_deltaEquivalenteRelation(np);
					}
					historynode* aux = np->hfirstchild;
					while (aux!=0 ){
						q.push(aux);
						aux = aux->hnext;
					}								
				}
			}
		}		
		//PARTE 2: Mark all node that will be removed of tree 			
	  
		if(bv::debug)	printf("PREPARANDO saida pra YUNG\n\n\n\n\n");
		
		this->makeClass();
		
		if(bv::debug)	 printf("IMPRIMINDO saida pra YUNG\n\n\n\n\n") ;
		 
		 
} // --historytree::pruneByOmega()

/******************************************************************************
*  VERSION 2: now delta incremental
* This function prune the tree by Collapsing Incremental Delta Equivalent Class. 
*		The input is the tree building by pruneByDeltaStepB (step 2).
*		This funciont will get each critical vertex v and mark alive=false
*		to all equivalent class that have root in v.
******************************************************************************/
/*First sort the elements of the set of numbers
  { L(x) - L(y) | x and y are critical vertices of T
                  such that x is a proper descendant of y
                  and  L(x) - L(y) <= delta }
into a strictly increasing sequence e[1] < e[2] < ... < e[n].

Then execute the following loop:
    FOR i = 1 TO n DO
    END FOR;
*/

 void pruneByDeltaStepC_interative(int omega,bool allDiff) {
		std::cout << ":: Simplification by  Interative Collapsing Delta Equivalent Class :: (delta=" << this->omega  << ")"<< std::endl; 
		set<long long> myset;
		
		if (allDiff){
		 	myset = differentLevelCVset();	
		}else{
			myset = differentLevelCVset((long long)omega);		 
		}
		//myTravesa(); 
		flattenClass(); 
		//myTravesa(); 
		set<long long>::iterator it;
		int tot = myset.size();
		int parcial = 0;
		myset.insert((long long)omega);
		 
		for (it=myset.begin() ; it != myset.end(); it++){		
			this->omega = static_cast<long long> (*it);
			parcial++;	
			float percent = (float)(((float)parcial/tot)*100);
    	//printf("STEP 3 = %.2f  [%d/%d] -(%lld)\n",percent,parcial,tot,this->omega);					
			std::queue<historynode*> q; 
			//unsigned lop =  static_cast<unsigned> (i);
			int leaff =1;
			if(!q.empty()) {      // should NEVER happen
				std::cerr << "huge error: stack s non-empty at beginning"
					<< " of historytree::prune()" << std::endl;
				throw std::runtime_error("stack nonempty at beginning of historytree::prune()");
			}
			if(hroot == 0) {
				 std::cerr << "in historytree::prune() hroot is null!!!" << std::endl;
				throw std::logic_error("in historytree::prune() hroot is null");
			}
			if(bv::debug)	printf("(*)-[root] [%lld, %d, %d] \n",hroot->gray, hroot->size,hroot->height);
			//if(bv::debug)	std::cout << "(*)-[root]  [" << static_cast<int> (hroot->gray ) << "," << static_cast<int> (hroot->size ) << "," << static_cast<int> (hroot->height )<<"]" <<std::endl;	
				
			//PARTE 1: Mark all (true) node that will be removed of tree 
			q.push(hroot);
			historynode* np; // node pointer
			//historynode* chp; // child pointer
			
			while(!(q.empty())) { 
				np = q.front();
				q.pop();
				if (np->hfirstchild==0){
					if(bv::debug)	printf("(*)-[leaf] [ %lld, %d, %d] \n",np->gray, np->size,np->height);
					leaff++;
					np->alive = false;
				}else{
					//the vertex is not critical. It have only one children.
					if (!isCriticalVertex(np) ){
						if(bv::debug)	printf("(*)-[not is critical] [ %lld, %d, %d] \n",np->gray, np->size,np->height);							
						historynode* aux = np->hfirstchild;
						q.push(aux);			
					}else{
						//the vertex is critical. If alive=true says that this vertex need be evaluate.
						if(bv::debug)	printf("(*)-[critical] [ %lld, %d, %d] \n",np->gray, np->size,np->height);															
						if (np->alive){					
							find_R_deltaEquivalenteRelation(np);
						}
						historynode* aux = np->hfirstchild;
						while (aux!=0 ){
							q.push(aux);
							aux = aux->hnext;
						}								
					}
				}
			}		
			//PARTE 2: Mark all node that will be removed of tree 			
			
			if(bv::debug)	printf("PREPARANDO saida pra YUNG\n\n\n\n\n");
	
			this->makeClass();
			
			if(bv::debug)	printf("IMPRIMINDO saida pra YUNG\n\n\n\n\n") ;
				edgytree* etp = this->generate_pruned_edgy_tree();								
				etp->flatten(); // once again, does not remove 'edge' from root down
				etp->setparents();
									//~ QString NewickFileName = QString(QString::number(this->omega))+".new";
									//~ QString XBMFileName =  QString(QString::number(this->omega))+".xbm";
									//~ QString jpgName = QString(QString::number(this->omega))+".jpg";
									//~ //etp->testando(etp->eroot);
									//~ if(true) {
									//~ if(true) {
										//~ std::cout  << "obtaining new Newick representation from edgytree and writing it to " << NewickFileName.toAscii().constData() << std::endl;
										//~ std::ofstream NewickFile (NewickFileName.toAscii().constData(), std::ios::out); // replace if exists!
										//~ if (NewickFile.is_open()) {
											//~ etp->NEXUS_out(NewickFile); // edgytree::NEXUS_out() writes in the simpler Newick format rather than NEXUS
											//~ NewickFile.close();
										//~ } else {
											//~ cout << "Error opening " << NewickFileName.toAscii().constData() << " for writing" << endl;
											//~ return;
										//~ }
									//~ }

                                    //~ QString command = "drawgram " +  // botanist warns user as first thing if no drawgram in path
									//~ NewickFileName + " " + // input file
									//~ XBMFileName;   // output file
									//~ system(command.toAscii().constData());

									//~ command ="";
									//~ command = "convert -flip " +  // botanist warns user as first thing if no drawgram in path (-rotate -180)
									//~ XBMFileName + " " + // input file
									//~ jpgName;   // output file
									//~ //QMessageBox::information(this,"",labelTab,"Fine");		
									//~ system(command.toAscii().constData());
									//~ }

	
		}//end		 
} 

/******************************************************************************
* This function will find a clas equivalent R'_Delta. It get a root of a
* subtre and set alive=false to root and all vertex between root and you leave
******************************************************************************/
void find_R_deltaEquivalenteRelation(historynode* v){
		if(bv::debug)	printf(" (*)-[find_R_deltaEquivalenteRelation] \n");
    if (v==hroot) return;	
		std::queue<historynode*> q;
		v->alive=false;
		q.push(v);
    historynode* np; // node pointer
    historynode* chp; // child pointer
    while(!(q.empty())){
			np = q.front();
			q.pop();	
			// if the vertex is not a leaf, is necessary look in all branch
			if (np->hfirstchild!=0){
				// part 1: firschild
				historynode* fch = np->hfirstchild;
				// Parameter (ancestor,node)	
				chp = find_R_delta(np,fch);
                //have some critical vertex. Mark all false and put a critical vertece at queue
				if (chp!=np){
					historynode* aux = chp;	
					while (aux!=np){
						aux->alive=false;
						aux = aux->hparent;
					}	
					q.push(chp);	
				}
				
				// part 2: another children	
				while (fch->hnext!=0){					
					chp = find_R_delta(np,fch->hnext);	
					//have some critical vertex. Mark all false
					if (chp!=np){
						historynode* aux = chp;	
						while (aux!=np){
							aux->alive=false;
							aux = aux->hparent;
						}							
						q.push(chp);		
					}						
					fch = fch->hnext;
				}				
			}
		}//end queue
}

/******************************************************************************
* This function find a R delta equivalence
* Parameter: ancestor and node (pai, filho)
******************************************************************************/
historynode* find_R_delta(historynode* u, historynode* v){
		
	historynode* aux = v;
	while ( aux->gray<= ( u->gray + this->omega) ){		
		if(bv::debug)	printf("(*)-[find_R_delta] [ %lld, %lld, %lld] \n",u->gray, aux->gray,this->omega);	
		if (isCriticalVertex(aux)){			
				return aux;				
		}
		aux = aux->hfirstchild;
			
	}		
	return u;
}


/******************************************************************************
* This function find a vertex that satisfy R_delta binary relation. 
* The function try find for a critical vertex in a branch indicate by anc. 
* If it dont find anything, it return itself u= ancentor v= vertex
******************************************************************************/
void flattenClass(){
	if(bv::debug)	std::cout << "flattenClass \n " << std::endl;		
	std::queue<historynode*> q;
	q.push(hroot);
	historynode* np; // node po
	while(!(q.empty())) {
		np = q.front();
		q.pop();	
		std::cout << "\n a - ";	
		if (isCriticalVertex(np)){
			std::cout << " b - ";		
			historynode* aux = np->hfirstchild;				
			while (aux!=0 ){
				if(bv::debug)	std::cout << " CV  ("<< static_cast<long long> (aux->gray ) << 	"," << static_cast<int> (aux->size )<< ")" <<  std::endl; 									
				q.push(aux);
				aux = aux->hnext;
			}
		}else{	
			std::cout << " c - ";		
			historynode* son = np->hfirstchild;		
			historynode* pai = np->hparent;
			if (np!=hroot){	
				std::cout << " d - ";		
				numnodes--;	
				for(int i=0;i<np->listCoord.size();i++){
					pai->listCoord.push_back(np->listCoord[i]);
				}		
				removeNode(np);	
			}
			q.push(son);	
		}
	}		
}

void makeClass(){
	if(bv::debug)	std::cout << "makeClass \n " << std::endl;		
	std::queue<historynode*> q;
	q.push(hroot);
	historynode* np; // node po
	while(!(q.empty())) {
		np = q.front();
		q.pop();	
		// is not a hroot	
		if (np!=hroot){
			// alive false <->. This node is a root of a class C
			if (!np->alive){
				compressClass(np);
				historynode* aux = np->hfirstchild;				
				while (aux!=0 ){
					if(bv::debug)	std::cout << ">>>  ("<< static_cast<long long> (aux->gray ) << 	"," << static_cast<int> (aux->size )<< ")" <<  std::endl; 									
					q.push(aux);
					aux = aux->hnext;
				}
			// you are not sinlge node 		
			}else{
				historynode* aux = np->hfirstchild;
				while (aux!=0 ){
					q.push(aux);
					aux = aux->hnext;
				}								
			  //removeNode(np);
			}	
		// is the hroot		
		}else{
			historynode* aux = np->hfirstchild;
			while (aux!=0 ){
				q.push(aux);
				aux = aux->hnext;
			}				
		}
	}
	//myTravesa();	
}

/******************************************************************************
* This function receive a critical vertex and remove all your descendant that are 
* in the same equivalence class.
******************************************************************************/

void compressClass(historynode* u){
	std::queue<historynode*> q;
	std::stack<historynode*> sRem;
		
	//if(bv::debug)	std::cout << "compressClass {"<< static_cast<int> (u->gray ) <<	"," << static_cast<int> (u->size )<< "}" <<  std::endl;	
	if(bv::debug)	printf("compressClass > { %lld, %d)\n",u->gray,u->size);	
	u->alive = true;
	if (u->hfirstchild!=0){
		q.push(u->hfirstchild);		
		historynode* net = u->hfirstchild;
		while (net->hnext!=0 ){
			q.push(net->hnext);
			net = net->hnext;
		}
	}
	historynode* np; // node po
	while(!(q.empty())) {
		np = q.front();
	  //if(bv::debug)	std::cout << "---- [a] {"<< static_cast<int> (np->gray ) << "," << static_cast<int> (np->size )<< "}" <<  std::endl;			
		if(bv::debug)	printf("compressClass >>{ %lld, %d)\n",np->gray,np->size);
		q.pop();	
		// alive false
		if(!np->alive){
			if ( np->gray<= ( np->hparent->gray + this->omega) ){	
				sRem.push(np);	
				removeNode(np);		
			}
			historynode* aux = np->hfirstchild;
			while (aux!=0 ){
	  		//if(bv::debug)	std::cout << "---- ---- [l] {"<< static_cast<int> (aux->gray ) << "," << static_cast<int> (aux->size )<< "}" <<  std::endl; 	
				if(bv::debug)	printf("---- ---- [l]  { %d, %d)\n",aux->gray,aux->size);
				if (!aux->alive){
					q.push(aux);
					if(bv::debug)	printf("---- ---- ---- [q] { %d, %d)\n",aux->gray,aux->size);	
	  			//if(bv::debug)	std::cout << "---- ---- ---- [q] {"<< static_cast<int> (aux->gray ) << " ," << static_cast<int> (aux->size )<< "}" <<  std::endl; 						
				}
				aux = aux->hnext;
			}
			if(bv::debug)	printf("--- ---- ---- ---- [r] { %d, %d)\n",np->gray,np->size);	
	  	//if(bv::debug)	std::cout << "---- ---- ---- ---- [r] {"<< static_cast<int> (np->gray ) << "," << static_cast<int> (np->size )<< "}" <<  std::endl;
						
		}				
	}
	//~ while(!(sRem.empty())) {	
		//~ historynode* aux = sRem.top();
		//~ sRem.pop();
		//~ if(bv::debug)	printf("---- ---- ---- ---- [reeeeeeeeeeeeeeeeeeeee] {%d, %d)\n",aux->gray,aux->size);	
	  //~ //if(bv::debug)	std::cout << "---- ---- ---- ---- [reeeeeeeeeeeeeeeeeeeee] {"<< static_cast<int> (aux->gray ) << "," << static_cast<int> (aux->size )<< "}" <<  std::endl;
		//~ removeNode(aux);			
	//~ }		
}


/** MATAV3
****************************************************************************
* Funtion that implement de FastPruning algotithm presented by Yung Kong
in march 2009. The description on this algorithm are in the paper 
FHT-fast_pruning.pdf.
// validAncestors[] is an array that will be used to hold pointers to vertex
// objects, for use in the procedure OutTree(). It is assumed to be
// an array of size at least K, where K is the maximum, over all leaves x of
// T, of the number of critical ancestors of x in T.
******************************************************************************/

 void fastPruningFHT(historynode** validAnc, int delta, int delta1) {
    printf("fastPruningFHT() :-) \n");
  
	if ( (delta==0) && (delta1==0) ){
		printf(" The valeu of delta is ZERO. The tree will not be pruned. \n");	
	}
	historynode* m = hroot; //the critical vertex of T that is an ancestor of every critical vertex of T
	validAncestors = validAnc;		 
	this->delta1 = delta1;
    this->omega = delta;
    // set the IPCD for all critical vertex.
    if (isCriticalVertex(m)){
        //computeIPCD(m,1);
        if(bv::debug) printf("the root is CV (%lld) \n",m->gray);
    }else{
        m=	immediateCriticalVertice(m);
        //computeIPCD(m,1);
        if(bv::debug) printf("the root is NOT a CV. The m is= (%lld) \n",m->gray);
    }

    /* Compute the farthestDescendent and DeltaLimite*/
    computeFarthestDescendantsAndDeltaLimits(m);
    if (!m->farDescendant)
         printf("m->farDescendant==0 \n");
    else
         printf(" Other \n");
    printf("Information about  m%lld   (%lld) ",m->gray, m->deltaLimit);
    printf("Information about  m->farDescendant  %lld   (%lld) ",m->farDescendant->gray, m->farDescendant->deltaLimit);
    m->farDescendant->deltaLimit = this->myPosInf;
	historynode* r = m;
	while (	(r->deltaLimit <= delta) && (r->hfirstchild!=0) ){
        if(bv::debug) printf("the root is NOT a CV. The m is= (%lld) \n",r->gray);
		r = r->successor;
	}
	printf("Informationa about m -->  gray Max %lld   (%lld) - %lld",m->farDescendant->gray, m->farDescendant->deltaLimit, r->gray);	 
	
	/* Call computeOutTree for all element of r.IPCD*/
	validAncestors = 	malloc((computeMaxK()+1)*sizeof(historynode*));  
	validAncestors[0] = r;
	computeMaxSizeNewIPCD(r); 
	r->label = this->myPosInf;
	allocatedNewIPCD(r);
	
	for (int i=0;i<r->sizeIPCD;i++){
			computeOutTree(r->IPCD[i],0,r);
	}
    outputNewFHTRootedAt(r,"outAqui");
 }


/** MAVAv3
---------------------------------
*/
void computeOutTree(historynode* thisVertex, int lastIndex, historynode* closestAnc){
    //printf(" computeOutTree (%lld,%d,%lld)\n",thisVertex->gray,lastIndex, closestAnc->gray);
	if (thisVertex->farDescendant->deltaLimit > this->omega){ 
		historynode* v = thisVertex;
		while ( (v->deltaLimit <= this->omega) && (v->hfirstchild!=0 ) ) {
			v = v->successor;
		} 
		int n = lastIndex + 1;
        do{
			n = n - 1;
			v->label = v->gray - validAncestors[n]->gray;
		}while ( !(v->label > this->delta1) && !(v->label <= validAncestors[n]->label));
		validAncestors[lastIndex + 1] = v;
        if (v->label > this->delta1){
            addNewIPCD(closestAnc,v);
			v->newIPCD=0;	
			for (int i=0; i< v->sizeIPCD; i++){
				computeOutTree(v->IPCD[i], lastIndex+1,v);
			}
		}else{
			for (int i=0; i< v->sizeIPCD;i++){
				computeOutTree(v->IPCD[i], lastIndex+1,closestAnc);
			}
		}
	}
	
}

/** MATAv3
----Informal descriptions of what the procedures computeFarthestDescendantsAndDeltaLimits()
NOTE:
 $  We write T_c to denote the subtree of T that is rooted at c.
 $  We write T(delta) to denote the tree of the FHT structure that is obtained when we prune (T,L,S) by branch
    length delta, using an ordering of the leaves of T such that, whenever a.graylevel = b.graylevel for two
    leaves a and b, leaf a precedes leaf b if preorder traversal of T would visit b before a. (Here preorder
    is defined with respect to the order of vertices in v.IPCD for each critical vertex v of T.) Note that
    T(0) = T.
  $ We write T_c(delta) to denote the subtree of T(delta) that is induced by the vertices of T(delta) which
    are vertices of T_c. Thus T_c(delta) will denote the empty tree if no vertex of T_c is a vertex of T(delta),
    but otherwise T_c(delta) will denote the subtree of T(delta) that is rooted at c. [In the latter case, the
    root c of T_c(delta) will only be a critical vertex of T_c(delta) if c is a critical vertex of T(delta).]
    Note that T_c(0) = T_c

Let c be any critical vertex of T. Then a call of computeFarthestDescendantsAndDeltaLimits(c) sets the values
of v.farthestDescendant, v.deltaLimit, and v.successor for each critical vertex v of T_c, except that
v.successor is not set if v is a leaf, and c.farthestDescendant.deltaLimit is not set by the call
computeFarthestDescendantsAndDeltaLimits(c) but will be set later.

For each critical vertex v of T_c, the values assigned by computeFarthestDescendantsAndDeltaLimits(c) to
v.farthestDescendant, v.deltaLimit, and v.successor satisfy the following conditions:
    (1) v.farthestDescendant is set to a descendant u of v in T for which u.graylevel is maximal.
        (This implies v.farthestDescendant will always be a leaf.) If there is more than one such
        descendant u of v, then v.farthestDescendant is set to whichever one of those descendants would
        be visited first in a preorder traversal of T.

    (2) If v != c.farthestDescendant, then v.deltaLimit is set to the LEAST value of delta for which v
        is not a critical vertex of T_c(delta); thus the value of v.deltaLimit will be such that v is a
        critical vertex of T_c(delta) if and only if v.deltaLimit > delta. Since T_c(0) = T_c, the value of
        v.deltaLimit will be strictly positive.

    (3) If v is not a leaf, then v.successor is set to the immediate proper critical descendant of v that
        lies on the path from v to v.farthestDescendant. If v is a leaf, then v.successor is not set.

Note that Algorithm 2 executes computeFarthestDescendantsAndDeltaLimits(m) m.farthestDescendant.deltaLimit=INFINITY
before entering its main loop, where m is the critical vertex of T that is an ancestor of every critical
vertex of T, and note that the deltaLimit fields of vertices are not changed during execution of the main loop.
Note also that the critical vertices of T_m and of T_m(delta) are just the critical vertices of T and of
T(delta), respectively. It follows from this and (2) that, for any positive delta, during execution of the main
loop the critical vertices v of T for which v.deltaLimit > delta are the critical vertices of T(delta).
**/
void computeFarthestDescendantsAndDeltaLimits(historynode* vertex){
    if(bv::debug) printf(" computeFarthestDescendantsAndDeltaLimits(%lld) \n	",vertex->gray);
    if ( vertex->sizeIPCD ==0 ){
        if(bv::debug) printf("LEAF CV %lld \n	",vertex->gray);
        vertex->farDescendant = vertex;
    }else{
        long long currentMax = this->myNegInf;
        if(bv::debug) printf("NORMAL CV gray=%lld (%d)) \n	",vertex->gray, vertex->sizeIPCD);

        for (int i=0;i< vertex->sizeIPCD; i++){
             if(bv::debug) std::cout << "  ---->Pai=" << static_cast<int> (vertex->gray ) << " filho" << static_cast<int> (vertex->IPCD[i]->gray ) << "] (" <<currentMax<< ")" <<std::endl;
            computeFarthestDescendantsAndDeltaLimits(vertex->IPCD[i]);
            if(bv::debug) printf("vertex->IPCD[i]->farDescendant %lld  [[[%lld]]](%lld) \n	",vertex->IPCD[i]->farDescendant->gray,currentMax,this->myNegInf);
            //compute the successor
            if (vertex->IPCD[i]->farDescendant->gray > currentMax){
                    currentMax = vertex->IPCD[i]->farDescendant->gray;
                    vertex->successor = vertex->IPCD[i];
                     if(bv::debug) printf(" i =%d  %lld is successor of %lld\n",i,vertex->successor->gray,vertex->gray );
            }
        }
        vertex->farDescendant = vertex->successor->farDescendant;
         if(bv::debug) if(bv::debug) printf("gray=%lld  vertex->successor->farDescendant->vertex->gray (%lld)) (%d)) \n	",vertex->gray,vertex->successor->farDescendant->gray);
        //compute deltaLimit
        vertex->deltaLimit = 0;
        for (int i=0;i<vertex->sizeIPCD; i++){
            if(bv::debug) std::cout << "      +++++DELTA LIMITE  >Pai=" << static_cast<int> (vertex->gray ) << " filho" << static_cast<int> (vertex->IPCD[i]->gray ) << "]" <<std::endl;
            if (vertex->IPCD[i] != vertex->successor){
                vertex->IPCD[i]->farDescendant->deltaLimit = vertex->IPCD[i]->farDescendant->gray - vertex->gray;
                vertex->deltaLimit = max(vertex->deltaLimit, vertex->IPCD[i]->farDescendant->deltaLimit);
                if(bv::debug) printf(" ++++++++++ i =%d  delta limite of %lld  is %lld\n",i,vertex->gray, vertex->deltaLimit );
                 if(bv::debug) printf(" ++++++++++ i =%d  delta limite of clild %lld  is %lld\n",i,vertex->IPCD[i]->farDescendant->gray , vertex->IPCD[i]->farDescendant->deltaLimit);
            }
        }
    }
}


void outputNewFHTRootedAt(historynode* r,QString name){
    printf("(*)-> outputNewFHTRootedAt (%lld) \n",r->gray);

    edgytree* etp = new edgytree(); 
    historynode* hnp;
    historynode* hchp;
    edgynode* enp;
    edgynode* echp;
    std::stack<historynode*> hns; // history node stack
    std::stack<edgynode*> ens;  // edgy node stack

    hns.push(r);

    etp->eroot = new edgynode();
    edgynode* enpF = new edgynode();
    enpF->eparent = etp->eroot;
    etp->eroot->efirstchild = enpF;
    etp->eroot->elastchild = enpF;
    enpF->edgetoparent = ABSDIF(r->gray, this->hroot->gray); // root (and only root) has no parent.
    etp->eroot->edgetoparent = 0;
    ens.push(enpF);
    // deep copy ONLY 'alive' nodes
    while(!hns.empty()) {
			// sanity check:
			if(ens.empty()) { 
				//std::cerr << "Non-empty hns but empty ens.  Returning null." << std::endl;
				return 0;
			}
			hnp = hns.top();
			hns.pop();
			enp = ens.top();
			ens.pop();

			if(hnp->sizeNewIPCD!=0) {
				echp = new edgynode();
				enp->elastchild = echp;
				enp->efirstchild = 	echp ;
				echp->eparent = enp;
					
				echp->edgetoparent =  hnp->newIPCD[0]->gray - hnp->gray ;	

				ens.push(echp);	
				hns.push(hnp->newIPCD[0]);			
				for (int i=1; i< hnp->sizeNewIPCD;i++){
					edgynode* aux = new edgynode();					
					echp->enext = aux; 	
					aux->eparent = enp;	
					enp->elastchild = aux;	
					aux->edgetoparent =  hnp->newIPCD[i]->gray - hnp->gray;
					ens.push(aux);	
					hns.push(hnp->newIPCD[i]);						
					echp=aux;	
				}											
			}else{
			  //printf("------------------------------------- \n");
			}
   }


    etp->flatten(); // once again, does not remove 'edge' from root down
    etp->setparents();

    QString NewickFileName = name+".new";
    QString XBMFileName =  name+".xbm";
    QString jpgName = name+".jpg";
    if(true) {
        std::ofstream NewickFile (NewickFileName.toAscii().constData(), std::ios::out); // replace if exists!
        if (NewickFile.is_open()) {
            etp->NEXUS_out(NewickFile); // edgytree::NEXUS_out() writes in the simpler Newick format rather than NEXUS
            NewickFile.close();
        } else {
            cout << "Error opening " << NewickFileName.toAscii().constData() << " for writing" << endl;
            return;
        }
    }

    QString command = "drawgram " +  // botanist warns user as first thing if no drawgram in path
    NewickFileName + " " + // input file
    XBMFileName +" 800 600 ";   // output file
    system(command.toAscii().constData());
    command ="";
    command = "convert -flip " +  // botanist warns user as first thing if no drawgram in path (-rotate -180)
    XBMFileName + " " + // input file
    jpgName;   // output file
    system(command.toAscii().constData());

    //printf("---------------------------------------");
   std::stack<edgynode*> s;

}	

/** MATAV3
This function compute the list of IPCD (immediate proper critical descendants)
of a node c in T. This field is never changed by the algorithm.

*/
void computeIPCD(historynode* np, int k){
	if (np->hfirstchild==0){
		np->sizeIPCD=0;
		np->maximumK = k;	
	}else{
		historynode* fcp = np->hfirstchild;
		historynode* icv = immediateCriticalVertice(fcp);	
		computeIPCD(icv,k+1);	
        alocateIPCD(np);
		addIPCD(np,icv);	
		if (fcp->hnext!=0){
			historynode* aux = fcp->hnext;	
			while(aux!=0){
				icv = immediateCriticalVertice(aux);	
				computeIPCD(icv,k+1);	
				addIPCD(np,icv);					
				aux=aux->hnext;
			}
		}
	}
}

/** MATAv3
 This function set IPCD for the vertice np.
*/
void computeMaxSizeNewIPCD(historynode* np){
	if (np->hfirstchild==0){
		np->maxSizeNewIPCD = 1;	
		//std::cout << ">  [externo]<" << static_cast<int> (np->gray )	<< ">maxSizeNewIPCD  <" << static_cast<int> (np->maxSizeNewIPCD )	<<std::endl;		
	}else{
		historynode* fcp = np->hfirstchild;
		historynode* icv = immediateCriticalVertice(fcp);	
			
		computeMaxSizeNewIPCD(icv);	
	  int size_inc = icv->maxSizeNewIPCD;
		if (fcp->hnext!=0){
			historynode* aux = fcp->hnext;	
			while(aux!=0){
				icv = immediateCriticalVertice(aux);	
				computeMaxSizeNewIPCD(icv);	
				size_inc = size_inc + icv->maxSizeNewIPCD;					
				aux=aux->hnext;
			}
		}
		np->maxSizeNewIPCD = size_inc;	
		//std::cout << ">  [interno]<" << static_cast<int> (np->gray )	<< ">maxSizeNewIPCD  <" << static_cast<int> (np->maxSizeNewIPCD )	<<std::endl;	
	}
}

int computeMaxK() {
    std::queue<historynode*> s;
	  s.push(hroot);	 
		int cont =0; 
    while(!(s.empty())) {
			historynode* np_; // node pointer
			historynode* chp_; // child pointer
			np_ = s.front();	
			s.pop();

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
				  if (cont < np_->maximumK) 
						cont = np_->maximumK;
				}
		}	
	///std::cout << "(*)-[ maximumK = " << cont <<  std::endl;	
	return cont;	
} 

/** MATAv3
This function add a critical node for the list of IPCD (immediate proper critical descendants)
*/
void addIPCD(historynode* vertex, historynode* ipcd){
	//std::cout << "-- <" << static_cast<int> (vertex->gray ) << ","<< static_cast<int> (ipcd	->gray ) << "," << static_cast<int> (vertex->sizeIPCD )<<"]" <<std::endl;	
	if (vertex->maxSizeIPCD < vertex->sizeIPCD){
		printf("Erro in addIPCD \n");
		exit(0);
	}else{
		vertex->IPCD[vertex->sizeIPCD]=	ipcd;
		vertex->sizeIPCD++;	
	}
}

/** MATAv3
Add ipcd in the vertex.newIPCD
*/
void addNewIPCD(historynode* vertex, historynode* ipcd){
    // std::cout << " ****addNewIPCD()<" << static_cast<int> (vertex->gray ) << ","<< static_cast<int> (ipcd	->gray ) << "," << static_cast<int> (vertex->sizeNewIPCD )<<"]" <<std::endl;
	if (vertex->maxSizeNewIPCD < vertex->sizeNewIPCD){
        printf("Erro in addNewIPCD (maxSizeNewIPCD= %d <	 sizeNewIPCD=%d) \n",vertex->maxSizeNewIPCD,vertex->sizeNewIPCD);
        exit(0);
	}else{		
		if (vertex->sizeNewIPCD==0){
            vertex->newIPCD = malloc((vertex->maxSizeNewIPCD+1)*sizeof(historynode*));
			vertex->newIPCD[0]=	ipcd;	
			vertex->sizeNewIPCD++;		
		}else{
			vertex->newIPCD[vertex->sizeNewIPCD]=	ipcd;
			vertex->sizeNewIPCD++;	
		}
	}
}
/** MATAv3
This function alocate memory to for the list of IPCD (immediate proper critical descendants) for a node c
*/
void alocateIPCD(historynode* np){
    int size =1;
    historynode* aux = np->hfirstchild;
    while(aux->hnext!=0){
        size++;
        aux=aux->hnext;
    }
    np->maxSizeIPCD=size;
    np->IPCD = 	malloc((size)*sizeof(historynode*));
    if (np->IPCD==NULL ){
         printf("Erro in addNewIPCD: Malloc excption!!\n");
         exit(0);
    }
}

/** MAVAv3
 clear memory for nerIPCD
 */
void allocatedNewIPCD(historynode* vertex){
	vertex->newIPCD = 	malloc(vertex->maxSizeNewIPCD*sizeof(historynode*));	
	vertex->sizeNewIPCD=0;	
}

int computeLeaveSubTree(historynode* nod) {
    std::queue<historynode*> s;
	  s.push(nod);	 
		int cont =0; 
    while(!(s.empty())) {
			historynode* np_; // node pointer
			historynode* chp_; // child pointer
			np_ = s.front();	
			s.pop();

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
				cont++;
			}
		}	
	///std::cout << "(*)-[ maximumK = " << cont <<  std::endl;	
	return cont;	
} 

/******************************************************************************
* Funtion that traverse a tree and show the gray and size of each node
******************************************************************************/
 void testando(historynode* up) {
 /*    std::stack<historynode*> s;
	 printf("(TETE] <no> [pnc:[%lld,%d] \n",up->gray,up->size);	  
	 s.push(up);	 
   
	int contador=0; 
    while(!(s.empty())) {
		historynode* np_; // node pointer
		historynode* chp_; // child pointer
		np_ = s.top();	
		std::cout << "(*) <no> [" << static_cast<int> (np_->gray ) << "," << static_cast<int> (np_->sizeIPCD )<<"]" <<std::endl;	
		s.pop();
		if (np_->hfirstchild!=0){	
			chp_ = (np_->hfirstchild);
			while(chp_ != 0) { // and add ptrs to all children to the stack
				if (isCriticalVertex(chp_)) s.push(chp_);		
					chp_ = (chp_->hnext);
				}
			}
	  }*/
} // --historytree::testando()

 edgytree* generate_pruned_edgy_treeFHCS() {
    if(bv::debug) {
       std::cout << "generate_pruned_edgy_tree() h is  "
         << static_cast<unsigned int> (prunedbyheight)
         << " in generate_pruned_edgy_tree()" << std::endl;
     }
     edgytree* etp = new edgytree();
     if(!(hroot->alive)) return etp;
     historynode* hnp;
     historynode* hchp;
     edgynode* enp;
     edgynode* echp;
     std::stack<historynode*> hns; // history node stack
     std::stack<edgynode*> ens;  // edgy node stack
     historynode* r = hroot;
     if (isCriticalVertex(hroot)){
         //hns.push(r);
         historynode* m = r;
         while (	(r->deltaLimit <= delta1) && (r->hfirstchild!=0) ){
         //while (	( (r->gray - m->gray) <= delta1) && (r->hfirstchild!=0) ){
             r = r->successor;
         }
         hns.push(r);
     }else{
         r = immediateCriticalVertice(hroot);
         while (	(r->deltaLimit <= delta1) && (r->hfirstchild!=0) ){
             r = r->successor;
         }
         hns.push(r);
     }
     etp->eroot = new edgynode;
     etp->eroot->edgetoparent = 0; // root (and only root) has no parent.

     etp->eroot->elastchild = etp->eroot->efirstchild = echp = new edgynode;
     printf(" root=%lld  r=%lld",hroot->gray,r->gray);
     echp->edgetoparent = ABSDIF(r->gray, hroot->gray);
     ens.push(echp);
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
        //printf("hnp=%lld \n",hnp->gray);
        if(hnp->sizeNewIPCD!=0) {
         edgynode* epp = 0; // ptr to edgy "predecessor" (left-sibling) node
         for (int i=0; i< hnp->sizeNewIPCD;i++){ // and add ptrs to selected children to stacks
                 hchp = (hnp->newIPCD[i]);
                 hns.push(hchp);
                 enp->elastchild = echp = new edgynode;
                 echp->edgetoparent = ABSDIF(hnp->newIPCD[i]->gray, hnp->gray);
                 ens.push(echp);
                 if(epp==0) {
                     enp->efirstchild = echp;
                 } else {
                     epp->enext = echp;
                 }
                 epp = echp;
                 //echp = (echp->enext);
         }
         //enp->elastchild = epp;
        }else{
        enp->elastchild = enp->efirstchild = 0;
        }
     }
    printf("generate_pruned_edgy_treeFHCS() END \n");
    return etp;
 } // --generate_pruned_edgy_treeFHCS()



 /*******************************************************************************
 * FUNCTION FOR THE NEW TREE STEP SIMPLIFICATION  VERSION 2. THIS
 This following function are used to prune the tree by branch and by collpsing
 equivalent class. The parameter used is delta (long). The tree could be buit
 using allgray values or a number os level between 1-255.
 *******************************************************************************/


 /******************************************************************************
 * This function prune the tree that is the output of step 1 (by component size)
 * by  simplification branch length delta (implements a BST to choice
 * the enumeration of leave visited).
 ******************************************************************************/
  void pruneByDeltaStepB(int omega) {
     std::queue<historynode*> q;
         this->omega = (omega);
         //std::cout << ":: Simplification by Branch Length Delta {" << this->omega << std::endl;


     if(!q.empty()) {      // should NEVER happen
       //std::cerr << "huge error: stack s non-empty at beginning"
                 //<< " of historytree::prune()" << std::endl;
       throw std::runtime_error("stack nonempty at beginning of historytree::prune()");
     }

          if(hroot == 0) {
       // std::cerr << "in historytree::prune() hroot is null!!!" << std::endl;
       throw std::logic_error("in historytree::prune() hroot is null");
     }

     q.push(hroot);
     historynode* np; // node pointer
     //historynode* chp; // child pointer

     while(!(q.empty())) {
             np = q.front();
             q.pop();
             //checking if the node is leaf
             if (np->hfirstchild == 0 && np->hlastchild == 0){
                     historynode* ppa= findProperAncestor(np);
                     //had a properAncestor
                     if (ppa != np){
                         int bef = ppa->listCoord.size();
                         int rem = moveVoxels_in_pruneByDeltaStepB(np,ppa);
                          printf("rem=%d hnp=%lld \n",rem,np->gray);
                         removedPath(np,ppa);
                     }

             // is not a leave.Then all children will be pushed in the queue
             }else{
                 historynode* chp = np->hfirstchild;
                 q.push(chp);
                 while(chp->hnext){
                     historynode* axu = chp->hnext;
                     q.push(axu);
                     chp = axu;
                 }
            }
         }
 } // --historytree::pruneByDeltaStepB()

  /******************************************************************************
  * This function move the voxe of node removed in the pruneByDeltaStepB to
  * the parent of x of the nodes between np and ppa
  *  for the nearest elements with the same level. Parameter:
  *  x	the frist children of ppa
  *  np leaf in the stream x--np
  *  ppa  -is the critical vertex of np
  ******************************************************************************/
  int moveVoxels_in_pruneByDeltaStepB(historynode* np, historynode* ppa){

    int paiA = ppa->listCoord.size();
      int paiD = 0;

      //loop for move components from leaf critical vertic np until x
      if (np == ppa) printf("Possivel erro in moveVoxels_in_pruneByDeltaStepB \n");
      int voxelsMoved =0;
      while (np != ppa) {
          for(int i=0;i<np->listCoord.size();i++){
              ppa->listCoord.push_back(np->listCoord[i]);
              voxelsMoved++;
              paiD++;
          }
          np->listCoord.clear();
          np = np->hparent;
      }
      return voxelsMoved;
  }

  /******************************************************************************
  * This funciont remove each edge between the properAncestor and a vertex
  * if the properAncestor have more than one children.  (np==node) ppa==properAncestor
  ******************************************************************************/
  void removedPath(historynode* np, historynode* ppa){

      // find the first children of 'x'.
      numnodes--;
      historynode* x = np->hparent;

      if (np->hparent == ppa){
       x = np;
      }else{
          while (x->hparent != ppa &&  x->hparent->hfirstchild==x->hparent->hlastchild){
              historynode* aux = x->hparent;
              x = aux;
              numnodes--;
          }
      }

      removeFromTree(x);
  }


  /******************************************************************************
  * This function move the componente size of the nodes between np and ppa
  *  for the nearest elements with the same level. Parameter:
  *  np   -is the removed nodes
  *  ppa  -is the critical vertex of np
  *  x is -the children of ppa (ancestor of X)
  ******************************************************************************/
  void moveComponentSize(historynode* x, historynode* np, historynode* ppa){

      //loop for move components from leaf critical vertic np until x
      while (np != ppa) {
          unsigned size_removed = np->size;
      //>> x is lastchild. Need looking for the preview
          if (x->hnext==0){
              historynode* tio = ppa->hfirstchild;
              while (tio->hnext != x ){
                  tio = tio->hnext;
              }
              while (tio->height != np->height){
                  tio = tio->hlastchild;
              }

              tio->size = tio->size	+ size_removed;
          //>> x is not the lastchild. Then get next
          }else{
              historynode* aux = x->hnext;
              //find a node with same level
              while (aux->height != np->height){
                  aux = aux->hfirstchild;
              }
              aux->size = aux->size	+ size_removed;
          }
          np = np->hparent;
      }
  }

  /******************************************************************************
  * The function update a gray value of when a node is removed.  NOT_CALLED
  *
  ******************************************************************************/
  void updateGray(historynode* up) {
    std::stack<historynode*> s;

    if(up->hfirstchild == 0) {
      up->gray--;
          return;
      }
      s.push(up);
      setalive = 0;
      while(!(s.empty())) {
          historynode* np_; // node pointer
          historynode* chp_; // child pointer
          np_ = s.top();

          s.pop();
          np_->gray--;
          chp_ = (np_->hfirstchild);
          while(chp_ != 0) { // and add ptrs to all children to the stack
              s.push(chp_);
              chp_ = (chp_->hnext);
          }
      }
  }

}; // --struct historytree

#endif
