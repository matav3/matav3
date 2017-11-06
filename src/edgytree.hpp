#ifndef EDGYNODEANDTREE
#define EDGYNODEANDTREE

#include <cassert>
#include <queue>     
#include <map>
#include <stack>
#include <list>
#include <iostream>

using namespace std;
/** @date last edited (modified comments) 2005.06.22 
 *        and again 2006.03.30
 *  @description struct edgynode is the node of an edgytree 
 */
struct edgynode {
  edgynode* enext;
  edgynode* eparent;
  edgynode* efirstchild;
  edgynode* elastchild;
  edgynode* einsertbefore;
  long long edgetoparent; // (length)
	int label;
	int isomorphic;	
	std::vector<int> code; 	
  edgynode() { /** no-arg ctor for edgynode */
    enext = eparent = efirstchild = elastchild = 0; // NULL
    edgetoparent = 0;
		label =isomorphic=0;
  }

  /** deep destructor for edgynode */
  ~edgynode() {
    edgynode* ch = efirstchild;
    while(ch != 0) { // while still have children
      edgynode* nextch = ch->enext; // remember next child
      delete ch; // get rid of subtree rooted at current child by means of recursion
      ch = nextch; // make next child current child
    }
  }
}; // -- struct edgynode


/** @description struct edgytree is a tree with edge lengths */
struct edgytree {
  edgynode* eroot;
  unsigned fakeLabelNo;
  unsigned int minheight;
  long long omega;	
  unsigned int properAncestor; 	
  edgytree() {
    eroot = 0; // initialization saved my life... again!
    fakeLabelNo = 1; // whatever
    minheight = 0;
  }

  /** 
   *  traverse and set parents properly in edgy subtree with root pointed to by arg 
   *  THIS WAS NOT DONE ANYWHERE! 
   *  (and our 2004 pruning method NEVER used the eparent field.  amazing!)
   */
  void setparents(edgynode* enp) {
    if(enp==0) {
      return;
    }
    for(edgynode *echp = enp->efirstchild; echp != 0; echp = echp->enext) {
      echp->eparent = enp;
      setparents(echp);
    }
  } // --edgytree::setparents()/1

	
  void setparents() {
    eroot->eparent = 0;
    setparents(eroot);
  } // --edgytree::setparents()/0

  /** traverse edgy subtree with root pointed to by arg */
  void traverse_subtree(edgynode* enp) {
    // maybe move this fnc (and historytree counterpart) into node classes?
    if(enp==0) {
      std::cerr << "null ptr passed to edgytree::traverse_subtree()" << std::endl;
      throw std::invalid_argument("null ptr passed to edgytree::traverse_subtree()");
    }
    std::cout << "Length of edge to parent: " 
	      << enp->edgetoparent  << std::endl;
    edgynode *echp = enp->efirstchild;
    if(echp) {
      std::cout << "E-Children: { " << std::endl;
      while(echp != 0) {
	traverse_subtree(echp);
	echp = echp->enext;
      }
      std::cout << "}" << std::endl;
    }
  } // --edgytree::traverse_subtree()

  void traverse() {
    if(eroot == 0) {
      std::cerr << "eroot is null in edgytree::traverse()" << std::endl;
      return;
    }
    std::cout << "Traversing a edgytree." << std::endl;
    traverse_subtree(eroot);    
  } // --edgytree::traverse()

  /** pre: enp points to a valid edgynode                              *
   * post: will CHANGE edgy subtree irreversibly (make it pruned)      *
   * returns height of subtree BEFORE the pruning, as needed by callee */
  unsigned int pruneSubtreeByHeight(edgynode* enp) {
    //    std::cout << "++edgytree::pruneSubtreeByHeight()" << std::endl;
    if(0==enp) std::cerr << "error in edgytree::pruneSubtreeByHeight: enp is NULL!" << std::endl;
    unsigned int heightOfThisNode = 0;
    // for each child
    // calculate distancetochild plus child's height (recursively)
    edgynode* chp = enp->efirstchild;
		if (chp==0)
		  	std::cout << "(*) [A]  <no> " << enp->edgetoparent << std::endl;
		else
		    std::cout << "(*) [B]  <pai><filho> " << enp->edgetoparent << ", " <<  chp->edgetoparent << std::endl;
	  
    edgynode* lastsurvivingchp = 0; // last KNOWN survivor, will be 0 if none.
    while(chp != 0) {
		// set pointers between children after taking some out!
		// combine them iteratively
		unsigned int subtreeheight = pruneSubtreeByHeight(chp); // recurse!
		
		edgynode* nextp = chp->enext; // remember next child
		if (nextp!=0)
			std::cout << "(*) [C] <next>  " <<  nextp->edgetoparent  << std::endl;
	
			
		if(subtreeheight <= minheight) {
			delete chp; // recursively, and don't do anything with lastsurvivingch yet
		} else {
			// assert: chp not null
			if(0==chp) {
				std::cerr << "Error in edgytree::pruneSubtreeByHeight(): chp cannot be NULL here!"
				<< std::endl;
			}
			if(0 == lastsurvivingchp) {
				enp->efirstchild = chp;
			} else {
				lastsurvivingchp->enext = chp; 
			}
			lastsurvivingchp = chp;
		} // -- else
		if(subtreeheight > heightOfThisNode) heightOfThisNode = subtreeheight;
		
		chp = nextp;
		    
	} // --while, gone through all children
	// std::cout << "heightOfThisNode=" 

	//   << (static_cast<unsigned>(heightOfThisNode)) << std::endl;
	if(0 == lastsurvivingchp) {
		enp->efirstchild = enp->elastchild = 0;
	} else {
		lastsurvivingchp->enext = 0;
		enp->elastchild = lastsurvivingchp;
	}
	long long retval = heightOfThisNode + (enp->edgetoparent);
	// if minheight>=retval, then WILL be deleted anyway.
	if((heightOfThisNode < minheight) && (minheight < retval)) {
		// pretend it's the next-highup logical if subtree due to this won't be deleted by parent, but not high enough.
		enp->edgetoparent = retval - minheight;
	}
    //    std::cout << "--edgytree::pruneSubtreeByHeight()" << std::endl;
    return retval; // returns the length of the subtree of parent that goes through "this"
} // --edgytree::pruneByHeight


  /**  pruning by height of edgytree */
  void pruneByHeight(unsigned int minheight) {
    if(eroot == 0) {
      std::cerr << "eroot is null in edgytree::pruneByHeight()" << std::endl;
      return;
    }
    this->minheight = minheight;
    //    std::cout << "in edgytree::pruneByHeight(), this->minheight is " 
    //	      << (static_cast<unsigned>(this->minheight)) << std::endl;
    //    setparents();
    flatten();
    long long rootsNewHeight = pruneSubtreeByHeight(eroot);
    std::cout << "edgytree::pruneByHeight() completes execution with root's height="
	      << rootsNewHeight << std::endl;
    flatten();
    //    setparents();
  } // --edgytree::pruneByHeight()

  /**
   *  set c to be a new child of p, right before p->einsertbefore 
   *  pre: p->einsertbefore cannot be null (it will be removed from tree eventually, some time after this step.)
   */
  void setNewChild(edgynode* c, edgynode* p, unsigned int cpdist) {
    if(bv::debug)
      std::cerr << "setNewChild()... c: " << c << ", p: " << p << ", cpdist: " << (unsigned)cpdist << std::endl;
    if((p->einsertbefore) == 0) {
      std::cerr << "setNewChild() SERIOUS ERROR - einsertbefore is NULL!  impossible!";
      return;
    }
    if((p->efirstchild) == 0) {
      std::cerr << "setNewChild() SERIOUS ERROR - efirstchild is NULL!  impossible!";
      return;
    }
    (c->enext) = (p->einsertbefore);    
    (c->eparent) = p;
    (c->edgetoparent) = cpdist;
    if((p->einsertbefore) == (p->efirstchild)) {
      c->enext = p->efirstchild;
      p->efirstchild = c;
    } else { // there is at least one node before p->einsertbefore, find it and call it ch
      edgynode* ch;
      for(ch = p->efirstchild; ; ch = ch->enext) {
	if((ch->enext) == (p->einsertbefore)) break;
	if((ch->enext) == 0) {
	  std::cerr << "SERIOUS ERROR - 'einsertbefore' is not a child of the parent!";
	  return;
	}
      }
      ch->enext = c;  // insert c after ch.
    }
  }

  /** 
   * emulates Yung Kong June 2006 formulation of Gabor Herman's March 2006 algorithm
   * assuming that this edgytree represents underlying tree with edge lengths 1 each
   * to handle: make sure that removeFromTree and setNewChild do the right thing
   *            need to set einsertbefore appropriately
   */
  void simplify(long long epsilon) {
    std::cout << ":: Prune STEP 2: Yung Kong June 2006 formulation of Gabor Herman's March 2006 algorithm :: " << std::endl;
    flatten();
    NEXUS_out_subtree(eroot, std::cerr); // for debug
    std::cerr << std::endl << std::endl;

    simplify(eroot, eroot, epsilon+1, 0);

    NEXUS_out_subtree(eroot, std::cerr); // for debug, and verifies "correctness" of tree
    std::cerr << std::endl << std::endl;
  }

  /** implements the above */
  void simplify(edgynode* v, edgynode* r, long long f, long long rv) {
    // foreach child c of v
    for(edgynode* c = (v->efirstchild); c != 0; ) {
      if(v==r) {
				v->einsertbefore = c;
      }
      edgynode* nextp = c->enext; // save it here!  (c may get deleted.) AQUI
      long long rc = rv + (c->edgetoparent);
      if(rc < f) {
				simplify(c, r, f, rc);
				removeFromTree(c);
      } else {
				long long offset = rc % f;
				if(offset == 0) {
	  			if((c->eparent) != r) {
	    			setNewChild(c, r, rc);
	  			}
				} else {
					// move c 'up' by offset without having to create a new node and deleting c!
					if((c->eparent) != r) {
						setNewChild(c, r, rc - offset);
					} else {
						(c->edgetoparent) -= offset;
					}
					for(edgynode* d = c->efirstchild; d != 0; d = d->enext) {
						(d->edgetoparent) += offset;
					}
				}
				simplify(c, c, f, 0);
      	}
      c = nextp;
    }
  }



  /** 
   *  pre: v is legitimate, and if parent non-null, then parent has it in its list
   *  post: completely take v out
   */
  void removeFromTree(edgynode *v) {
    //   v->edgetoparent = ((unsigned int)v % 253); // for debug
    if(bv::debug) std::cerr << "removeFromTree()... v: " << v << ", v->edgetoparent set to " << v->edgetoparent << std::endl;
    if (v==eroot) {
      std::cerr << "in removeFromTree(): ATTEMPT TO REMOVE ROOT FROM TREE!" << std::endl;
      //   eroot = 0;
      return;
    }
    edgynode* p = (v->eparent);
    if(p != 0) {
      edgynode* c = (p->efirstchild);
      if(c==0) {
				// std::cerr << "parent has no children!" << std::endl;
				return;
      }
      if(v == c) { // v is the first child
				if(v == (p->elastchild)) { // v is the only child
	  			//std::cerr << "removing some only child from tree" << std::endl;
	  			p->efirstchild = 0;
	  			p->elastchild = 0;
				} else {
	  			p->efirstchild = v->enext;
				}
      } else { // there is at least one child that is not v
				while((c->enext) != v) {
	  			c = (c->enext);
	  			if(c==0) {
	    			// hmmm....
	    			//std::cerr << "fatal error: parent does not know about v!" << std::endl;
	    			break;
	    			//return;
	  			}
				}
				if(c != 0) {
				  if(bv::debug) std::cerr << "removeFromTree()... c: " << c << ", v: " << v << std::endl;
			 		if((p->elastchild) == v) {
						(p->elastchild) = c;
						(c->enext) = 0;
			  	} else {
						(c->enext) = (v->enext);
			  	}
				}
      }
    } else {
      //      std::cerr << "parent is null in removeFromTree() ??" << std::endl;
    }
    v->eparent = 0;
    v->enext = 0;
    v->efirstchild = 0; // these are NECESSARY
    v->elastchild = 0;
    delete(v); // calls deep deep destructor
  } // --removeFromTree()

  /** 
   *  to be called from pruneByHeightGTH0603
   *  pre: r is an ancestor of v  and this->minheight is "epsilon"
   *  post: all children of v are merged with r (merge and pruning them recursively)
   *        v is killed without damage to surrounding tissue
   */
  void mergeAndPruneGTH0603(edgynode* v, edgynode* r, unsigned int vrdist) {
    if(bv::debug)
      std::cerr << "mergeAndPruneGTH0603()... v: " << v << ", r: " << r 
		<< ", vrdist: " << (unsigned)vrdist << std::endl;
    edgynode* nextp = 0;
    for(edgynode* c = (v->efirstchild); c != 0; ) {
      long long crdist = vrdist + (c->edgetoparent);
      if(crdist == 0) { // should never get here
	std::cerr << "mergeAndPruneGTH0603() crdist: " << (unsigned)crdist
		  << ", vrdist: " << (unsigned)vrdist << ", (c->edgetoparent): " 
		  <<  (unsigned)(c->edgetoparent) << std::endl; 
	return;
      }
      nextp = (c->enext);
      if(bv::debug)
	std::cerr << "mergeAndPruneGTH0603()... c = " << c << ", nextp = " << nextp << std::endl;
      if(crdist > minheight) {
	setNewChild(c, r, crdist);
	pruneSubtreeByHeightGTH0603(c);
      } else {
	mergeAndPruneGTH0603(c, r, crdist);
	// will destroy c!... good thing we saved next ptr
      }
      c = nextp;
    }
    removeFromTree(v);
  } // --mergeAndPruneGTH0603()


  /** pre: enp points to a valid edgynode.  this->minheight is "epsilon"  
   * post: will CHANGE edgy subtree irreversibly (make it pruned)         
   * returns height of subtree BEFORE the pruning, as needed by callee 
   */
  void pruneSubtreeByHeightGTH0603(edgynode* v) {
    if(bv::debug)
      std::cerr << "pruneSubtreeByHeightGTH0603()... v: " << v << std::endl;
    if(v==0) return;
    edgynode* nextp = 0;
    for(edgynode* c = (v->efirstchild); c != 0; ) {
      long long cvdist = (c->edgetoparent);
      if(bv::debug)
	std::cerr << "pruneSubtreeByHeightGTH0603() c: " << c << ", cvdist: " << cvdist
		  << ", (c->edgetoparent): " <<  c->edgetoparent << std::endl; 
      if(cvdist == 0) {
	return;
      }
      nextp = c->enext; // remember next child (this won't be necessary anymore)
      if(cvdist > minheight) { // keep the edge from c to v as is, recurse on c
	pruneSubtreeByHeightGTH0603(c);
      } else { // merge edge from c to v with ones below, taking c out of tree
	v->einsertbefore = c; // insert everything before this.
	mergeAndPruneGTH0603(c, v, cvdist);
	// mergeAndPrune deletes c
	//   so we save "the next guy" above!
      }
      c = nextp;
    }
  } // --edgytree::pruneSubtreeByHeightGTH0603

  void pruneByHeightGTH0603(long long minheight) {
    if(eroot == 0) {
      std::cerr << "eroot is null in edgytree::pruneByHeightGTH0603()" << std::endl;
      return;
    }
    this->minheight = minheight;
    //    setparents();
    flatten();
    NEXUS_out_subtree(eroot, std::cerr); // for debug
    std::cerr << std::endl << std::endl;
    //    std::cout << "in edgytree::pruneByHeight(), this->minheight is " 
    //	      << (static_cast<unsigned>(this->minheight)) << std::endl;
    pruneSubtreeByHeightGTH0603(eroot);
    flatten(); // newly added /here/
    //    setparents();
    NEXUS_out_subtree(eroot, std::cerr); // for debug, and verifies "correctness" of tree
    std::cerr << std::endl << std::endl;
    // will need to flatten

    // don't really need to know root's height
    //    unsigned rootsNewHeight = static_cast<unsigned>(pruneSubtreeByHeight(eroot));
    //    std::cout << "edgytree::pruneByHeight() completes execution with root's height="
    //	      << rootsNewHeight << std::endl;
  }

  // now, verify too.
  void NEXUS_out_subtree(edgynode* enp, std::ostream & myout = std::cout) {
		
		edgynode* chp = enp->efirstchild;
    if(0!=chp) {
      if((chp->eparent) != enp) {
	    std::cerr << "MALFORMED tree (children don't know parents) rooted at " << eroot << std::endl;	
      }
      while(true){
				if(chp->efirstchild) {
	  			myout << "(";
	  			NEXUS_out_subtree(chp, myout);
	  			myout << ")";
				} else { // child is a leaf
	  			myout << fakeLabelNo++;	  
				}
				myout << ":" << chp->edgetoparent;
				if(chp->enext) {
	  			myout << ",";
	  			chp = chp->enext;
				} else { // done
	  			if(chp != (enp->elastchild)) {
	    			std::cerr << "MALFORMED tree (last child does not have null exnext) rooted at " << eroot << std::endl;
	    			exit(1);
	  			} else {
	    			//std::cerr << "[OK] ";
	    			break;
	  			}
				}
      } // --while
    } // --if
  } // --edgynode::NEXUS_out_subtree()


  /** output edgytree in nexus format.  please call flatten() before this... */
  void NEXUS_out(std::ostream & myout = std::cout) {
    if(eroot) {
      fakeLabelNo = 1; // initialize this instance variable
      //      std::cout << "(";
      // std::cout << "(";
      myout << "(";
      NEXUS_out_subtree(eroot, myout);
      // std::cout << ");";
      myout << ");";
    } else {
      std::cerr << "No root specified in NEXUS_out." << std::endl;
    }
    //    std::cout << std::endl;
    myout << std::endl;
  } // --edgynode::NEXUS_out()


/******************************************************************************
* Function used to read a new file and build a FHT
******************************************************************************/
  char NEXUS_in_subtree(edgytree* enp,  edgynode* pai, std::istream & myin, QString nameFile,int *lab) {
   		printf(" calling NEXUS_in_subtree\n");		
			char ch;
			long long lgt;
			myin >> ch; //is: (  => read: ( or N1
			
			/******* create the first clidren *******/
			/* After this the myin point to , or )  */
			edgynode* n1 = new edgynode();	
			//[1]---create internal node
			if (ch=='('){
				printf("first children is internal node\n");
				n1->label = (*lab)++;	
				std::cout << ch<< std::endl;	
				pai->efirstchild = n1;
				n1->eparent = pai;		
				ch = NEXUS_in_subtree(enp,n1,myin,nameFile,lab);	
				printf("closing a  interal node (saindo NEXUS_in_subtree FC) \n");
				std::cout << "tenho: " << ch<< std::endl;						
			//[1]---create a leave
			}else{ 
				printf("first  node is N1\n");
				n1->label = (*lab)++;		
				pai->efirstchild = n1;
				n1->eparent = pai;	
				std::cout << ".. checking if N1  >9- ["<<  std::endl;		
				std::cout << ch;						
				myin >> ch;	
				std::cout << ch;	
				while(isdigit(ch)){
					myin >> ch; // edge
					std::cout << ch;		
				}	
				std::cout << std::endl;	
				if (ch==':'){
					myin >> ch;
					std::cout << ".. checking if LEVEL >9: "<<  std::endl;		
					std::cout << ch;	
					int add = 0;	
					int var = 0;	
					while(isdigit(ch)){
						add = atoi(&ch); 
						if (var){
							var = var*10 + add;	
						}else{
							var = add;
						}		
						myin >> ch; // edge
						std::cout << ch;		
					}	
					std::cout << std::endl;	
					//var = var + static_cast<unsigned int> ( (n1->eparent)->edgetoparent);
					//std::cout << ">>pai :" << (n1->eparent)->label << "["<< static_cast<unsigned int> ( (n1->eparent)->edgetoparent)<<"]" <<
					//				" filho :  " << static_cast<unsigned int>(n1->edgetoparent) << std::endl;			
					lgt = var;
					n1->edgetoparent = lgt;					
				}else{
					std::cout << ch<<  " *******hgs gfhdsfg******************hdfghd***** dfghdf*************************dgh*******fghdfgh*****g******************************************possivel erro. Era pra ler Numero e leu ooutra coisas<<"<<std::endl;
					int var = atoi(&ch); 
					lgt= var;
					n1->edgetoparent = lgt;	
					return NULL;
				}					
			}// end of first children
			
			
			/****** create the second clidren *******/
			printf("second children:\n");				
			std::cout << ch<< std::endl;
			
			edgynode* nr = new edgynode();	
			edgynode* aux = n1;
			while(true){	
				// [A]: read ,	
				if(ch==','){
					printf(">>>\n");	
					myin >> ch;	//=> is , --  read: ( or Nn
					printf("%c\n",ch);
					//[A1] read (
					if (ch=='('){
						printf("next	 children is interal node\n");		
						nr = new edgynode();	
						nr->label = (*lab)++;		
						nr->eparent = pai;
						aux->enext = nr;	
						ch = NEXUS_in_subtree(enp,nr,myin,nameFile,lab);
						// after close a level
						aux = nr;		
						printf("closing a  interal node (saindo NEXUS_in_subtree) (OUTRO) \n");		
						std::cout << "tenho: " << ch<< std::endl;	
					//[A2] read Nn
					}else{
						printf("next children is a N1 node\n");;	
						edgynode* nn = new edgynode();	
						nn->label = (*lab)++;		
						nn->eparent = pai;
						aux->enext = nn;	

						std::cout << ".. checking if N1>9-["<< ch <<  std::endl;		
						std::cout << ch;						
						myin >> ch;	
						std::cout << ch;	
						while(isdigit(ch)){
							myin >> ch; // edge
							std::cout << ch;		
						}	
						std::cout << std::endl;							
						std::cout << ".. geting the level"<<  std::endl;		
						if (ch==':'){
							myin >> ch;	
							std::cout << ".. checking if N1>9 in the level "<<  std::endl;		
							std::cout << ch;	
							int add = 0;	
							int var = 0;	
							while(isdigit(ch)){
								add = atoi(&ch); 
								if (var){
									var = var*10 + add;	
								}else{
									var = add;
								}			
								myin >> ch; // edge
								std::cout << ch;		
							}	
							std::cout << std::endl;	
							//var = var + static_cast<unsigned int> ( (nn->eparent)->edgetoparent);
							//std::cout << ">>pai :" << (nn->eparent)->label << "["<<static_cast<unsigned int> ( (nn->eparent)->edgetoparent)<<"]" <<
							//		" filho : " << static_cast<unsigned int> (nn->edgetoparent) << std::endl;										
							lgt= var;													
							nn->edgetoparent = lgt;					
						}
						aux = nn;	
						printf("get the next one \n");;	
						std::cout << ch<< std::endl;							
					}
				//end of [A]
				}else{
					//[B]: read )	
					if(ch==')'){
						printf("closing a NODE !!!.  \n");		
						myin >> ch; // edge
						std::cout << ch<< std::endl;	
						if (ch==':'){
							myin >> ch; // edge
							std::cout << ch<< std::endl;								
							std::cout << ".. checking if N1>9 in the level: "<<  std::endl;		
							std::cout << ch;	
							int add = 0;	
							int var = 0;	
							while(isdigit(ch)){
								add = atoi(&ch); 
								if (var){
									var = var*10 + add;	
								}else{
									var = add;
								}		
								myin >> ch; // edge
								std::cout << ch;		
							}	
							std::cout << std::endl;	
							//var = var + static_cast<unsigned int> ( (pai->eparent)->edgetoparent);	
							//std::cout << ">>pai : " << (pai->eparent)->label << "["<<  static_cast<unsigned int> ( (pai->eparent)->edgetoparent) <<"]" <<
							//		" filho : " << static_cast<unsigned int> (pai->edgetoparent) << std::endl;										
							lgt= var;
							pai->edgetoparent = lgt;					
						}
						pai->elastchild = aux;
						if(ch==';'){
							printf("find DE TUDO.  \n");			
							//QString NewickFile = "ateste.new";

							std::ofstream var (nameFile.toAscii().constData(), std::ios::out); // replace if exists!
							if (var.is_open()) {
								enp->NEXUS_out(var); // edgytree::NEXUS_out() writes in the simpler Newick format rather than NEXUS
								var.close();
							} else {
								std::cout << "Error writing " << nameFile.toAscii().constData() << " for writing" << std::endl;
								return ";";	
							}								
							return ";";		
						}else{
							printf("find of recursion...  \n");		
							return ch;
						}		
					}else{
						printf("ERROR: the input is not ESPERADO. Leu outra coisa que na era ) ,  (\n");
						return NULL;
					}
				}//end of [B]
			}
  } // --edgynode::NEXUS_out_subtree()


  /** receive a buffer in nexus format and create a edgytree.*/
  void NEXUS_in(edgytree *enp, std::istream & myin,QString nameFile) {    
			eroot = new edgynode;
			char ch;
			printf("NEXUS_in\n");
			int lab=1;
			if(myin){
				myin >> ch;
				std::cout << ch<< std::endl;
				if (ch=='(' ){
					eroot->edgetoparent = 0;	
					eroot->label= lab++;	
					enp->eroot = eroot;	
					NEXUS_in_subtree(enp,eroot,myin,nameFile,&lab);
				}else{
					std::cerr << "MALFORMED newick file." << std::endl;
				}
			}else{
				std::cout << "error";	
			}
			//QString NewickFile = "ateste.new";

			//~ std::ofstream var (nameFile, std::ios::out); // replace if exists!
			//~ if (var.is_open()) {
				//~ enp->NEXUS_out(var); // edgytree::NEXUS_out() writes in the simpler Newick format rather than NEXUS
				//~ var.close();
			//~ } else {
				//~ std::cout << "Error opening " <<nameFile << " for writing" << std::endl;
				//~ return;
			//~ }			
  } // --edgynode::NEXUS_out()


  /** flatten this "edgy" subtree with root pointed by enp 
      post: (*enp) still points to same object
  */
  void flatten_subtree(edgynode* enp) {
    // first, iteratively flatten this one line-graph subgraph of subtree, if at all
    if(enp == 0) {
      std::cerr << "enp is null in flatten_subtree()" << std::endl;
      return;
    }
    while((enp->efirstchild) == (enp->elastchild)) {
      if(0 == (enp->efirstchild)) break; // leaf: done
      // number of children == 1.  "absorb" child, don't lose next.
      edgynode* chp = enp->efirstchild;
      (enp->edgetoparent) += (chp->edgetoparent);  // increment edgelength by that much
      (enp->efirstchild) = (chp->efirstchild);    // adopt grandchildren
      (enp->elastchild) = (chp->elastchild);
      (chp->enext) = 0;
      (chp->eparent) = 0;
      (chp->efirstchild) = 0;
      (chp->elastchild) = 0; 
      // prepared child for destruction
      delete chp; // free mem
    }
    // else/finally have not exactly 1 child
    edgynode* chp = enp->efirstchild;
    while(0 != chp) { // if enp is not a leaf,
      chp->eparent = enp;
      flatten_subtree(chp); // recurse!
      chp = chp->enext;
    }
  } // --edgytree::flatten_subtree()

  /** removes all non-root nodes of degree 2
      while maintaining appropriate height information
      in addition, if root has only one child, makes that child the new root and so on
  */
  void flattenDontKeepRoot() {
    if(0 == eroot) {
      std::cerr << "In edgytree::flatten() eroot is null!" << std::endl;
      return;
    }
    while((eroot->efirstchild) == (eroot->elastchild)) {
      if(0 == (eroot->efirstchild)) break; // leaf: done
      edgynode* oldroot = eroot;
      eroot = eroot->efirstchild;
      oldroot->efirstchild = 0;
      oldroot->elastchild = 0;
      delete oldroot;
    }
    eroot->eparent = 0;
    eroot->edgetoparent = 0;
    assert((eroot->enext) == 0);
    // now eroot cannot have exactly 1 child
    edgynode* enp = eroot->efirstchild;
    while(enp != 0) {
      enp->eparent = eroot;
      flatten_subtree(enp);
      enp = enp->enext;
    }
  } // --edgytree::flattenDontKeepRoot()

  /** removes all non-root nodes of degree 2
      while maintaining appropriate height information */
  void flatten() {
    // don't worry about it if root has or ends up having only one child.
    if(0 != eroot) {
      edgynode* enp = eroot->efirstchild;
      while(0 != enp) {
	flatten_subtree(enp);
	enp = enp->enext;
      }
      //      flatten_subtree(eroot); // might leave root with a positive "edgetoparent" value!
      // (which NEXUS wouldn't mind, but we should be cautious)
    } else {
      std::cerr << "In edgytree::flatten() eroot is null!" << std::endl;
    }
  } // --edgytree::flattenKeepRoot()

  /** deep destructor for edgytree */
  ~edgytree() {
    delete eroot;
  } // --edgytree::~edgytree() dtor


  
   /** This function implements a BST to visit each leaf in the tree.
   * pre: enp points to a valid edgynode                              *
   * post: will CHANGE edgy subtree irreversibly (make it pruned)      *
   * returns height of subtree BEFORE the pruning, as needed by callee */
unsigned int pruneSubtreeByOmega(edgynode* enp) {
  if(0==enp) std::cerr << "error in edgytree::pruneSubtreeByOmega: enp is NULL!" << std::endl;
    
	std::queue<edgynode*> myNodes;  
	myNodes.push(enp);
	edgynode* oneedge;
	 
	while(!myNodes.empty()) {  
		oneedge = myNodes.front();
		std::cout << "(*)-[A]  " << oneedge->edgetoparent <<	 std::endl;
		myNodes.pop();
		//checking if the node is leaf
		if (oneedge->efirstchild == 0){
		  	std::cout << "(*)-[B]- <no> " <<  oneedge->edgetoparent  <<	 std::endl;
	  		//edgynode* pa = findProperAncestor();
			//had a properAncestor
			//if (pa != 0){	
				//removedPath(oneedge,pa);
			//}
		}else{
			edgynode* chp = oneedge->efirstchild;	
			std::cout << "(*)-[C]- <pai><filho>  " << oneedge->edgetoparent <<  ", " <<  chp->edgetoparent  <<	 std::endl;
			myNodes.push(chp);
			edgynode* nextn = chp->enext;
			while (nextn!=0){
				myNodes.push(nextn);
				std::cout << "(*)-[D]- <next> " <<  nextn->edgetoparent <<	 std::endl;
				chp = nextn->enext;
				nextn = chp;
			}
		}
		
	}
    return 1; // returns the length of the subtree of parent that goes through "this"
} // --edgytree::pruneByHeight

  /**  pruning by height of edgytree */
  void pruneByOmega(long long omega, unsigned int pa) {
    if(eroot == 0) {
      std::cerr << "eroot is null in edgytree::pruneByHeight()" << std::endl;
      return;
    }
    this->omega = omega;
		this->properAncestor = pa;  
    //    std::cout << "in edgytree::pruneByHeight(), this->minheight is " 
    //	      << (static_cast<unsigned>(this->minheight)) << std::endl;
    //    setparents();
    flatten();
    long long rootsNewHeight = pruneSubtreeByOmega(eroot);
    std::cout << "edgytree::pruneByOmega() completes execution with root's height="
	      << rootsNewHeight << std::endl;
    flatten();
    //    setparents();
  } // --edgytree::pruneByHeight()



/*********************************************************************
THIS NEXT FUNCTION IS FOR TEST THE ISOMORPHISM BETWEEN TREE.
**********************************************************************/

void updateCodeIsomorphism(edgynode* enp) {
		
    if(enp==0) {
      std::cerr << "null ptr passed to edgytree::traverse_subtree()" << std::endl;
      throw std::invalid_argument("null ptr passed to edgytree::traverse_subtree()");
    }

    edgynode *echp = enp->efirstchild;
    if(echp) {
      while(echp != 0) {
				updateCodeIsomorphism(echp);
				echp = echp->enext;
    	}
    }else{
			// set the code to leave					
			enp->code.push_back(1);
		}
		
		
		// set the code to internal nodes		
		if(enp->efirstchild!=0){
			std::map<int, std::vector<int> > list_codes;
			list_codes.clear();	
			edgynode *ef = enp->efirstchild;
			int ord =0;
			int codePar=1;
			cout << " VERTEX (" << enp->label <<")"<< endl;	
			// get code for his children	
			if(ef) {
				while(ef!= 0) {
					codePar += ef->code.front();
					cout << "		filhos (" << ef->label <<") with size "<< ef->code.size() << endl;
					list_codes.insert( pair<int,std::vector<int> >(ord,ef->code));
					ef = ef->enext;
					ord++;
				}
				// set parent code 	
				enp->code.push_back(codePar);
				cout << "<1> " << endl; 
				for (std::multimap<int, std::vector<int> >::iterator it = list_codes.begin(); it!= list_codes.end() ; it++) {
					std::vector<int> aux = 	(*it).second;
					 for (int i=0;i<aux.size();i++)  cout << " " << aux[i]; cout << endl; 	
				}	
					
				// set ordered children code	
				std::vector<int> child_codes =  radix_sort2(&list_codes);
				vector<int>::iterator it;	
				for (int i=0;i<child_codes.size();i++){
						enp->code.push_back(child_codes[i]);	
				}	
				cout << "<2> " << enp->label; for (int i=0;i<enp->code.size();i++)  cout << " " << enp->code[i]; cout << endl;	
				cout << "------------------------------------------------------"<< endl;		
					
			}				
		}
}

std::vector<int>  radix_sort2(std::map<int, std::vector<int> >* list_codes){
	std::multimap<int, std::vector<int> >::iterator it;
	std::multimap<int, std::vector<int> >::iterator xt;	
	int item = 0;
	cout << " ----- radix_sort2 -----"<< endl;		
	// get the vector max size in list_codes
	std::vector<int> aInt;	
	for (it = list_codes->begin(); it!= list_codes->end() ; it++) {
		std::vector<int> aux = 	(*it).second;
		aInt.push_back(a2int(aux, aux.size()));
	}
	cout << " *antes  (int) "; for (int i=0;i<aInt.size();i++)  cout << " " << aInt[i]; cout << endl;		
	// ordering
	for (int i=0;i< (aInt.size()-1);i++){
		for (int j=i; j< aInt.size();j++){
			if (aInt[i]>aInt[j]){
				int aux = aInt[j];
				aInt[j]	=aInt[i];
				aInt[i]	= aux;	
			}
		}
	}
	cout << " *depois  (int) "; for (int i=0;i<aInt.size();i++)  cout << " " << aInt[i]; cout << endl;		
	// setting the return vector	
	std::vector<int> ret;	
	for (int i=0;i<aInt.size();i++){				
		std::vector<int> ref = int2a(aInt[i]);			
		for (int i=0;i<ref.size();i++){
			ret.push_back(ref[i]);
			cout << " put "	<< ref[i] ;
		}
		cout << endl;	
	}
	return ret;
}

int a2int(std::vector<int> arr, int sz) {
	int ret = 0;
	for (int i = 0; i < sz; ++i) {
	  ret *= 10;
	  ret += arr[i];
 	}
	return ret;
}

std::vector<int> int2a(int myInt) {
	std::vector<int> ret;	
	std::vector<int> fin;		
	while (myInt){
		int aux = (int) floor(myInt/10);
		ret.push_back(abs(myInt-aux*10));
		myInt =aux;	
	}
	for (int i = ret.size()-1; i>=0; --i) {
			fin.push_back(ret[i]);
	}
	return fin;
}

std::vector<int>  radix_sort(std::map<int, std::vector<int> >* list_codes){

	std::multimap<int, std::vector<int> >::iterator it;
	std::multimap<int, std::vector<int> >::iterator xt;	
	int max_size_list = 0;
		
	// get the vector max size in list_codes
	for (it = list_codes->begin(); it!= list_codes->end() ; it++) {
		std::vector<int> aux = 	(*it).second;
		if ( aux.size() > max_size_list)
				max_size_list	= aux.size() ;
	}
	// ordering
	for (int i=0;i< max_size_list;i++){
		for (it = list_codes->begin(); it!= list_codes->end() ; it++) { //pega o vector referencia
			std::vector<int> ref = (*it).second;	
			if (ref.size()>i){	
				for (xt = list_codes->begin(); xt!= list_codes->end() ; xt++) {// pecorre os outros vetore para comparar
					std::vector<int> comp = (*xt).second;
					if (ref.size()>i) break;	
					if ( (comp.size()<i) || (comp[i]<ref[i]) ){
						std::vector<int> au  = (*it).second;
						(*it).second = comp;
						(*xt).second = au;
					}
				}
			}
		}				
	}
	// setting the return vector	
	std::vector<int> ret;	
	for (it = list_codes->begin(); it!= list_codes->end() ; it++) {	
		std::vector<int> ref = (*it).second;	
		cout << "ret in radix. Size "<< ref.size();
		for (int i=0;i<ref.size();i++){
			ret.push_back(ref[i]);
			cout << " put "	<< ref[i] ;
		}
		cout << endl;	
	}
	return ret;
			
}

/* Set the edgetoparent of each leave zero */
 int  setLeavesZero(edgynode* root) {
    std::queue<edgynode*> s;
		printf(" FUNCTION  setLeavesZero \n");  
	  
		s.push(root);	 
		long long h=0;
    while(!(s.empty())) {
			edgynode* np_; // node pointer
			edgynode* chp_; // child pointer
			np_ = s.front();
			s.pop();
			if (np_->efirstchild!=0){
					chp_ = (np_->efirstchild);
					std::cout << "(*)-[First C]  [ " << (chp_->label ) << "," << static_cast<int> (chp_->edgetoparent )<<"]" <<std::endl;		
					s.push(chp_);		

					while(chp_->enext != 0) {				
						edgynode* aux = chp_->enext;				
						if ( aux != np_){	
							std::cout << "(*)-[Next C]  [" <<  (aux->label ) << "," << static_cast<int> (aux->edgetoparent )<<"]" <<std::endl;
							s.push(aux);	
						}
						chp_ = aux;	
					}						
			}else{
				long long val = 	np_->edgetoparent;
				std::cout << "(*)-[NODE]   [" << (np_	->label ) << " ," << static_cast<int> (np_->edgetoparent )<<"]" << "val:" << val <<std::endl;	
				if (val > h ) h = val;
				np_->isomorphic = 0;
				//np_->edgetoparent = static_cast<unsigned char> (0);								
			}
	}
	std::cout << "(*)-[H :"<< h<<std::endl;		 
	return h;
}


/* Set the edgetoparent of each node for a real distance between root and node*/
 int  setLabels() {
		printf(" FUNCTION setLabels \n"); 
    std::queue<edgynode*> s;
				//std::cout << "(*)-[setLeavesZero] <no> [root:" <<  (eroot->label ) << "," << static_cast<int> (eroot->edgetoparent )<<"]" <<std::endl;
	  
		s.push(eroot);	 
		int h=0;
    while(!(s.empty())) {
			edgynode* np_; // node pointer
			edgynode* chp_; // child pointer
			np_ = s.front();
			s.pop();
				
			
			if ( np_->eparent!=0){
						//int var = static_cast<unsigned int> ((np_->eparent)->edgetoparent) + static_cast<unsigned int> (np_->edgetoparent);						
						int var =  ((np_->eparent)->edgetoparent) + 1;						
						std::cout << ">> >>>  pai : [" << (np_->eparent)->label << " , "<<   (np_->eparent)->edgetoparent  <<"]" <<
									" eu : ["<< np_->label << " , "<< np_->edgetoparent << "]" << " var:" << var << std::endl;	
						np_->edgetoparent = 	var;	
						
			}
				
			if (np_->efirstchild!=0){
					chp_ = (np_->efirstchild);
					//std::cout << "(*)-[First C]  [ " << (chp_->label ) << "," << static_cast<int> (chp_->edgetoparent )<<"]"<< std::endl;	
					s.push(chp_);		

					while(chp_->enext != 0) {				
						edgynode* aux = chp_->enext;				
						if ( aux != np_){	
							//std::cout << "(*)-[Next C]  [" <<  (aux->label ) << "," << static_cast<int> (aux->edgetoparent )<<"]" << std::endl;	
							s.push(aux);	
						}
						chp_ = aux;	
					}						
			}
		}
		std::cout << "(*)-[H :"<< h<<std::endl;		 
		return h;
}


std::list<edgynode*>  getLeavesInLevel(int h){
		printf(" FUNCTION getLeavesInLevel \n");
		std::list<edgynode*>  list;
		std::queue<edgynode*> s;
	  
		s.push(eroot);	 
    while(!(s.empty())) {
			edgynode* np_; // node pointer
			edgynode* chp_; // child pointer
			np_ = s.front();				
			s.pop();
			if (np_->efirstchild!=0){
				chp_ = (np_->efirstchild);
				std::cout << "(*)-[First C]  [ " << (chp_->label ) << "," << static_cast<int> (chp_->edgetoparent )<<"]" <<std::endl;		
				s.push(chp_);		

				while(chp_->enext != 0) {				
					edgynode* aux = chp_->enext;				
					if ( aux != np_){	
						std::cout << "(*)-[Next C]  [" <<  (aux->label ) << "," << static_cast<int> (aux->edgetoparent )<<"]" <<std::endl;
						s.push(aux);	
					}
					chp_ = aux;	
				}						
			}else{
				std::cout << "(*)-[NO ]  [" <<  (np_->label ) << "," << static_cast<int> (np_->edgetoparent )<<" , " << np_->isomorphic <<  "]" <<std::endl;	
				if ( (static_cast<int> (np_->edgetoparent )) == h ){
					std::cout << "h:" << h<< " iso:"<< (static_cast<int> (np_->isomorphic))<< std::endl;
					list.push_back(np_);
				}
			}
		}
		std::cout << "h:" << h<< " size list:"<< list.size()<< std::endl;
		return list;
}

void builTuple( std::list<edgynode*>* list, std::multimap<int, int>* mymultimap){		
	printf(" FUNCTION builTuple \n");	
	int pai = 300;
	int filho;
	while(!(list->empty())) {
		edgynode* no = list->back();
		filho = no->isomorphic;
		edgynode* np = no->eparent;	
		pai = static_cast<int> (np->edgetoparent);
		mymultimap->insert ( pair<int,int>(pai,filho) );
		list->pop_back();	
	}
}



int distintcElements(std::multimap<int, int>* map){
	printf(" FUNCTION distintcElements \n");
	
	std::multimap<int,int>::iterator it;
	int f;
	int s;
	int count =0;
	int dist= 1;	
	for (it = map->begin(); it!= map->end() ; it++) {
		std::cout << (*it).first << " => " << (*it).second << std::endl;
		if (count==0){
			f = (*it).first;
			s = (*it).second;		
			count++;
		}else{
			if (  (f!=(*it).first)&& (s!=(*it).second) ){
				dist++;
			}
			count++;
		}
	}
	std::cout << " distintcElements==" << dist	<< std::endl;	
	return dist;	
}
/******************************************************************************
* Funtion that traverse a tree and show the gray and size of each node
******************************************************************************/
 void testando(edgynode* up) {
    std::stack<edgynode*> s;
		if(bv::debug) printf("\n (*)-[teste EDGE] <root> [%d] \n",up->label);	 
	 //if(bv::debug)	std::cout << "(*)-[teste] <no> [pnc:" << static_cast<int> (up->gray ) << "," << static_cast<int> (up->size )<<"]" <<std::endl;
	 
		long long a =0;
		long long b = 0;
	  s.push(up);	 
		cout << "root: "<< up->label << " It contains:";  for (int i=0;i<up->code.size();i++)  cout << " " << up->code[i]; cout << endl; 
    while(!(s.empty())) {
			edgynode* np_; // node pointer
			edgynode* chp_; // child pointer
			np_ = s.top();	

			std::cout << "(*)-[A] <no> [pnc:" << static_cast<int> (np_->edgetoparent ) << "" <<std::endl;	
			s.pop();
			chp_ = (np_->efirstchild);
			while(chp_ != 0) { // and add ptrs to all children to the stack
				s.push(chp_);		
				std::cout << "(*)-[B] <no> [pnc:" << static_cast<int> (np_->edgetoparent ) << "" <<std::endl;		
				chp_ = (chp_->enext);
			}
		}

} // --historytree::testando()


}; // -- struct edgytree
#endif
