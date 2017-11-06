#include <QtGui>
#include "TreeManager.h"
#include "TreeScene.h"
#include "Node.h"
#include "Edge.h"
#include <QDebug>

using namespace std;

TreeManager::   TreeManager() {
    const int defaultSceneWidth = 700;
    const int defaultSceneHeight = 400;
    m_scene = new TreeScene(-defaultSceneWidth / 2, -Node::maxSize().height(),
                            defaultSceneWidth, defaultSceneHeight);
    connect(m_scene, SIGNAL(itemsMoving()), this, SLOT(updateSceneRect()));
    m_scene->setBackgroundBrush(QColor(Qt::white));
    m_root = 0;
    m_treeSize = 0;
    showBox = true;
    showLabel = true;
    showVoxel = false;
    showCV = false;
    connect(this, SIGNAL(treeChanged()), this, SLOT(updateScene()));
}

TreeScene* TreeManager::scene() {
    return m_scene;
}

int TreeManager::treeSize() {
    return m_treeSize;
}

void TreeManager::setShowBox(bool b){
    showBox =b;
}

void TreeManager::setShowLabel(bool b){
    showLabel=b;
}

bool TreeManager::getShowLabel(){
    return showLabel;
}

bool TreeManager::getShowBox(){
    return showBox;
}

void TreeManager::setShowVoxel(bool b){
    showVoxel=b;
}

bool TreeManager::getShowVoxel(){
    return showVoxel;
}

void TreeManager::setShowCV(bool b){
    showCV=b;
}

bool TreeManager::getShowCV(){
    return showCV;
}

bool TreeManager::parseTree(historytree *myTree, int lambda){
    if(bv::debug) printf(" TreeManager::parseTree() \n");
    deleteNode(m_root);
    int cc=0;
    QTime myTimer;
    myTimer.start();

    m_root = new Node(myTree->hroot);
    m_scene->addItem(m_root);
    m_scene->addItem(m_root->parentEdge());
    m_scene->addItem(m_root->horizontalEdge());

    Node* enp;
    Node* echp;
    std::stack<historynode*> hns; // history node stack
    std::stack<Node*> ens;  // edgy node stack
    if (lambda!=0){
        if(bv::debug) printf(" parsing simplified tree \n");
        historynode* hnp;
        historynode* hchp;
        historynode* r = myTree->hroot;

        if(bv::debug) printf("root=%lld (%lld,%lld) \n",r->gray,r->deltaLimit,myTree->delta1);
        if (myTree->isCriticalVertex(myTree->hroot)){
            historynode* m = r;
            while (	(r->deltaLimit <= myTree->delta1) && (r->hfirstchild!=0) ){
                 if(bv::debug) printf(".. %lld (%lld,%lld)\n",r->gray,r->deltaLimit,myTree->delta1);
                r = r->successor;
            }
        }else{
            r = myTree->immediateCriticalVertice(myTree->hroot);
            while (	(r->deltaLimit <= myTree->delta1) && (r->hfirstchild!=0) ){
                if(bv::debug) printf(",,  %lld (%lld,%lld)\n",r->gray,r->deltaLimit,myTree->delta1);
                r = r->successor;
            }             
        }
        hns.push(r);
        Node* n_icv = new Node(r);
        if(bv::debug) printf("r=%lld \n",r->gray);
        m_root->addChild(n_icv);
        m_scene->addItem(n_icv);
        m_scene->addItem(n_icv->parentEdge());
        m_scene->addItem(n_icv->horizontalEdge());

        ens.push(n_icv);
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
            m_treeSize++;
            if(bv::debug) printf("$--$ hnp=%lld %d \n",hnp->gray,hnp->listCoord.size());
            if(hnp->sizeNewIPCD!=0) {
                Node* epp = 0; // ptr to edgy "predecessor" (left-sibling) node
                for (int i=0; i< hnp->sizeNewIPCD;i++){ // and add ptrs to selected children to stacks
                     hchp = (hnp->newIPCD[i]);
                     hns.push(hchp);
                     echp = new Node(hchp);
                     enp->addChild(echp);
                     m_scene->addItem(echp);
                     m_scene->addItem(echp->parentEdge());
                     m_scene->addItem(echp->horizontalEdge());
                     ens.push(echp);
                     epp = echp;
                }
            }else{
                enp->setMeLeaf();
            }
         }
    }else{
        if(bv::debug) printf(" parsing usimplified tree \n");
         hns.push(myTree->hroot);
         ens.push(m_root);
         historynode* np_; // node pointer
         historynode* chp_; // child pointer
         Node* enp;
         Node* echp;

         while(!(hns.empty())) {
             cc++;
             np_ = hns.top();
             hns.pop();
             enp = ens.top();
             ens.pop();
             m_treeSize++;
             historynode* aux;
             if (np_->hfirstchild!=0){
                 chp_ = (np_->hfirstchild);
                 if(bv::debug) printf(" FF child\n");
                 if (getShowCV()){
                    aux = myTree->immediateCriticalVertice(chp_);
                 }else{
                    aux = chp_;
                 }

                 echp = new Node(aux);
                 enp->addChild(echp);
                 m_scene->addItem(echp);
                 m_scene->addItem(echp->parentEdge());
                 m_scene->addItem(echp->horizontalEdge());
                 ens.push(echp);
                 hns.push(aux);

                 historynode* pnc_= chp_->hnext;
                 while(pnc_!= 0) {
                     if(bv::debug) printf(" SS child\n");

                     if (getShowCV()){
                        aux = myTree->immediateCriticalVertice(pnc_);
                     }else{
                        aux =pnc_;
                     }
                     hns.push(aux);
                     Node* ec_ = new Node(aux);
                     enp->addChild(ec_);
                     m_scene->addItem(ec_);
                     m_scene->addItem(ec_->parentEdge());
                     m_scene->addItem(ec_->horizontalEdge());
                     ens.push(ec_);

                     chp_ = pnc_->hnext;
                     pnc_ = chp_;
                 }
             }
         }

    }
    cout << " .....parseTree: A [" << myTimer.elapsed() <<" ms]  cc=" <<cc << std::endl;
    int tim = myTimer.elapsed();
    emit treeChanged();
    cout << " .....treeChanged: A [" << (myTimer.elapsed() -tim) <<" ms]" << std::endl;
    return true;
}

bool TreeManager::removeValue(int value) {
    Node *node = searchValue(value);
    if (node) {
        deleteNode(node);
        return true;
    } else
        return false;
}

bool TreeManager::deleteSelected() {
    QList<QGraphicsItem *> selection = m_scene->selectedItems();
    if (selection.isEmpty())
        return false;

    foreach(QGraphicsItem *item, selection) {
        Node *node = dynamic_cast<Node *>(item);
        if (node) deleteNode(node);
    }

    return true;
}

historynode* TreeManager::getHistoryNodeSelected() {
    historynode* n;
    QList<QGraphicsItem *> selection = m_scene->selectedItems();
    if (selection.isEmpty()){
        return 0;
    }else{
        if (selection.size()!=1){
            return 0;
        }
    }

    foreach(QGraphicsItem *item, selection) {
        Node *node = dynamic_cast<Node *>(item);
        n = node->getHistoryNode();
    }

    return n ;
}

Node* TreeManager::getNodeSelected() {
    historynode* n;
    QList<QGraphicsItem *> selection = m_scene->selectedItems();
    if (selection.isEmpty()){
        return 0;
    }else{
        if (selection.size()!=1){
            return 0;
        }
    }
    Node *node;
    foreach(QGraphicsItem *item, selection) {
        node = dynamic_cast<Node *>(item);
    }

    return node;
}

void TreeManager::setRoot(Node *node) {
    m_root = node;
    node->removeParentEdge();
}

Node* TreeManager::searchValue(int value) {
    if (m_root)
        return searchDFS(m_root, value);
    else return 0;
}
Node* TreeManager::searchDFS(Node *cur, int val) {
/*    if (cur->value() == val)
        return cur;

    if (val < cur->value()) {
        if (Node *left = cur->leftChild())
            return searchDFS(left, val);
        else return 0;
    }

    if (cur->value() < val) {
        if (Node *right = cur->rightChild())
            return searchDFS(right, val);
        else return 0;
    }
*/
    return 0;
}

/*void TreeManager::deleteNode(Node *node) {
    if (!node) return;

    Node* p = node->parent();
    Node* l = node->leftChild();
    Node* r = node->rightChild();
    int v = node->value();
    bool root = (node == m_root);

    m_scene->removeItem(node);
    if (node->parentEdge())
        m_scene->removeItem(node->parentEdge());
    delete node;

    Node *lowestNode = 0;       //subTreeSize update start node

    if (!l && !r) {             //node has no children
        if (root)
            m_root = 0;
        else
            lowestNode = p;     //lowest = node's parent
    } else {                    //node has children
        Node *x;
        if (l && r) {           //node has both children
            x = r;
            if (x->leftChild()) {       //node's right child has left child
                while (x->leftChild())  //x = the most left node in node's right child subtree
                    x = x->leftChild();

                lowestNode = x->parent();   //lowest = parent of the x node

                if (x->rightChild()) {      //x has no left child, but maybe it has right
                    x->parent()->setLeftChild(x->rightChild());
                }

                x->setRightChild(r);
            } else {
                lowestNode = x;     //node's right child has no left child; lowest = node's right child
            }
            x->setLeftChild(l);
        } else {                //node has only one child
            x = l ? l : r;
            lowestNode = x;     //lowest = node's child which takes node's stand
        }

        if (root) {
            setRoot(x);
        } else {
            if (v < p->value())
                p->setLeftChild(x);
            else
                p->setRightChild(x);

        }
    }

    if (lowestNode)
        lowestNode->updateSubTreeSize();
    m_treeSize--;
    emit treeChanged();

}*/




int TreeManager::subTreeHeight(Node *node) {
    return subTreeHeightDFS(node);
}
int TreeManager::subTreeHeightDFS(Node *cur) {
    int l = 0, r = 0;
/*    if (Node *left = cur->leftChild())
        l = subTreeHeightDFS(left);
    if (Node *right = cur->rightChild())
        r = subTreeHeightDFS(right);
*/
    return qMax(l, r) + 1;
}

int TreeManager::subTreeSlots(Node *node) {
    return subTreeSlotsDFS(node);
}

int TreeManager::subTreeSlotsDFS(Node *cur) {
    if (cur->getFistChild()){//non leave
       if (cur->getFistChild()==cur->getLastChild()){//no CV
           cur->setTreeSlots(subTreeSlotsDFS(cur->getFistChild()));
           return  cur->subTreeSlots() ;
       }else{//is CV
           Node *chp_ = cur->getFistChild();
           int c = 0;
           c = subTreeSlotsDFS(chp_); // first child
           Node* pnc_ = chp_->getNextChild();
           while(pnc_!= 0) {
               c = c + subTreeSlotsDFS(pnc_);// others child
               chp_ = pnc_->getNextChild();
               pnc_ = chp_;
           }
           cur->setTreeSlots(c);
           return c;
       }
    }else{//leaf
        cur->setTreeSlots(1);
        return 1;
    }
}

void TreeManager::printTree() {
     std::stack<Node*> s;
     //printf("(***************************)-[printTree] \n");
     s.push(m_root);
     int mnv = m_root->getHistoryNode()->numberVoxel;
     m_root->setShowBox(showBox);
     m_root->setShowLabel(showLabel);
     m_root->setShowVoxel(showVoxel);
     while(!(s.empty())) {
         Node* np_; // node pointer
         Node* chp_; // child pointer
         np_ = s.top();
         np_->setMaximumNumberVoxels(mnv);
         //printf(" Node (%d %d  %d) \n",np_->value(),np_->subTreeSlots(),np_->subTreeSize());
         s.pop();
         if (m_root!=np_){
             np_->plotNewickEdge();
             np_->setShowBox(showBox);
             np_->setShowLabel(showLabel);
             np_->setShowVoxel(showVoxel);
         }else{
             np_->plotNewickEdgeRoot();
             np_->setShowBox(showBox);
             np_->setShowLabel(showLabel);
             np_->setShowVoxel(showVoxel);
         }
         if (np_->getFistChild()!=0){
             chp_ = (np_->getFistChild());
             //printf("< [FC] (%d %d %d)  ",chp_->value(),chp_->subTreeSlots(),chp_->subTreeSize());
             s.push(chp_);
             Node* pnc_ = chp_->getNextChild();
             while(pnc_!= 0) {
                 s.push(pnc_);
                 //printf("[NC] (%d %d %d)",pnc_->value(),pnc_->subTreeSlots(),pnc_->subTreeSize());
                 chp_ = pnc_->getNextChild();
                 pnc_ = chp_;
             }
             //printf(">\n");
         }else{
             //printf(" [L] (%d %d  %d)\n",np_->value(),np_->subTreeSlots(),np_->subTreeSize());
         }
        //printf(" ]\n");
     }
}
void TreeManager::updateScene() {
    if(bv::debug) printf("updateScene \n");
    subTreeSlots(m_root);
    recalcRelPos_FLAT();
    recalcAbsPos_FLAT();
    updateSceneRect();
    printTree();
}

void TreeManager::updateSceneRect() {
    if(bv::debug) printf("updateSceneRect \n");
    qreal maxY = 0;
    qreal maxX = 0;

    QList<QGraphicsItem *> items = m_scene->items();

    foreach(QGraphicsItem *item, items) {
        maxX = qMax(maxX, qAbs(item->sceneBoundingRect().right()));
        maxX = qMax(maxX, qAbs(item->sceneBoundingRect().left()));
        maxY = qMax(maxY, qAbs(item->sceneBoundingRect().bottom()));
    }

    maxX += Node::maxSize().width() / 2;
    maxY += 2 * Node::maxSize().height();

    m_scene->setSceneRect(-maxX, -Node::maxSize().height(), 2 * maxX, maxY);
}

void TreeManager::recalcRelPos() {
    if(bv::debug) printf("recalcRelPos \n");
    if (!m_root) return;
    qreal l = m_root->boundingRect().width() / 2;
    qreal r = m_root->boundingRect().width() / 2;
    m_root->setPos((l - r) / 2, 0);
}

qreal TreeManager::recalcRelPosDFS(Node *node, ChildType type) {
    const qreal margin = 5;
    qreal leftSubTreeWidth = 0, rightSubTreeWidth = 0;

    return leftSubTreeWidth + rightSubTreeWidth;
}

void TreeManager::recalcRelPos_FLAT(){
    recalcRelPosDFS_FLAT(m_root);
}

void TreeManager::recalcRelPosDFS_FLAT(Node *node) {
    const qreal margin = 2;
    int unit = 10;
    int val=0;
    if (node->getFistChild()){ //non leave
       if (node->getFistChild()==node->getLastChild()){//no CV
           Node *np_ = node->getFistChild();
           recalcRelPosDFS_FLAT(np_);
           if (np_!=m_root){
               if (getShowCV()){
                   //val = 1;
                    val =np_->getHistoryNode()->gray- np_->parent()->getHistoryNode()->gray;
               }else{
                   val =np_->getHistoryNode()->gray- np_->parent()->getHistoryNode()->gray;
               }
           }
           np_->setRelPos(QPointF(0,val*unit));
           if(bv::debug) printf("[[[[HC]]]] (%d)  \n",val);
       }else{//is CV
           qreal dv = 0; //width allocated to  node
           qreal dc = 0; //width allocated to a child of node
           qreal pc = 0; //relative position of a child of node
           qreal g = 0;
           Node *chp_ = node->getFistChild();
           //if(bv::debug) printf("[[[[FC]]]] (%d)  \n",val);
           recalcRelPosDFS_FLAT(chp_); // first child
           Node* pnc_ = chp_->getNextChild();
           while(pnc_!= 0) {
               //if(bv::debug) printf("[[[[NC]]]] (%d)  \n",val);
               recalcRelPosDFS_FLAT(pnc_);// others child
               chp_ = pnc_->getNextChild();
               pnc_ = chp_;
           }

           //compute the relative position for the childrem
           dv = (node->boundingRect().width()*node->subTreeSlots()) + (node->subTreeSlots()-1)*margin;
           chp_ = node->getFistChild();
           dc = (chp_->boundingRect().width()*chp_->subTreeSlots()) + (chp_->subTreeSlots()-1)*margin;
           pc = -(dv/2) + dc/2;

           if (chp_!=m_root){
               if (getShowCV()){
                   //val = 1;
                   val =chp_->getHistoryNode()->gray- chp_->parent()->getHistoryNode()->gray;
               }else{
                   val =  chp_->getHistoryNode()->gray -chp_->parent()->getHistoryNode()->gray ;
               }
           }

           chp_->setRelPos(QPointF(pc,val*unit));
           if(bv::debug) printf("<FC> (%d) {%f} [%d] \n",chp_->value(),pc,val);
           pnc_ = chp_->getNextChild();
           while(pnc_!= 0) {
               qreal dn = (pnc_->boundingRect().width()*pnc_->subTreeSlots()) + (pnc_->subTreeSlots()-1)*margin;
               qreal pn = pc + dc/2 + margin + dn/  2;
               if (chp_!=m_root){
                   if (getShowCV()){
                       //val = 1;
                        val = pnc_->getHistoryNode()->gray- pnc_->parent()->getHistoryNode()->gray;
                   }else{
                        val = pnc_->getHistoryNode()->gray -pnc_->parent()->getHistoryNode()->gray;
                   }
               }
               pnc_->setRelPos(QPointF(pn,val*unit));
               pc = pn;
               dc = dn;
               chp_ = pnc_->getNextChild();
               pnc_ = chp_;
           }
       }
    }else{//leaf
        if (getShowCV()){
            //val = 1;
            val = node->getHistoryNode()->gray - node->parent()->getHistoryNode()->gray;
        }else{
            val = node->getHistoryNode()->gray - node->parent()->getHistoryNode()->gray;
        }
        node->setRelPos(QPointF(0,val*unit));
    }
}

void TreeManager::recalcAbsPos_FLAT() {
    if (!m_root) return;
    if(bv::debug) printf("..(NODE) (%d)(%f)(%f,%f) \n",m_root->value(),m_root->boundingRect().width(),m_root->pos().x(),m_root->pos().y());
    //m_root->setPos(0, 0);
    if (m_root->getFistChild()){
        Node *chp_ = m_root->getFistChild();
        recalcAbsPosDFS_FLAT(chp_);
        Node* pnc_ = chp_->getNextChild();
        while(pnc_!= 0) {
            recalcAbsPosDFS_FLAT(pnc_);
            chp_ = pnc_->getNextChild();
            pnc_ = chp_;
        }
    }

}
void TreeManager::recalcAbsPosDFS_FLAT(Node *node) {

    if(bv::debug) printf("(NO) (%d)<%d>(%f)(%f,%f)(%f,%f)\n",node->value(),node->subTreeSlots(),node->boundingRect().width(),node->relPos().x(),node->relPos().y(),node->parent()->pos().x(),node->parent()->pos().y());
    node->setPos(node->parent()->pos()+node->relPos()); // first child
    if (node->getFistChild()){
        Node *chp_ = node->getFistChild();
        if(bv::debug) printf("..(FC) (%d)  \n",chp_->value());
        recalcAbsPosDFS_FLAT(chp_);
        Node* pnc_ = chp_->getNextChild();
        while(pnc_!= 0) {
            if(bv::debug) printf("..(NC) (%d)  \n",pnc_->value());
            recalcAbsPosDFS_FLAT(pnc_);// others child
            chp_ = pnc_->getNextChild();
            pnc_ = chp_;
        }
    }
}


void TreeManager::recalcAbsPos() {
    if (!m_root) return;
/*
    //m_root->setPos(0, 0);
    if (Node *left = m_root->leftChild())
        recalcAbsPosDFS(left);
    if (Node *right = m_root->rightChild())
        recalcAbsPosDFS(right);*/
}
void TreeManager::recalcAbsPosDFS(Node *node) {
/*    node->setPos(node->parent()->pos() + node->relPos());
    if (Node *left = node->leftChild())
        recalcAbsPosDFS(left);
    if (Node *right = node->rightChild())
        recalcAbsPosDFS(right);*/
}


bool TreeManager::isEmpty() {
    return m_treeSize == 0;
}

int TreeManager::valueByIndex(int index) {
    if ((index > m_treeSize) || (index < 0) || (!m_root))
        return 0;

    return valueByIndexDFS(m_root, index);
}

int TreeManager::valueByIndexDFS(Node *node, int index) {
/*    int nodeIndex = (node->leftChild() ? node->leftChild()->subTreeSize() : 0) + 1;

    if (index < nodeIndex)
        return valueByIndexDFS(node->leftChild(), index);
    else if (index > nodeIndex)
        return valueByIndexDFS(node->rightChild(), index - nodeIndex);
*/    return node->value();
}

void TreeManager::deleteNode(Node *node) {
    if (!node) return;
    //if(bv::debug) printf(" To delete [%d] ",node->value());
    if (node->getFistChild()){
        Node *chp_ = node->getFistChild();
        Node *pnc_ = chp_->getNextChild();
        //if(bv::debug) printf("    to delete FC (%d) ",chp_->value());
        deleteNode(chp_);
        while(pnc_!= 0) {
            chp_ = pnc_->getNextChild();
            //if(bv::debug) printf("    to delete NC (%d) ",pnc_->value());
            deleteNode(pnc_);
            pnc_ = chp_;
        }
    }

    bool root = (node == m_root);

    m_scene->removeItem(node);
    if (node->parentEdge()){
        m_scene->removeItem(node->parentEdge());
        m_scene->removeItem(node->horizontalEdge());
    }
    if(bv::debug) printf(" deleting [%d] ",node->value());
    delete node;
    if (root)
        m_root = 0;

    m_treeSize--;
    //emit treeChanged();
}

bool TreeManager::plotDescendant(historynode *myNod, historytree *myTree, Node *nod){
     printf(" TreeManager::plotDescendant() \n");

    historynode *hnod = myNod->hparent;
    Node *newNode;
    while (!myTree->isCriticalVertex(hnod)){
        //create a new Node
        newNode = new Node(hnod);
        //get the parent of nod and remove from nod
        Node *pai = nod->parent();
        // add nod into newNode
        newNode->addChild(nod);
        newNode->setParent(pai);
        m_scene->addItem(newNode);
        m_scene->addItem(newNode->parentEdge());
        m_scene->addItem(newNode->horizontalEdge());
        //find the position of the NOD in PAI
        Node * fic = pai->getFistChild();
        if (fic==nod){
            Node * sib = fic->getNextChild();
            pai->setFistChild(newNode);
            newNode->setNextChild(sib);
        }else{
            Node * sib = fic->getNextChild();
            while (sib!=nod){
                fic = sib;
                sib = sib->getNextChild();
            }
            fic->setNextChild(newNode);
            if (sib->getNextChild()!=0){
                newNode->setNextChild(sib->getNextChild());
            }else{
                pai->setLastChild(newNode);
            }

        }
        hnod = hnod->hparent;
        nod = newNode;
    }
    emit treeChanged();
}
