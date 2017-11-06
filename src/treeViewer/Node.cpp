#include <QtGui>
#include "Node.h"
#include "Edge.h"
//#include "../src/CPPHeaders.hpp"
//#include "../src/constants.hpp"
//#include "../src/deftypes.hpp"
//#include "../historytree.hpp"
#include <QDebug>

//typedef int valtype;

using namespace std;

const qreal Node::hoverScaleFactor = 1.1;
const int Node::padding = 2;
const int Node::glowWidth = 4;
const int Node::margin = Node::glowWidth + 1;

Node::Node(historynode *node) {
    histNode = node;
    m_value =  node->gray;
    m_parentEdge = new Edge(this, this);
    m_horizontalEdge = new Edge(this, this);
    //m_leftChild = m_rightChild = 0;
    m_subTreeSize = 1;
    m_next = m_firstchild = m_lastchild=0;

    m_textColor.setRgb(0x744f0e);
    m_outlineColor.setRgb(0x65502c);
    m_backColor.setRgb(0xffdda2);
    m_backGradient.setStart(QPointF(0, -10));
    m_backGradient.setFinalStop(QPointF(0, 30));
    m_backGradient.setColorAt(0, QColor(0xffebc8));
    m_backGradient.setColorAt(1, QColor(0xdaa64b));   
    setFlags(ItemIsMovable | ItemIsSelectable);
    setAcceptsHoverEvents(true);
    voxelNumberMessage = new QGraphicsSimpleTextItem(this);
    m_hover = false;

#if QT_VERSION >= 0x040600
    setFlag(ItemSendsGeometryChanges, true);
#endif
}

Node::~Node() {
/*    if (parent())
        parent()->removeChild(this);

    if (m_leftChild)
        m_leftChild->removeParent();

    if (m_rightChild)
        m_rightChild->removeParent();

    delete m_parentEdge;*/
}

//--------------------------------------------------------------------

int Node::value() const {
    return m_value;
}

QColor Node::backColor() const {
    return m_backColor;
}

void Node::setBackColor(const QColor &color) {
    m_backColor = color;
    update();
}

//--------------------------------------------------------------------

Edge* Node::parentEdge() const {
    return m_parentEdge;
}

Edge* Node::horizontalEdge() const {
    return m_horizontalEdge;
}
Node* Node::parent() const {
    Node* parent = m_parentEdge ? m_parentEdge->fromNode() : 0;
    if (parent == this)
        parent = 0;
    return parent;
}

void Node::setParent(Node *node) {
    if(bv::debug) std::cout << "setParent()  \n";
    if (parent() == node)
        return;
    if (m_parentEdge) {
        m_horizontalEdge->setFromNode(node);
        m_parentEdge->setFromNode(node);
    } else {
        m_parentEdge = new Edge(node, this);
        m_horizontalEdge = new Edge(node, this);
    }
}

void Node::setMeLeaf() {
    m_firstchild = m_lastchild = 0;
}

void Node::addChild(Node* chp) {
    if(bv::debug) std::cout << "addChild()  \n";
    if(chp==0) {
        std::cerr << "error: attempt to pass null pointer as child!" << std::endl;
        throw std::invalid_argument("tried to add null child");
        return;
    }

    if(m_lastchild) { // add new sibling
        if(bv::debug) std::cout << " add new sibling .. \n";
        m_lastchild->m_next = chp;
        m_lastchild = chp;
    } else {
        if(bv::debug) std::cout << " add first child.. \n";
        m_firstchild = m_lastchild = chp;
    }
    chp->setParent(this);
    chp->parent()->updateSubTreeSize();
} // --addChild()


void Node::removeParent() {
    if (m_parentEdge){
        m_parentEdge->setFromNode(0);
        m_horizontalEdge->setFromNode(0);
    }
}
void Node::removeChild(Node *node) {   
/*    if (m_leftChild == node)
        m_leftChild = 0;
    if (m_rightChild == node)
        m_rightChild = 0;*/
}
void Node::removeParentEdge() {
/*    if (m_parentEdge) {
        if (parent())
            parent()->removeChild(this);

        m_parentEdge->setFromNode(0);
    }*/
}

QPointF Node::relPos() const {
    return m_relPos;
}

void Node::setRelPos(const QPointF &relPos) {
    m_relPos = relPos;
}

//--------------------------------------------------------------------

QVector<Edge *> Node::edges() {
    QVector<Edge *> edges;
    if (m_parentEdge)
        edges.append(m_parentEdge);
    if (m_firstchild)
        edges.append(m_firstchild->parentEdge());
    Node *aux = m_next;
    while(aux){
        Node *n = aux;
        edges.append(n->parentEdge());
        aux = n->m_next;
    }
    return edges;
}

//--------------------------------------------------------------------

QRectF Node::boundingRect() const {
    return outlineRect().adjusted(-margin, -margin, margin, margin);
}

QPainterPath Node::shape() const {
    QRectF rect = outlineRect();
    QPainterPath path;
    path.addRoundRect(rect, roundness(rect.width()), roundness(rect.height()));

    return path;
}

QRectF Node::outlineRect() const {
    //QFontMetrics fm = qApp->font();
    QFont font("Verdana", 5);
    QFontMetrics fm(font);
    QRectF rect = fm.boundingRect(5);
    //rect = fm.boundingRect(QString::number(value()));
    rect = fm.boundingRect(10);
    rect.adjust(-padding, -padding, padding, padding);
    rect.translate(-rect.center());
    return rect;
}

int Node::roundness(double size) const {
    const int Diameter = 120;
    return 100 * Diameter / int(size);
}

void Node::paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
        QWidget *widget) {

    if (m_hover) {
        QRect r = outlineRect().toRect();
        r.adjust(-glowWidth, -glowWidth, glowWidth, glowWidth);
        QPixmap glow = QPixmap(":/img/nodeglow.png").scaled(r.size(),
                Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
        painter->drawPixmap(-glow.size().width() / 2,
               -glow.size().height() / 2, glow);
    }

    QPen pen(m_outlineColor.darker(350));
    if (option->state & QStyle::State_Selected) {
        pen.setStyle(Qt::DotLine);
        pen.setWidth(2);
    }

    painter->setPen(pen);
    if(showBox){
        if(showVoxel){
            m_backGradient.setStart(QPointF(0, -10));
            m_backGradient.setFinalStop(QPointF(0, 30));
            float range = 255.0;
            float myvalue = (this->getHistoryNode()->numberVoxel*range)/this->getMaximumNumberVoxels();
            int myvalueInt = (int)floor(myvalue);
            int H = myvalueInt/255;
            int S = myvalueInt%255;
            //  S = 126;
            myvalueInt = abs(myvalueInt-255);
            //printf("[%d,%d,%d],(%d,%d,%d) \n",this->value(),this->getHistoryNode()->numberVoxel,this->getMaximumNumberVoxels(),myvalueInt,H,S);
            QColor q = new QColor();
            q.setRgb(myvalueInt,myvalueInt,myvalueInt,255);
            m_backGradient.setColorAt(0, q);
            m_backGradient.setColorAt(1, q);
            painter->setBrush(m_backGradient);
        }else{
            painter->setBrush(m_backGradient);
        }
    }

    QRectF rect = outlineRect();

    //write the box
    if(showBox){
        if (showVoxel){
            QFont font("times", 5);
            QFontMetrics fm(font);
            //QFontMetrics fm = qApp->font();
            rect = fm.boundingRect(QString::number(this->getHistoryNode()->numberVoxel));
            rect.adjust(-padding, -padding, padding, padding);
            rect.translate(-rect.center());
            painter->drawRoundedRect(rect, roundness(rect.width()), roundness(rect.height()), Qt::RelativeSize);
        }else{
            painter->drawRoundedRect(rect, roundness(rect.width()), roundness(rect.height()), Qt::RelativeSize);
        }
    }

    painter->setPen(m_textColor);

    //Write the Labels
    if(showLabel){
        QFont font("times", 5);
        QFontMetrics fm(font);
        painter->setFont(font);
        painter->drawText(rect, Qt::AlignCenter, QString::number(m_value));
    }else{
        if(showVoxel){
            painter->drawText(rect, Qt::AlignCenter, QString::number(this->getHistoryNode()->numberVoxel));
        }
    }
}

QVariant Node::itemChange(GraphicsItemChange change, const QVariant &value) {
    if (change == ItemPositionHasChanged || change == ItemPositionChange) {
        QVector<Edge *> edges = this->edges();
        foreach (Edge *edge, edges)
                edge->trackNodes();
    }
    return QGraphicsItem::itemChange(change, value);
}

void Node::plotNewickEdge() {
     if (this->getFistChild()!=0){
         if (this->getFistChild()==this->getLastChild()){
             if (this->parent()->getFistChild()==this->parent()->getLastChild()){
                this->parentEdge()->trackNodesNewich_VV(this,this->parent());
             }else{
                this->parentEdge()->trackNodesNewich_VN(this,this->parent());
             }
         }else{
            if (this->parent()->getFistChild()==this->parent()->getLastChild()){
               this->parentEdge()->trackNodesNewich_VV(this,this->parent());
            }else{
               this->parentEdge()->trackNodesNewich_VN(this,this->parent());
            }
            this->horizontalEdge()->trackNodesNewich_H(this->getFistChild(), this->getLastChild(),this);
         }
    }else{
         if (this->parent()->getFistChild()==this->parent()->getLastChild()){
            this->parentEdge()->trackNodesNewich_VV(this,this->parent());
         }else{
            this->parentEdge()->trackNodesNewich_VN(this,this->parent());
         }
    }
}

void Node::plotNewickEdgeRoot() {
     if (this->getFistChild()!=0){
         if (this->getFistChild()==this->getLastChild()){
            //this->parentEdge()->trackNodesNewich_VV(this,this);
         }else{
            this->horizontalEdge()->trackNodesNewich_H(this->getFistChild(), this->getLastChild(),this);
         }
    }
}

QSizeF Node::maxSize() {
    QFontMetrics fm = qApp->font();
    //QFont font("times", 4);
    //QFontMetrics fm(font);
    QRectF rect = fm.boundingRect(QString::number(-2147483647));
    rect.adjust(-padding, -padding, padding, padding);
    rect.adjust(-margin, -margin, margin, margin);

    return rect.size();
}

qreal Node::levelSpacing() {
    const int num = 2;
    return maxSize().height() * num;
}

qreal Node::normalHeight() {
    return maxSize().height() - 2 * margin;
}

void Node::hoverEnterEvent(QGraphicsSceneHoverEvent *event) {
    scale(hoverScaleFactor, hoverScaleFactor);
    voxelNumberMessage->setText(QString::number(this->getHistoryNode()->numberVoxel));
    voxelNumberMessage->setPen(QPen(QColor(0,0,0)));
    voxelNumberMessage->setPos(10,10); // it's just for this sample
    m_hover = true;
    //QRectF rect = outlineRect();
    //painter->drawRoundedRect(rect, roundness(rect.width()), roundness(rect.height()), Qt::RelativeSize);
    //painter->drawText(rect, Qt::AlignCenter, QString::number(this->getHistoryNode()->numberVoxel));
    QGraphicsItem::hoverEnterEvent(event);
}

void Node::hoverLeaveEvent(QGraphicsSceneHoverEvent *event) {
    scale(1 / hoverScaleFactor, 1 / hoverScaleFactor);
    voxelNumberMessage->setText("");
    voxelNumberMessage->setPen(QPen(QColor(255,255,255)));
    m_hover = false;
    QGraphicsItem::hoverLeaveEvent(event);
}

void Node::updateSubTreeSize() {
    m_subTreeSize = 1;

    if (m_firstchild){
        m_subTreeSize += m_firstchild->subTreeSize();

        Node *aux = m_firstchild->m_next;
        while(aux){
            Node *n = aux;
            m_subTreeSize += n->subTreeSize();
            aux = n->m_next;
        }
    }
    if (parent())
        parent()->updateSubTreeSize();
}

int Node::subTreeSize() const {
    return m_subTreeSize;
}


int Node::subTreeSlots() const {
    return m_subTreeSlots;
}

void Node::setTreeSlots(int c) {
    m_subTreeSlots=c;
}


void Node::setFistChild(Node *node){
    m_firstchild = node;
}

Node* Node::getFistChild(){
    return m_firstchild;
}

void Node::setLastChild(Node *node){
    m_lastchild = node;
}

Node* Node::getLastChild(){
    return m_lastchild;
}

void Node::setNextChild(Node *node){
    m_next = node;
}

Node* Node::getNextChild(){
    return m_next;
}

void Node::setShowBox(bool b){
    showBox =b;
}

void Node::setShowLabel(bool b){
    showLabel=b;
}

void Node::setShowVoxel(bool b){
    showVoxel=b;
}

bool Node::getShowVoxel(){
    return showVoxel;
}

bool Node::getShowLabel(){
    return showLabel;
}

bool Node::getShowBox(){
    return showBox;
}

historynode* Node::getHistoryNode(){
    return histNode;
}

void Node::setMaximumNumberVoxels(int b){
   maximumNumberVoxels =b;
}

int Node::getMaximumNumberVoxels(){
   return maximumNumberVoxels;
}
