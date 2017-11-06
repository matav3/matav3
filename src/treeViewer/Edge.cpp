#include <QtGui>
#include "Edge.h"
#include "Node.h"
#include <QDebug>

Edge::Edge(Node *fromNode, Node *toNode) {
    m_fromNode = fromNode;
    m_toNode = toNode;

    setZValue(-1);
    setColor(QColor(Qt::black));
    trackNodes();
}

void Edge::setColor(const QColor &color){
    setPen(QPen(color, 2));
}

QColor Edge::color() const {
    return pen().color();
}

Node* Edge::fromNode() const {
    return m_fromNode;
}

void Edge::setFromNode(Node *node) {
    m_fromNode = node;
    setVisible(m_fromNode ? true : false);
}

Node* Edge::toNode() const {
    return m_toNode;
}

void Edge::trackNodes() {
    if (m_fromNode && m_toNode){
        //setLine(QLineF(m_fromNode->pos(), m_toNode->pos()));    //mapFromScene()
    }
}

void Edge::trackNodesNewich_H(Node *f, Node *t,Node *p) {
    //if(bv::debug) printf("[NO_H (%d) (%f,%f)(%f,%f)] ",f->value(),f->pos().x(),p->pos().y(),t->pos().x(),p->pos().y());
    setLine(QLineF(QPointF(f->pos().x(),p->pos().y()), QPointF(t->pos().x(),p->pos().y())));
}

void Edge::trackNodesNewich_VN(Node *f,Node *p) {
    //if(bv::debug) printf("[NO_V (%d) (%f,%f)(%f,%f)] ",f->value(),f->pos().x(),p->pos().y(),p->pos().x(),p->pos().y());
    setLine(QLineF(QPointF(f->pos().x(),f->pos().y()), QPointF(f->pos().x(),p->pos().y())));
}

void Edge::trackNodesNewich_VV(Node *f,Node *p) {
    //if(bv::debug) printf("[NO_V (%d) (%f,%f)(%f,%f)] ",f->value(),f->pos().x(),f->pos().y(),p->pos().x(),p->pos().y());
    setLine(QLineF(QPointF(f->pos().x(),f->pos().y()), QPointF(p->pos().x(),p->pos().y())));
}
