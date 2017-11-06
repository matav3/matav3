#ifndef EDGE_H_
#define EDGE_H_

#include <QGraphicsLineItem>

class Node;

class Edge : public QGraphicsLineItem {
public:
    Edge(Node *fromNode, Node *toNode);

    void setColor(const QColor &color);
    QColor color() const;

    Node* fromNode() const;
    void setFromNode(Node *node);
    Node* toNode() const;

    void trackNodesNewich_H(Node *f, Node *t,Node *p);
    void trackNodesNewich_VN(Node *t,Node *p);
    void trackNodesNewich_VV(Node *t,Node *p);
    void trackNodes();
private:
    Node *m_fromNode, *m_toNode;
    qreal thickness;
};

#endif /* EDGE_H_ */
