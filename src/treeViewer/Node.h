#ifndef NODE_H_
#define NODE_H_

#include <QGraphicsItem>
#include <QVector>
#include "../CPPHeaders.hpp"
#include "../constants.hpp"
#include "../deftypes.hpp"
#include <cstdio>
#include "../edgytree.hpp"
#include "../historytree.hpp"
#include "../historytree_bgnd.hpp"
#include "../perbin.hpp"

class Edge;

class Node : public QGraphicsItem {
public:
    //Node(int value);
    Node(historynode* node);

    ~Node();

    int value() const;
    QColor backColor() const;
    void setBackColor(const QColor &color);

    QVector<Edge *> edges();

    QRectF boundingRect() const;
    QPainterPath shape() const;
    void paint(QPainter *painter,
               const QStyleOptionGraphicsItem *option, QWidget *widget);

    static QSizeF maxSize();
    static qreal levelSpacing();
    static qreal normalHeight();

    QPointF relPos() const;
    void setRelPos(const QPointF &relPos);

    Edge* parentEdge() const;
    Edge* horizontalEdge() const;
    Node* parent() const;

    void setFistChild(Node *node);
    Node* getFistChild();
    void setLastChild(Node *node);
    Node* getLastChild();
    void setNextChild(Node *node);
    Node* getNextChild();

    void setShowBox(bool b);
    void setShowLabel(bool b);
    void setShowVoxel(bool b);
    bool getShowVoxel();
    bool getShowLabel();
    bool getShowBox();

    void setMaximumNumberVoxels(int b);
    int getMaximumNumberVoxels();

    void addChild(Node* chp);
    void removeParentEdge();
    void setMeLeaf();

    void updateSubTreeSize();
    int subTreeSize() const;
    int subTreeSlots() const;
    void setTreeSlots(int c);
    historynode* getHistoryNode();
    void plotNewickEdge();
    void plotNewickEdgeRoot();
    void setParent(Node *node);

protected:
    QVariant itemChange(GraphicsItemChange change,
                        const QVariant &value);

    void hoverEnterEvent(QGraphicsSceneHoverEvent *event);
    void hoverLeaveEvent(QGraphicsSceneHoverEvent *event);

private:
    QRectF outlineRect() const;
    int roundness(double size) const;


    void removeParent();
    void removeChild(Node *node);

    static const qreal hoverScaleFactor;
    static const int padding;
    static const int glowWidth;
    static const int margin;

    int m_value;
    int m_subTreeSize;
    int m_subTreeSlots;
    QPointF m_relPos;
    QColor m_backColor;
    QColor m_textColor;
    QColor m_outlineColor;
    QLinearGradient m_backGradient;
    QGraphicsSimpleTextItem* voxelNumberMessage;
    //Node *m_leftChild, *m_rightChild, *f_child, *l_child, *_n_clild;
    Node *m_next, *m_firstchild, *m_lastchild;
    Edge *m_parentEdge;
    Edge *m_horizontalEdge;

    historynode* histNode;
    bool m_hover;
    bool showBox;
    bool showLabel;
    bool showVoxel;
    int maximumNumberVoxels;
};

#endif /* NODE_H_ */
