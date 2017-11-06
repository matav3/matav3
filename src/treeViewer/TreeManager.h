#ifndef TREEMANAGER_H_
#define TREEMANAGER_H_
#include <QObject>
#include "../CPPHeaders.hpp"
#include "../constants.hpp"
#include "../deftypes.hpp"
#include <cstdio>
#include "../edgytree.hpp"
#include "../historytree.hpp"
#include "../historytree_bgnd.hpp"
#include "../perbin.hpp"

class TreeScene;
class Node;

class TreeManager : public QObject{
    Q_OBJECT;
public:
    TreeManager();

    bool parseTree(historytree *myTree, int lambda);
    bool removeValue(int value);
    bool deleteSelected();
    void rotateSelectedLeft();
    void rotateSelectedRigth();
    int valueByIndex(int index);
    int subTreeSlots(Node *node);
    int subTreeSlotsDFS(Node *cur);
    int treeSize();
    bool isEmpty();
    TreeScene* scene();
    void printTree();
    void setShowBox(bool b);
    void setShowLabel(bool b);
    void setShowVoxel(bool b);
    bool getShowVoxel();
    bool getShowLabel();
    bool getShowBox();
    void setShowCV(bool b);
    bool getShowCV();
    historynode* TreeManager::getHistoryNodeSelected();
    bool TreeManager::plotDescendant(historynode *myNod, historytree *myTree, Node *nod);
    Node* TreeManager::getNodeSelected();
private slots:
    void updateScene();
    void updateSceneRect();

signals:
    void treeChanged();

private:
    enum ChildType { LeftChild, RightChild };

    void deleteNode(Node *node);
    Node* searchValue(int value);
    int subTreeHeight(Node *node);

    Node* searchDFS(Node *cur, int val);
    int subTreeHeightDFS(Node *cur);
    void subTreeHeightDFS_FLAT(Node *cur);
    void recalcRelPos_FLAT();
    void recalcRelPosDFS_FLAT(Node *node);
    void TreeManager::recalcAbsPosDFS_FLAT(Node *node);
    void TreeManager::recalcAbsPos_FLAT();
    void recalcRelPos();
    void recalcAbsPos();
    qreal recalcRelPosDFS(Node *node, ChildType type);
    void recalcAbsPosDFS(Node *node);
    bool showBox;
    bool showLabel;
    bool showVoxel;
    bool showCV;
    int valueByIndexDFS(Node *node, int index);

    void setRoot(Node *node);

    TreeScene *m_scene;
    Node *m_root;
    int m_treeSize;
};

#endif /* TREEMANAGER_H */
