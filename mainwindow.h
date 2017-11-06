#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QString>
#include "ui_mainwindow.h"
#include "src/CPPHeaders.hpp"
#include "src/constants.hpp"
#include "src/deftypes.hpp"
#include <cstdio>
#include "src/DIGFile.h"
#include "src/matrix3d.hpp"
#include "src/edgytree.hpp"
#include "src/historytree.hpp"
#include "src/historytree_bgnd.hpp"
#include "src/perbin.hpp"
#include "src/dsf_for_tree_gen.hpp"
#include "src/dsf_for_tree_bgnd_gen.hpp"
#include "src/image3d.hpp"
#include "src/treeViewer/TreeManager.h"
#include "src/treeViewer/TreeView.h"
#include "src/treeViewer/Node.h"
#include "src/treeViewer/TreeManager.h"
#include "src/treeViewer/TreeScene.h"

using namespace std;
namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    void parameterLine(QString name);
    ~MainWindow();

private:
    Ui::MainWindow *ui;

    //list of objects loaded
    QMap<QString,QString> mapVolumeName; // <base_name, complete_name>
    allmydesciptors_t allmyDescpt; //std::map<QString, allmybins_t*> structuresMap (the key is the complete_name
    QVector<bool> treeVolumeVector;
    QVector<int> maxBinsInVolumeVector;
    //actual objects
    image3d<float>* acImage3dFloaP;
    image3d<double>* acImage3dDoubleP;

    //validAncestors is an array that will be used to hold pointers to vertex
    //objects, for use in the procedure computeOutTree(). It is assumed to be
    //an array of size at least K, where K is the maximum over all leaves x of
    //T, of the number of critical ancestors of x in T
    historynode** validAncestors; // I'm not sure why I insert this pointer here. Verify
    historytree* fcts_tree; // I'm not sure why I insert this pointer here. Verify
    QString acVolumeName;  //Name for the actual volume used in the matav3
    allmybins_t acAllmybeansP; //Allmybeans for the actual volume used in the matav3
    int acBeans;
    params_t acParams; //Parametrs for the actual volume used in the matav3
    DIGFile digvorfil;
    bool openNewVolume;
    double scaleFactor;
    void putDIGFileIntoI3D(QString name);
    template <class codomain> void showDescForCodomain(int max_numbins,int numbins_i, int minsize_i,  int lambda_i, image3d<codomain>& myI3D);
    void plotTreeUsingPhylip(QString volumeName, int numbins, params_t params,int X, int Y);
    void showTreeJpg(QString volumeName, int bin, int k, int lambda,int X, int Y);
    QString MainWindow::getJpgTree(int X, int Y);
    void createQplotTree(QString volumeName, int bin, int k, int lambda);


    TreeManager* m_tree;
    TreeView *m_view;
    int treeViewOption; // 0-jpg, 1-qtPlot, 2-Pdf/Ps
    bool showBoxC;
    bool showLabelC;
    bool showVoxelC;
    bool showCV;
    bool showAllVoxels;
    QString msgTime;

private slots:
    void numberChanged();
    bool openDIGFile();
    void updateVolumeList();
    void updateTreeList(int bin, params_t par);
    void computeTree();
    void displayTree();
    void zoomPlus();
    void zoomLess();
    void saveVolume();
    void treeViewOptionJPG();
    void treeViewOptionQPlot();
    void treeViewOptionPSPDF();
    bool verifyTocomputeTree();
    bool verifyToDisplayTree();
    void showQTreeLabel();
    void showQTreeBox();
    void showQTreeVoxels();
    void showQTreeCV();
    void showDescendant();
    void computeOcuppance();


};

#endif // MAINWINDOW_H
