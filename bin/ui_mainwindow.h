/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created: Thu Nov 15 19:16:55 2012
**      by: Qt User Interface Compiler version 4.8.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QDockWidget>
#include <QtGui/QGraphicsView>
#include <QtGui/QGroupBox>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QListWidget>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QPushButton>
#include <QtGui/QScrollArea>
#include <QtGui/QSpacerItem>
#include <QtGui/QStatusBar>
#include <QtGui/QTabWidget>
#include <QtGui/QToolBar>
#include <QtGui/QTreeWidget>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QAction *actionOpen;
    QAction *actionZoomPlus;
    QAction *actionZoomLess;
    QAction *actionJPG_Newick;
    QAction *actionQPlot_Tree;
    QAction *actionPS_PDF;
    QAction *actionShowLabel;
    QAction *actionShowBox;
    QAction *actionShowVoxels;
    QAction *actionShowCV;
    QAction *actionSaveVolume;
    QWidget *centralWidget;
    QHBoxLayout *horizontalLayout_2;
    QScrollArea *scrollAreaTree;
    QWidget *scrollAreaWidgetContents;
    QVBoxLayout *verticalLayout_11;
    QLabel *labelImageTree;
    QMenuBar *menuBar;
    QMenu *menuFile;
    QMenu *menuTreeViewOption;
    QMenu *menuPlot_Tree_Option;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;
    QDockWidget *dockWidget;
    QWidget *dockWidgetContents;
    QVBoxLayout *verticalLayout_7;
    QWidget *widget;
    QVBoxLayout *verticalLayout_6;
    QGroupBox *groupBoxVolume;
    QVBoxLayout *verticalLayout_5;
    QListWidget *listWidgetVolume;
    QTabWidget *tabWidgetOption;
    QWidget *tabInfo;
    QWidget *layoutWidget;
    QVBoxLayout *verticalLayout_4;
    QVBoxLayout *verticalLayout_2;
    QHBoxLayout *horizontalLayout_3;
    QLabel *label;
    QSpacerItem *horizontalSpacer;
    QLineEdit *lineEdit;
    QGraphicsView *graphicsView;
    QWidget *tabTree;
    QVBoxLayout *verticalLayout_8;
    QVBoxLayout *verticalLayout_3;
    QGroupBox *groupBox;
    QVBoxLayout *verticalLayout_9;
    QVBoxLayout *verticalLayout;
    QTreeWidget *treeWidgetTree;
    QHBoxLayout *horizontalLayout;
    QLabel *labelBins;
    QLineEdit *lineEditBins;
    QLabel *labelLambda;
    QLineEdit *lineEditLambda;
    QLabel *labelK;
    QLineEdit *lineEditK;
    QPushButton *pushButtonComputeTree;
    QGroupBox *groupBox_2;
    QVBoxLayout *verticalLayout_10;
    QListWidget *listWidgetSegments;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(660, 968);
        actionOpen = new QAction(MainWindow);
        actionOpen->setObjectName(QString::fromUtf8("actionOpen"));
        actionOpen->setCheckable(false);
        actionOpen->setEnabled(true);
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/img/addvalue.png"), QSize(), QIcon::Normal, QIcon::Off);
        icon.addFile(QString::fromUtf8("img/open.ico"), QSize(), QIcon::Normal, QIcon::On);
        actionOpen->setIcon(icon);
        actionOpen->setIconVisibleInMenu(true);
        actionZoomPlus = new QAction(MainWindow);
        actionZoomPlus->setObjectName(QString::fromUtf8("actionZoomPlus"));
        QIcon icon1;
        icon1.addFile(QString::fromUtf8(":/img/zoomin.png"), QSize(), QIcon::Normal, QIcon::On);
        actionZoomPlus->setIcon(icon1);
        actionZoomLess = new QAction(MainWindow);
        actionZoomLess->setObjectName(QString::fromUtf8("actionZoomLess"));
        QIcon icon2;
        icon2.addFile(QString::fromUtf8(":/img/zoomout.png"), QSize(), QIcon::Normal, QIcon::On);
        actionZoomLess->setIcon(icon2);
        actionJPG_Newick = new QAction(MainWindow);
        actionJPG_Newick->setObjectName(QString::fromUtf8("actionJPG_Newick"));
        actionJPG_Newick->setCheckable(true);
        actionJPG_Newick->setChecked(false);
        actionQPlot_Tree = new QAction(MainWindow);
        actionQPlot_Tree->setObjectName(QString::fromUtf8("actionQPlot_Tree"));
        actionQPlot_Tree->setCheckable(true);
        actionPS_PDF = new QAction(MainWindow);
        actionPS_PDF->setObjectName(QString::fromUtf8("actionPS_PDF"));
        actionPS_PDF->setCheckable(true);
        actionShowLabel = new QAction(MainWindow);
        actionShowLabel->setObjectName(QString::fromUtf8("actionShowLabel"));
        actionShowLabel->setCheckable(true);
        actionShowBox = new QAction(MainWindow);
        actionShowBox->setObjectName(QString::fromUtf8("actionShowBox"));
        actionShowBox->setCheckable(true);
        actionShowBox->setChecked(false);
        actionShowVoxels = new QAction(MainWindow);
        actionShowVoxels->setObjectName(QString::fromUtf8("actionShowVoxels"));
        actionShowVoxels->setCheckable(true);
        actionShowCV = new QAction(MainWindow);
        actionShowCV->setObjectName(QString::fromUtf8("actionShowCV"));
        actionShowCV->setCheckable(true);
        actionSaveVolume = new QAction(MainWindow);
        actionSaveVolume->setObjectName(QString::fromUtf8("actionSaveVolume"));
        QIcon icon3;
        icon3.addFile(QString::fromUtf8(":/img/arrow.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionSaveVolume->setIcon(icon3);
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        horizontalLayout_2 = new QHBoxLayout(centralWidget);
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        scrollAreaTree = new QScrollArea(centralWidget);
        scrollAreaTree->setObjectName(QString::fromUtf8("scrollAreaTree"));
        scrollAreaTree->setWidgetResizable(true);
        scrollAreaWidgetContents = new QWidget();
        scrollAreaWidgetContents->setObjectName(QString::fromUtf8("scrollAreaWidgetContents"));
        scrollAreaWidgetContents->setGeometry(QRect(0, 0, 81, 861));
        verticalLayout_11 = new QVBoxLayout(scrollAreaWidgetContents);
        verticalLayout_11->setSpacing(6);
        verticalLayout_11->setContentsMargins(11, 11, 11, 11);
        verticalLayout_11->setObjectName(QString::fromUtf8("verticalLayout_11"));
        labelImageTree = new QLabel(scrollAreaWidgetContents);
        labelImageTree->setObjectName(QString::fromUtf8("labelImageTree"));

        verticalLayout_11->addWidget(labelImageTree);

        scrollAreaTree->setWidget(scrollAreaWidgetContents);

        horizontalLayout_2->addWidget(scrollAreaTree);

        MainWindow->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(MainWindow);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 660, 25));
        menuFile = new QMenu(menuBar);
        menuFile->setObjectName(QString::fromUtf8("menuFile"));
        menuTreeViewOption = new QMenu(menuBar);
        menuTreeViewOption->setObjectName(QString::fromUtf8("menuTreeViewOption"));
        menuPlot_Tree_Option = new QMenu(menuBar);
        menuPlot_Tree_Option->setObjectName(QString::fromUtf8("menuPlot_Tree_Option"));
        MainWindow->setMenuBar(menuBar);
        mainToolBar = new QToolBar(MainWindow);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        mainToolBar->setLayoutDirection(Qt::LeftToRight);
        MainWindow->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(MainWindow);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        MainWindow->setStatusBar(statusBar);
        dockWidget = new QDockWidget(MainWindow);
        dockWidget->setObjectName(QString::fromUtf8("dockWidget"));
        dockWidget->setMinimumSize(QSize(280, 600));
        dockWidgetContents = new QWidget();
        dockWidgetContents->setObjectName(QString::fromUtf8("dockWidgetContents"));
        dockWidgetContents->setMinimumSize(QSize(280, 0));
        verticalLayout_7 = new QVBoxLayout(dockWidgetContents);
        verticalLayout_7->setSpacing(6);
        verticalLayout_7->setContentsMargins(11, 11, 11, 11);
        verticalLayout_7->setObjectName(QString::fromUtf8("verticalLayout_7"));
        widget = new QWidget(dockWidgetContents);
        widget->setObjectName(QString::fromUtf8("widget"));
        verticalLayout_6 = new QVBoxLayout(widget);
        verticalLayout_6->setSpacing(6);
        verticalLayout_6->setContentsMargins(11, 11, 11, 11);
        verticalLayout_6->setObjectName(QString::fromUtf8("verticalLayout_6"));
        groupBoxVolume = new QGroupBox(widget);
        groupBoxVolume->setObjectName(QString::fromUtf8("groupBoxVolume"));
        verticalLayout_5 = new QVBoxLayout(groupBoxVolume);
        verticalLayout_5->setSpacing(6);
        verticalLayout_5->setContentsMargins(11, 11, 11, 11);
        verticalLayout_5->setObjectName(QString::fromUtf8("verticalLayout_5"));
        listWidgetVolume = new QListWidget(groupBoxVolume);
        listWidgetVolume->setObjectName(QString::fromUtf8("listWidgetVolume"));

        verticalLayout_5->addWidget(listWidgetVolume);


        verticalLayout_6->addWidget(groupBoxVolume);

        tabWidgetOption = new QTabWidget(widget);
        tabWidgetOption->setObjectName(QString::fromUtf8("tabWidgetOption"));
        tabInfo = new QWidget();
        tabInfo->setObjectName(QString::fromUtf8("tabInfo"));
        layoutWidget = new QWidget(tabInfo);
        layoutWidget->setObjectName(QString::fromUtf8("layoutWidget"));
        layoutWidget->setGeometry(QRect(0, 10, 261, 231));
        verticalLayout_4 = new QVBoxLayout(layoutWidget);
        verticalLayout_4->setSpacing(6);
        verticalLayout_4->setContentsMargins(11, 11, 11, 11);
        verticalLayout_4->setObjectName(QString::fromUtf8("verticalLayout_4"));
        verticalLayout_4->setContentsMargins(0, 0, 0, 0);
        verticalLayout_2 = new QVBoxLayout();
        verticalLayout_2->setSpacing(6);
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setSpacing(6);
        horizontalLayout_3->setObjectName(QString::fromUtf8("horizontalLayout_3"));
        label = new QLabel(layoutWidget);
        label->setObjectName(QString::fromUtf8("label"));

        horizontalLayout_3->addWidget(label);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_3->addItem(horizontalSpacer);

        lineEdit = new QLineEdit(layoutWidget);
        lineEdit->setObjectName(QString::fromUtf8("lineEdit"));

        horizontalLayout_3->addWidget(lineEdit);


        verticalLayout_2->addLayout(horizontalLayout_3);


        verticalLayout_4->addLayout(verticalLayout_2);

        graphicsView = new QGraphicsView(layoutWidget);
        graphicsView->setObjectName(QString::fromUtf8("graphicsView"));

        verticalLayout_4->addWidget(graphicsView);

        tabWidgetOption->addTab(tabInfo, QString());
        tabTree = new QWidget();
        tabTree->setObjectName(QString::fromUtf8("tabTree"));
        verticalLayout_8 = new QVBoxLayout(tabTree);
        verticalLayout_8->setSpacing(6);
        verticalLayout_8->setContentsMargins(11, 11, 11, 11);
        verticalLayout_8->setObjectName(QString::fromUtf8("verticalLayout_8"));
        verticalLayout_3 = new QVBoxLayout();
        verticalLayout_3->setSpacing(6);
        verticalLayout_3->setObjectName(QString::fromUtf8("verticalLayout_3"));
        groupBox = new QGroupBox(tabTree);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        verticalLayout_9 = new QVBoxLayout(groupBox);
        verticalLayout_9->setSpacing(6);
        verticalLayout_9->setContentsMargins(11, 11, 11, 11);
        verticalLayout_9->setObjectName(QString::fromUtf8("verticalLayout_9"));
        verticalLayout = new QVBoxLayout();
        verticalLayout->setSpacing(6);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        treeWidgetTree = new QTreeWidget(groupBox);
        QTreeWidgetItem *__qtreewidgetitem = new QTreeWidgetItem();
        __qtreewidgetitem->setText(0, QString::fromUtf8("1"));
        treeWidgetTree->setHeaderItem(__qtreewidgetitem);
        treeWidgetTree->setObjectName(QString::fromUtf8("treeWidgetTree"));

        verticalLayout->addWidget(treeWidgetTree);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setSpacing(6);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        labelBins = new QLabel(groupBox);
        labelBins->setObjectName(QString::fromUtf8("labelBins"));

        horizontalLayout->addWidget(labelBins);

        lineEditBins = new QLineEdit(groupBox);
        lineEditBins->setObjectName(QString::fromUtf8("lineEditBins"));
        lineEditBins->setEnabled(true);

        horizontalLayout->addWidget(lineEditBins);

        labelLambda = new QLabel(groupBox);
        labelLambda->setObjectName(QString::fromUtf8("labelLambda"));
        QFont font;
        font.setFamily(QString::fromUtf8("Lucida Grande"));
        font.setPointSize(14);
        labelLambda->setFont(font);
        labelLambda->setTextFormat(Qt::RichText);

        horizontalLayout->addWidget(labelLambda);

        lineEditLambda = new QLineEdit(groupBox);
        lineEditLambda->setObjectName(QString::fromUtf8("lineEditLambda"));
        lineEditLambda->setEnabled(true);
        lineEditLambda->setReadOnly(false);

        horizontalLayout->addWidget(lineEditLambda);

        labelK = new QLabel(groupBox);
        labelK->setObjectName(QString::fromUtf8("labelK"));

        horizontalLayout->addWidget(labelK);

        lineEditK = new QLineEdit(groupBox);
        lineEditK->setObjectName(QString::fromUtf8("lineEditK"));
        lineEditK->setEnabled(true);

        horizontalLayout->addWidget(lineEditK);


        verticalLayout->addLayout(horizontalLayout);

        pushButtonComputeTree = new QPushButton(groupBox);
        pushButtonComputeTree->setObjectName(QString::fromUtf8("pushButtonComputeTree"));

        verticalLayout->addWidget(pushButtonComputeTree);


        verticalLayout_9->addLayout(verticalLayout);


        verticalLayout_3->addWidget(groupBox);

        groupBox_2 = new QGroupBox(tabTree);
        groupBox_2->setObjectName(QString::fromUtf8("groupBox_2"));
        verticalLayout_10 = new QVBoxLayout(groupBox_2);
        verticalLayout_10->setSpacing(6);
        verticalLayout_10->setContentsMargins(11, 11, 11, 11);
        verticalLayout_10->setObjectName(QString::fromUtf8("verticalLayout_10"));
        listWidgetSegments = new QListWidget(groupBox_2);
        listWidgetSegments->setObjectName(QString::fromUtf8("listWidgetSegments"));

        verticalLayout_10->addWidget(listWidgetSegments);


        verticalLayout_3->addWidget(groupBox_2);


        verticalLayout_8->addLayout(verticalLayout_3);

        tabWidgetOption->addTab(tabTree, QString());

        verticalLayout_6->addWidget(tabWidgetOption);


        verticalLayout_7->addWidget(widget);

        dockWidget->setWidget(dockWidgetContents);
        MainWindow->addDockWidget(static_cast<Qt::DockWidgetArea>(1), dockWidget);

        menuBar->addAction(menuFile->menuAction());
        menuBar->addAction(menuTreeViewOption->menuAction());
        menuBar->addAction(menuPlot_Tree_Option->menuAction());
        menuTreeViewOption->addAction(actionJPG_Newick);
        menuTreeViewOption->addAction(actionQPlot_Tree);
        menuTreeViewOption->addAction(actionPS_PDF);
        menuPlot_Tree_Option->addSeparator();
        menuPlot_Tree_Option->addAction(actionShowLabel);
        menuPlot_Tree_Option->addAction(actionShowBox);
        menuPlot_Tree_Option->addAction(actionShowVoxels);
        menuPlot_Tree_Option->addAction(actionShowCV);
        mainToolBar->addAction(actionOpen);
        mainToolBar->addSeparator();
        mainToolBar->addAction(actionZoomPlus);
        mainToolBar->addAction(actionZoomLess);
        mainToolBar->addSeparator();

        retranslateUi(MainWindow);

        tabWidgetOption->setCurrentIndex(1);


        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "MainWindow", 0, QApplication::UnicodeUTF8));
        actionOpen->setText(QApplication::translate("MainWindow", "Open", 0, QApplication::UnicodeUTF8));
        actionZoomPlus->setText(QApplication::translate("MainWindow", "zoomPlus", 0, QApplication::UnicodeUTF8));
        actionZoomLess->setText(QApplication::translate("MainWindow", "zoomLess", 0, QApplication::UnicodeUTF8));
        actionJPG_Newick->setText(QApplication::translate("MainWindow", "JPG Newick ", 0, QApplication::UnicodeUTF8));
        actionQPlot_Tree->setText(QApplication::translate("MainWindow", "QPlot Tree", 0, QApplication::UnicodeUTF8));
        actionPS_PDF->setText(QApplication::translate("MainWindow", "PS / PDF", 0, QApplication::UnicodeUTF8));
        actionShowLabel->setText(QApplication::translate("MainWindow", "Show Label", 0, QApplication::UnicodeUTF8));
        actionShowBox->setText(QApplication::translate("MainWindow", "Show Box", 0, QApplication::UnicodeUTF8));
        actionShowVoxels->setText(QApplication::translate("MainWindow", "Show Voxel", 0, QApplication::UnicodeUTF8));
        actionShowCV->setText(QApplication::translate("MainWindow", "Show CV", 0, QApplication::UnicodeUTF8));
        actionSaveVolume->setText(QApplication::translate("MainWindow", "SaveVolume", 0, QApplication::UnicodeUTF8));
        labelImageTree->setText(QString());
        menuFile->setTitle(QApplication::translate("MainWindow", "File", 0, QApplication::UnicodeUTF8));
        menuTreeViewOption->setTitle(QApplication::translate("MainWindow", "TreeViewOption", 0, QApplication::UnicodeUTF8));
        menuPlot_Tree_Option->setTitle(QApplication::translate("MainWindow", "Plot Tree Option", 0, QApplication::UnicodeUTF8));
        groupBoxVolume->setTitle(QApplication::translate("MainWindow", "Volumes", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("MainWindow", "Volume Size", 0, QApplication::UnicodeUTF8));
        tabWidgetOption->setTabText(tabWidgetOption->indexOf(tabInfo), QApplication::translate("MainWindow", "Info", 0, QApplication::UnicodeUTF8));
        groupBox->setTitle(QApplication::translate("MainWindow", "Componet Trees", 0, QApplication::UnicodeUTF8));
        labelBins->setText(QApplication::translate("MainWindow", "Bins", 0, QApplication::UnicodeUTF8));
        lineEditBins->setText(QApplication::translate("MainWindow", "0", 0, QApplication::UnicodeUTF8));
        labelLambda->setText(QApplication::translate("MainWindow", "&lambda", 0, QApplication::UnicodeUTF8));
        lineEditLambda->setInputMask(QString());
        lineEditLambda->setText(QApplication::translate("MainWindow", "0", 0, QApplication::UnicodeUTF8));
        labelK->setText(QApplication::translate("MainWindow", "k", 0, QApplication::UnicodeUTF8));
        lineEditK->setText(QApplication::translate("MainWindow", "0", 0, QApplication::UnicodeUTF8));
        pushButtonComputeTree->setText(QApplication::translate("MainWindow", "Run", 0, QApplication::UnicodeUTF8));
        groupBox_2->setTitle(QApplication::translate("MainWindow", "Segmented Volumes", 0, QApplication::UnicodeUTF8));
        tabWidgetOption->setTabText(tabWidgetOption->indexOf(tabTree), QApplication::translate("MainWindow", "Tree", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
