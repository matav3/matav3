#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QMessageBox>
#include <QFileDialog>
#include <QStatusBar>
#include <QTreeWidget>
#include <QTreeWidgetItem>
#include <QProcess>
#include <QTime>

using namespace std;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    connect(ui->lineEdit,SIGNAL(textChanged(QString)),this,SLOT(numberChanged()));
    connect(ui->actionOpen,SIGNAL(triggered()),this,SLOT(openDIGFile()));
    connect(ui->actionZoomPlus,SIGNAL(triggered()),this,SLOT(zoomPlus()));
    connect(ui->actionZoomLess,SIGNAL(triggered()),this,SLOT(zoomLess()));
    connect(ui->actionSaveVolume,SIGNAL(triggered()),this,SLOT(saveVolume()));
    connect(ui->listWidgetVolume,SIGNAL(itemSelectionChanged()),this,SLOT(updateVolumeList()));
    connect(ui->pushButtonComputeTree,SIGNAL(clicked()),this,SLOT(computeTree()));
    connect(ui->treeWidgetTree,SIGNAL(itemSelectionChanged()),this,SLOT(displayTree()));  
    connect(ui->actionPS_PDF,SIGNAL(triggered()),this,SLOT(treeViewOptionPSPDF()));
    connect(ui->actionQPlot_Tree,SIGNAL(triggered()),this,SLOT(treeViewOptionQPlot()));
    connect(ui->actionJPG_Newick,SIGNAL(triggered()),this,SLOT(treeViewOptionJPG()));
    connect(ui->actionShowLabel,SIGNAL(triggered()),this,SLOT(showQTreeLabel()));
    connect(ui->actionShowBox,SIGNAL(triggered()),this,SLOT(showQTreeBox()));
    connect(ui->actionShowVoxels,SIGNAL(triggered()),this,SLOT(showQTreeVoxels()));
    connect(ui->actionShowCV,SIGNAL(triggered()),this,SLOT(showQTreeCV()));
    connect(ui->actionShowDescendant,SIGNAL(triggered()),this,SLOT(showDescendant()));
    connect(ui->actionComputeOcuppance,SIGNAL(triggered()),this,SLOT(computeOcuppance()));


    acImage3dFloaP=0;
    acImage3dDoubleP=0;
    allmyDescpt.structuresMap.clear();
    acAllmybeansP.binsMap.clear();
    DIGFile digvorfil;
    openNewVolume=true;
    ui->treeWidgetTree->setColumnCount(4);
    ui->treeWidgetTree->setSelectionMode(QAbstractItemView::SingleSelection);
    QTreeWidgetItem *compTree = new QTreeWidgetItem(ui->treeWidgetTree);
    QStringList headers;
    headers << tr("Name") << tr("Bins") << tr("L") << tr("k") ;
    ui->treeWidgetTree->header()->resizeSection(0, 200);
    ui->treeWidgetTree->header()->resizeSection(1, 50);
    ui->treeWidgetTree->header()->resizeSection(2, 50);
    ui->treeWidgetTree->header()->resizeSection(3, 50);
    ui->treeWidgetTree->setHeaderLabels(headers);
    scaleFactor =1.0;
    ui->labelImageTree->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding );
    ui->scrollAreaTree->setWidgetResizable(true);
    ui->actionZoomLess->setDisabled(true);
    ui->actionZoomPlus->setDisabled(true);
    m_tree = new TreeManager();
    m_view = new TreeView();
    treeViewOption =0;
    showBoxC=1;
    showLabelC=1;
    showVoxelC=0;
    showCV=1;
    m_tree->setShowCV(showCV);
    showAllVoxels=0;
    QActionGroup *plotTreeGroup = new QActionGroup(this);
    plotTreeGroup->addAction(ui->actionQPlot_Tree);
    plotTreeGroup->addAction(ui->actionPS_PDF);
    plotTreeGroup->addAction(ui->actionJPG_Newick);
    msgTime="Information Time:\n";
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::numberChanged() {
    int numOne=ui->lineEdit->text().toInt();
    ui->lineEditLambda->setText(QString::number(numOne));
}

/** The following 5 function control de plot tree option. **/
void MainWindow::treeViewOptionJPG(){
    if(bv::debug) std::cout << "**************************************treeViewOptionJPG()"<<"\n";
    if (!verifyToDisplayTree()) return;
    ui->actionJPG_Newick->setChecked(true);
    ui->actionShowBox->setChecked(false);
    ui->actionShowLabel->setChecked(false);
    treeViewOption =0;
    displayTree();
}

void MainWindow::treeViewOptionQPlot(){
    if(bv::debug) std::cout << "**************************************treeViewOptionQPlot()"<<"\n";
    if (!verifyToDisplayTree()) return;
    ui->actionQPlot_Tree->setChecked(true);
    treeViewOption =1;
    m_tree->setShowLabel(showLabelC);
    m_tree->setShowBox(showBoxC);
    ui->actionShowBox->setChecked(true);
    ui->actionShowLabel->setChecked(true);
    ui->actionShowCV->setChecked(true);
    displayTree();
}

void MainWindow::treeViewOptionPSPDF(){
    if(bv::debug) std::cout << "**************************************treeViewOptionPSPDF()"<<"\n";
    if (!verifyToDisplayTree()) return;
    treeViewOption =2;
    ui->actionPS_PDF->setChecked(true);
    QMessageBox::warning(this,"MATAV3","Ops! No code yet!",
                   "I know. I will do soon...");
}

void MainWindow::showQTreeLabel(){
    if(bv::debug) std::cout << "**************************************showQTreeLabel()"<<"\n";
    if (!verifyToDisplayTree()) return;
    if(showLabelC){
       showLabelC = false;
       m_tree->setShowLabel(showLabelC);
    }else{
       showLabelC = true;
       m_tree->setShowLabel(showLabelC);
       m_tree->setShowVoxel(!showLabelC);
    }
    treeViewOption =1;
    ui->actionQPlot_Tree->setChecked(true);
    displayTree();
}

void MainWindow::showQTreeBox(){
    if(bv::debug) std::cout << "**************************************showQTreeBox()"<<"\n";
    if (!verifyToDisplayTree()) return;
    if(showBoxC){
       showBoxC = false;
       m_tree->setShowBox(showBoxC);
    }else{
       showBoxC = true;
       m_tree->setShowBox(showBoxC);
    }
    ui->actionQPlot_Tree->setChecked(true);
    treeViewOption =1;
    displayTree();
}

void MainWindow::showQTreeVoxels(){
    if(bv::debug) std::cout << "**************************************showQTreeVoxels()"<<"\n";
    if (!verifyToDisplayTree()) return;
    if(showVoxelC){
       showVoxelC = false;
       m_tree->setShowVoxel(showVoxelC);
       ui->actionShowVoxels->setChecked(showVoxelC);
    }else{
       showVoxelC = true;
       m_tree->setShowVoxel(showVoxelC);
       ui->actionShowVoxels->setChecked(showVoxelC);
    }
    treeViewOption =1;
    displayTree();
}

void MainWindow::showQTreeCV(){
    if(bv::debug) std::cout << "**************************************showQTreeVoxels()"<<"\n";
    if (!verifyToDisplayTree()) return;
    if(showCV){
       showCV = false;
       m_tree->setShowCV(showCV);
       ui->actionShowCV->setChecked(showCV);
    }else{
       showCV = true;
       m_tree->setShowCV(showCV);
       ui->actionShowCV->setChecked(showCV);

    }
    treeViewOption =1;
    displayTree();
}
/**
This function get the Tree selected and retrivel the History tree. Then
we simplify this tree and create a TreeView with the History tree nodes.
The treeView has the coordenates information for the voluem
**/
void MainWindow::createQplotTree(QString volumeName, int bin, int k, int lambda){
    QTime myTimer;
    myTimer.start();

    if(bv::debug) std::cout << "createQplotTree ="<<acVolumeName.toAscii().constData()<<"\n";
    int num_bins = bin;
    if (ui->lineEditBins->text().toInt()==0){
       QListWidgetItem *item =ui->listWidgetVolume->currentItem();
       int i = ui->listWidgetVolume->row(item);
       num_bins= maxBinsInVolumeVector.at(i);
    }

    allmybins_t *acAllmybeansP = allmyDescpt[acVolumeName];
    /** There is no allmyDescpt for this volume**/
    if (!acAllmybeansP){
        std::cout << "The FHTS for selected volume is computed for the first time" << std::endl;
        return;
    }

    std::cout << "The FHTS for selected volume is ALREADY computed. num_bins=" << num_bins<< std::endl;

    perbinsize_t* curBinP = (*acAllmybeansP)[num_bins]; // don't leave curBinP as null !

    if (!curBinP){
        std::cout << "The perbinsize_t for selected allmybins_t is NOT computed yet" << std::endl;
        return;
    }

    std::cout << "going to be parsed  [m_tree->treeSize()=" << m_tree->treeSize() << std::endl;



    if (lambda!=0 || k!=0){
        historytree * http = acAllmybeansP->get_historyTree(bin,k,lambda);
        m_tree->parseTree(http,lambda);
    }else{
        m_tree->parseTree(curBinP->htfgp,lambda);
    }

    std::cout << "createQplotTree-[" << myTimer.elapsed() <<" ms]" <<std::endl;
}

/** This function will receive the args from command line and pass to the mainwindown**/

void MainWindow::parameterLine(QString name){
    openNewVolume=false;
    acVolumeName = name;
    openDIGFile();
}


bool MainWindow::openDIGFile() {
    QTime myTimer;
    myTimer.start();

    ui->statusBar->showMessage("Open DIGFile...",1000);
    if(bv::debug) std::cout << "openDIGFile ="<<acVolumeName.toAscii().constData()<<"\n";
    QString digvorfilname;
    digvorfil.Close();
    if (openNewVolume){
        if(bv::debug) std::cout << "open NEW \n";
        digvorfilname = QFileDialog::getOpenFileName(this, tr("Open DIG_Voronoi File")/*currentMovieDirectory*/);
    }else{
        digvorfilname = acVolumeName;
        if(bv::debug) std::cout << "open OLD ="<<digvorfilname.toAscii().constData()<<"\n";
        openNewVolume=true;
    }

    QMap<QString,QString>::const_iterator it =  mapVolumeName.find(digvorfilname);
    if (it != mapVolumeName.end()){
        QMessageBox::warning(this,"MATAV3","The DIGFile "+digvorfilname+" is already opened",
                       "I know. Fail in open");
    }

    if (digvorfilname.isEmpty()){
        QMessageBox::warning(this,"MATAV3","Selection cancelled",
                       "I know. Fail in open");
        return false;
    }else{
        if(0 != digvorfil.Open(digvorfilname.toAscii().constData())) {
            digvorfil.Close(); // this enables us to read another file successfully
            QString err = "Erro to open `" + digvorfilname + "' \ sure that it is a readable AND writable valid DIGFile.";
            QMessageBox::information(this,
                     "openDIGVoronoi error",
                     err,
                     "I shall choose wisely in the future.");
            return false;
        }else{//success...

            // Get information about the file names
            QStringList list = digvorfilname.split("/");
            QString outdir="";
            for (int var=0 ; var < (list.size()-1); var++){
                outdir.append(list.at(var));
                outdir.append("/");
            }
            QString digvorfilname_nopath = list.at(list.size()-1);// get part after last '/'

            QStringList nameComplete = digvorfilname_nopath.split(".");
            QString base = nameComplete.at(0);

            // Creat a Dir to put the save the images
            QString outputbasedir_str = outdir+"/"+"output";
            QDir outputbasedir(outputbasedir_str);
            if(! outputbasedir.exists ()) {
                QDir d;
                if ( d.mkdir(outputbasedir_str) ) { // this feels weird
                    std::cerr << "Folder 'output SUCCESSFULLY created." << outputbasedir_str.toAscii().constData() << std::endl;
                } else {
                    std::cerr << "FAILED to create the folder 'output'." << outputbasedir_str.toAscii().constData()<< std::endl;
                }
            } else {
                std::cerr << "The folder 'output' alread exist." << outputbasedir_str.toAscii().constData() << std::endl;
            }

            // if directory ./output/digvornoext does not exist, create
            outputbasedir_str = outputbasedir_str + "/" + base;
            QDir outputbasedir2(outputbasedir_str);

            if(! outputbasedir2.exists() ) {
                QDir d;
                if (d.mkdir(outputbasedir_str)) {
                    std::cerr << "Folder 'output SUCCESSFULLY created." << outputbasedir_str.toAscii().constData() << std::endl;
                } else {
                    std::cerr << "FAILED to create the folder 'output'." << outputbasedir_str.toAscii().constData() << std::endl;
                }
            } else {
                std::cerr << "The folder 'output' alread exist." << outputbasedir_str.toAscii().constData() << std::endl;
            }

            if(bv::debug)
                std::cout << "Loaded File" << std::endl;

            if (openNewVolume){
                mapVolumeName.insert(base,digvorfilname);
                QListWidgetItem *item = new QListWidgetItem(base);
                ui->listWidgetVolume->addItem(item);
                ui->listWidgetVolume->setCurrentItem(item);
                ui->statusBar->showMessage("DIGFlipe "+digvorfilname_nopath+ " successfully opened!",5000);
            }
            return true;

        }
    }
    std::cout << "openDIGFile-[" << myTimer.elapsed() <<" ms]" <<std::endl;
}


void MainWindow::computeTree(){
    QTime myTimer;
    myTimer.start();
    msgTime= msgTime+ "### computeTree ### \n";
    if(bv::debug) std::cout << "computeTree()"<<"\n";
    treeViewOption =0;
    ui->actionJPG_Newick->setChecked(true);
    if (verifyTocomputeTree()){
        if(bv::debug) std::cout << "acVolumeName="<<acVolumeName.toAscii().constData()<<"\n";
        openNewVolume=false;
        if (!openDIGFile()){
            return 0;
        }

        ui->statusBar->showMessage("Computing a component tree for: "+acVolumeName+ " !",5000);
        putDIGFileIntoI3D(acVolumeName);
        msgTime= msgTime+"      putDIGFileIntoI3D: [" + QString::number(myTimer.elapsed()) +" ms]\n";
        cout << "    putDIGFileIntoI3D: [" << myTimer.elapsed()  <<"ms]"<< std::endl;

        int tim = myTimer.elapsed();
        int numbins_i,max_numbins;
        int k = ui->lineEditK->text().toInt();
        int lambda = ui->lineEditLambda->text().toInt();
        openNewVolume=true;
        /** No tree for this volume **/
        if (treeVolumeVector.at(ui->listWidgetVolume->currentRow())){
            allmybins_t * mb= new allmybins_t;
        }


        cout << "\n\n ********  ATENCAO!!  THE MAXIMUM OF BEANS IS 1000000 ******** \n\n";
        max_numbins = 1000000;
        /*if (acImage3dFloaP){
            max_numbins = 	acImage3dFloaP->calculeFloatBins();
        }else{
            max_numbins = 	acImage3dDoubleP->calculeDoubleBins();
        }*/
         msgTime= msgTime+"      calculeFloatBins: [" + QString::number(myTimer.elapsed()-tim) +" ms]\n";
         msgTime= msgTime+"      max_numbins= " + QString::number(max_numbins) +"\n";
        cout << " .....calculeFloatBins: [" << myTimer.elapsed() <<" ms]" << std::endl;
        if(ui->lineEditBins->text().toInt()==0){
            numbins_i = max_numbins;
        }else{
            numbins_i = ui->lineEditBins->text().toInt();
        }
        if(acImage3dFloaP) {
            showDescForCodomain<float>(max_numbins, numbins_i, k, lambda, (*acImage3dFloaP));
        } else{
            if(acImage3dDoubleP) {
                showDescForCodomain<double>(max_numbins, numbins_i, k, lambda, (*acImage3dDoubleP));
            }else{
                QMessageBox::critical(this,"MATAV3","No image3d in in MATAV3::computeTree()!!","ok");
            }
        }

    }
    msgTime= msgTime+ "Final: [ " + QString::number(myTimer.elapsed()) +" ms] \n";

    cout << msgTime.toAscii().constData() << std::endl;

}
/******************************************************************************
* Function dynamically compute desired descriptor as required, retain in RAM once computed
*   @param [type param: implementationwise 'codomain' of floating-point image,]
*   number of bins, true if bg, to prune: minimum size (voxels), minimum height,
*		image3d of specified 'codomain' (passed by reference) and stick the cascade of
*		parameters as a "visual pointer" into the UI (QListView) as necessary
*
*struct perbinsize_t {
*  historytree*      htfgp;
*  historytree_bgnd* htbgp;
*  edgies_t fg_edgies;
*  edgies_t bg_edgies;
*
* ps:DENIS FUNCTION
******************************************************************************/

template <class codomain>
void MainWindow::showDescForCodomain(int max_numbins, int numbins_i, int minsize_i,  int lambda_i, image3d<codomain>& myI3D) {
    QTime myTimer;
    myTimer.start();
    int tim=0;
    if(bv::debug) printf( " \n Enter in showDescForCodomain. \n");

    unsigned int numbins_ui = static_cast<unsigned int> (numbins_i); // this kind of thing is in bad taste...
    unsigned int minsize_ui = static_cast<unsigned int> (minsize_i); // this kind of thing is in bad taste...


    allmybins_t *allmybins = allmyDescpt[acVolumeName];

    /** There is no allmyDescpt for this volume**/
    if (!allmybins){
        std::cout << "The FHTS for selected volume is computed for the first time" << std::endl;
        allmyDescpt.makeNewStructures(acVolumeName); //create a new allmybeans in my map descriptor
    }

    std::cout << "The FHTS for selected volume is ALREADY computed" << std::endl;
    allmybins = allmyDescpt[acVolumeName]; //get the allmybeans in my map descriptor
    if (!allmybins){
        allmybins->makeNewBin(numbins_i); // create a perbinsize_t for bin number
    }

    perbinsize_t* curBinP = (*allmybins)[numbins_i]; // don't leave curBinP as null !

    /** There is no allMyBeans for this volume. Then, I will create one and set binnedImgP for the specific bin**/
    if (!curBinP){
        std::cout << "The FHTS for selected number of bins (" <<
                            numbins_i << ") not calculated before in this run" << std::endl;


        (*allmybins).makeNewBin(numbins_i);
        curBinP = (*allmybins)[numbins_i];
        curBinP->htfgp = 0;
        maxBinsInVolumeVector.push_back(max_numbins);
        image3d<long long>* binnedImgP;
        if (ui->lineEditBins->text().toInt()==0){
            printf( "Call myI3D.codomain_to_long()  \n");
            binnedImgP = myI3D.codomain_to_long(); // new image3d
        }else{
            printf( "Call myI3D.codomain_to_long_Nbins()  \n");
            binnedImgP = myI3D.codomain_to_long_Nbins(numbins_i);
        }
        msgTime= msgTime+"      codomain_to_long_[" + QString::number(myTimer.elapsed()) +" ms]\n";
        cout <<"............................codomain_to_long: [" << (myTimer.elapsed()) <<" ms]" << std::endl;
        cout <<"............................sizeof(long long)= [" << sizeof(long long) <<" sizeof(float) =" << sizeof(float) <<std::endl;
        myTimer.restart();
        // In this point is construct the history tree
         curBinP->htfgp = binnedImgP->construct_history_tree_long_bins();
         msgTime= msgTime+"      construct_history_tree_long_bins-[" + QString::number(myTimer.elapsed()) +" ms]\n";
         cout <<"............................construct_history_tree_long_bins: [" << myTimer.elapsed() <<" ms]" << std::endl;
         //curBinP->htfgp->updateAllVoxesTree(curBinP->htfgp->hroot); // ATTENTION!! This is temporary. The voxel should not be copied for all simplified tree.
         //exit(0);
         myTimer.restart();
         curBinP->htfgp->updateNumberVoxelTree(curBinP->htfgp->hroot);
         msgTime= msgTime+"      updateNumberVoxelTree-[" + QString::number((myTimer.elapsed())) +" ms]\n";
         cout <<" ............................updateNumberVoxelTree: [" << (myTimer.elapsed()) <<" ms]" << std::endl;

        // Get information about the file names
        QStringList list = acVolumeName.split("/");
        QString outdir="";
        for (int var=0 ; var < (list.size()-1); var++){
            outdir.append(list.at(var));
            outdir.append("/");
        }
        printf("(GIL) a=%lld, b=%lld, c=%lld, d=%lld \n",myI3D.minOfImage,myI3D.maxOfImage,curBinP->htfgp->myNegInf,curBinP->htfgp->myPosInf );
        printf("(GIG - e=%lld f=%lld  \n",binnedImgP->maxOfImage,binnedImgP->minOfImage );
        curBinP->htfgp->myPosInf= binnedImgP->maxOfImage;
        curBinP->htfgp->myNegInf = binnedImgP->minOfImage;


        printf("g=%lld,h=%lld, i=%lld, j=%lld \n",binnedImgP->maxOfImage,binnedImgP->minOfImage,curBinP->htfgp->myNegInf,curBinP->htfgp->myPosInf );

        //Preparando Edgies for the simplification
        params_t zeroparam = params_t(0,0);
        edgytree* zefgp = curBinP->htfgp->generate_pruned_edgy_tree();
        zefgp->flatten();
        zefgp->setparents();
        curBinP->fg_edgies[zeroparam] = zefgp;
        delete binnedImgP; // no use hanging on to this binned image once we have its fg and bg history trees
        plotTreeUsingPhylip(acVolumeName,numbins_i,zeroparam,800,600);        
        //treeViewOption =0;
        //ui->actionJPG_Newick->setChecked(true);
        updateTreeList(numbins_i,zeroparam);
        treeViewOptionJPG();
        //showTreeJpg(acVolumeName,numbins_i,0, 0,800,600);
        ui->actionZoomLess->setDisabled(false);
        ui->actionZoomPlus->setDisabled(false);
        ui->actionJPG_Newick->setChecked(true);
    }else{
         std::cout << "The bg and fg trees for selected number of bins (" <<
            numbins_i << ") WAS BEFORE CALCULATED before in this run" << std::endl;
    }

    tim= myTimer.elapsed();
    /** There is FCTS for this volume and this beans **/
    if( (minsize_i==0) && (lambda_i==0) ){
        return; // we must have been done already
    }

    params_t selectedTuple = params_t(minsize_i,lambda_i);
    historytree* run_htfgp = new historytree();
    //*run_htfgp = *curBinP->htfgp;
    myTimer.start();;;
    curBinP->htfgp->cloneTree(run_htfgp);
    fcts_tree = run_htfgp;
    msgTime= msgTime+"      cloneTree-[" + QString::number((myTimer.elapsed())) +" ms]\n";
    cout <<"............................ cloneTree: [" << (myTimer.elapsed()) <<" ms]" << std::endl;

    //Prune by K
    tim= myTimer.elapsed();
    run_htfgp->prune(minsize_i, 0);
    //run_htfgp->ClearPrunedNode();
    //if (minsize_i!=0) run_htfgp->pruneByDeltaStepB(2); // IWCIA
    run_htfgp->ClearPrunedNode();
    //curBinP->htfgp->prune(minsize_i, 0);
    //curBinP->htfgp->ClearPrunedNode();


    //Apply pruning version 4. Implementation of the new Yung Version dated in 03/10/09
    if( lambda_i!=0 ){        
        historynode* m = run_htfgp->hroot;
        if (run_htfgp->isCriticalVertex(m)){
            run_htfgp->computeIPCD(m,1);
            if(bv::debug) printf("the root is CV (%lld) \n",m->gray);
        }else{
            m=	run_htfgp->immediateCriticalVertice(m);
            run_htfgp->computeIPCD(m,1);
            if(bv::debug) printf("the root is NOT a CV. The m is= (%lld) \n",m->gray);
        }

        /* Compute the farthestDescendent and DeltaLimite*/
        run_htfgp->computeFarthestDescendantsAndDeltaLimits(m);

        run_htfgp->fastPruningFHT(validAncestors,lambda_i,lambda_i);
    }
     msgTime= msgTime+"      pruningTree-[" + QString::number((myTimer.elapsed()-tim))+" ms]\n";
    cout <<" .....pruningTree: [" << (myTimer.elapsed()-tim) <<" ms]" << std::endl;
    tim= myTimer.elapsed();
    //Create the Edge
    edgytree* efgp;
    if ( lambda_i!=0 ){
        efgp =  run_htfgp->generate_pruned_edgy_treeFHCS();
    }else{
        efgp =  run_htfgp->generate_pruned_edgy_tree();
    }
    efgp->setparents();
    msgTime= msgTime+"      generate_pruned_edgy_treeFHCS-[" + QString::number((myTimer.elapsed()-tim)) +" ms]\n";

    edgies_t& edgies = curBinP->fg_edgies;
    edgies[selectedTuple] = efgp;
    histories_t& hors = curBinP->fg_histories; // new! I will store all the simplification for the HistoryTree
    hors[selectedTuple] = run_htfgp;
    //efgp->testando(efgp->eroot);
    plotTreeUsingPhylip(acVolumeName,numbins_i,selectedTuple,800,600);
    showTreeJpg(acVolumeName,numbins_i,minsize_i, lambda_i,800,600);
    ui->actionJPG_Newick->setChecked(true);
    updateTreeList(numbins_i,selectedTuple);
    ui->actionZoomLess->setDisabled(false);
    ui->actionZoomPlus->setDisabled(false);

    if( edgies.find(selectedTuple) != edgies.end()) {
    if(bv::debug)std::cout << " TEM o parameters (bins " <<	numbins_i << ", minsize_i= " << minsize_i<< " lambda="<<  lambda_i << std::endl;
    } else {
    if(bv::debug)std::cout << " NAO TEM o parameters  (bins " <<	numbins_i << ", minsize_i= " << minsize_i<< " lambda="<<  lambda_i << std::endl;
    }
    if(bv::debug) printf(" node CurBinP= %d and run_htfgp = %d  delta1=%lld\n",curBinP->htfgp->numnodes,run_htfgp->numnodes,run_htfgp->delta1);
    //allmybinsP = &allmybins;

} // --botanist::showDescForCodomain<codomain>() definition

/**
This function update the value for the acVolumeName when the user change the selection int the VolumeList.
**/
void MainWindow:: updateVolumeList(){
    if(bv::debug) printf( " \n updateVolumeList \n");
    QListWidgetItem *item =ui->listWidgetVolume->currentItem();
    QMap<QString,QString>::const_iterator it =  mapVolumeName.find(item->text());
    if (it != mapVolumeName.end()){//There is tree for this volume
        acVolumeName= it.value();
        openNewVolume = true;
    }else{
        QString s = item->text();
        int i = ui->listWidgetVolume->row(item);
        QMessageBox::warning(this,"MATAV3","There's NO tree for this volume ->"+i,
                       "ok");
    }    
}

/**
This function create a new item in the TreeList after a new params_t or beans for a
tree be computed. This function is called after  showDescForCodomain()
**/
void MainWindow:: updateTreeList(int bins, params_t param){
    if(bv::debug) printf( " \n updateTreeList \n");
    //get base name for the volume
    QStringList list = acVolumeName.split("/");
    QString base = list.at(list.size()-1).split(".").at(0);

    QList<QTreeWidgetItem*> items = ui->treeWidgetTree->findItems(base,Qt::MatchExactly,0);
    if (items.count()==0){
        QTreeWidgetItem *item = new QTreeWidgetItem(ui->treeWidgetTree->invisibleRootItem());
        item->setText(0,base);
        QTreeWidgetItem *item2 = new QTreeWidgetItem(item);
        item2->setText(0,base);
        item2->setText(1,QString::number(bins));
        item2->setText(2,QString::number(param.delta));
        item2->setText(3,QString::number(param.minsize));
        ui->treeWidgetTree->setCurrentItem(item2);
        //item2->setSelected(true);

    }else{
        QTreeWidgetItem *item = items.first();
        QTreeWidgetItem *item2 = new QTreeWidgetItem(item);
        item2->setText(0,base);
        item2->setText(1,QString::number(bins));
        item2->setText(2,QString::number(param.delta));
        item2->setText(3,QString::number(param.minsize));
        ui->treeWidgetTree->setCurrentItem(item2);
        //item2->setSelected(true);

    }
}

void MainWindow::putDIGFileIntoI3D(QString name) {
     ui->statusBar->showMessage("Creating image3D",1000);
    if(bv::debug) std::cout << "putDIGFileIntoI3D("<<name.toAscii().constData()<<")";
    if(acImage3dFloaP) {
        delete acImage3dFloaP;
        acImage3dFloaP = 0;
    }
    if(acImage3dDoubleP) {
        delete acImage3dDoubleP;
        acImage3dDoubleP = 0;
    }

    //  skip error checking (to be implemented later)
    digvorfil.SelectArraySet(0);
    digvorfil.SelectArray(0);
    DIGDataType DataType; // DIG_Float, DIG_Double, etc.
    digvorfil.GetDataType(&DataType);
    const DIGDimensions* Dimensions;
    digvorfil.GetDimensions(&Dimensions);

    dimSizeType zdim = static_cast<dimSizeType> (Dimensions->z);
    dimSizeType ydim = static_cast<dimSizeType> (Dimensions->y);
    dimSizeType xdim = static_cast<dimSizeType> (Dimensions->x);

    switch(DataType) {
        case DIGDataType_FLOAT:
            std::cout << "I'm made of FLOAT  " << std::endl;
            acImage3dFloaP = new image3d<float>;
            acImage3dFloaP->set_dims(zdim, ydim, xdim);
            acImage3dFloaP->load(digvorfil);
            digvorfil.Close();
            break;
         case DIGDataType_DOUBLE:
            std::cout << "I'm made of DOUBLE" << std::endl;
            acImage3dDoubleP = new image3d<double>;
            acImage3dDoubleP->set_dims(zdim, ydim, xdim);
            acImage3dDoubleP->load(digvorfil);
            digvorfil.Close();
           break;
        default:
           std::cerr << "ERROR: opened DIGVoronoi file has UNSUPPORTED DataType" << std::endl;
           std::cerr << "(must be DIGDataType_FLOAT or DIGDataType_DOUBLE)" << std::endl;
    }
    //listWidget->setSortingEnabled(true);
}

/******************************************************************************
* Function to get a eggy tree for a  FHT and build a tree using open the programa
* 'drawgram' (part of PHYLIP, v3.6 then modified by Deniz) that need be  in the path
******************************************************************************/
void MainWindow::plotTreeUsingPhylip(QString volumeName, int numbins, params_t params,int X, int Y){
    if(bv::debug) printf( " \n plotTreeUsingPhylip \n");
    allmybins_t *allmybinsP = allmyDescpt[volumeName];

    edgytree* etp = allmybinsP->get_edgy(numbins, params.minsize, params.delta);
    QStringList list = volumeName.split("/");
    QString base = list.at(list.size()-1).split(".").at(0);

    QString labImg=base+"_B"+QString::number(numbins,10);
    labImg+="_K"+QString::number(params.minsize,10);
    labImg+="_L"+QString::number(params.delta,10);
    labImg+="_"+QString::number(X,10);
    labImg+="x"+QString::number(Y,10);

    QString outdir="";
    for (int var=0 ; var < (list.size()-1); var++){
        outdir.append(list.at(var));
        outdir.append("/");
    }

    QString NewickFileName = outdir + "output"+"/"+ base +"/"+labImg+ ".new";
    QString XBMFileName = outdir + "output"+"/"+ base +"/"+labImg+ ".xbm";
    QString jpgName = outdir + "output"+"/"+ base +"/"+labImg+".jpg";

    std::cout  << "obtaining new Newick representation from edgytree and writing it to " << NewickFileName.toAscii().constData() << std::endl;
    std::ofstream NewickFile (NewickFileName.toAscii().constData(), std::ios::out); // replace if exists!
    if (NewickFile.is_open()) {
        etp->NEXUS_out(NewickFile); // edgytree::NEXUS_out() writes in the simpler Newick format rather than NEXUS
        NewickFile.close();
    } else {
        cout << "Error opening " << NewickFileName.toAscii().constData() << " for writing" << endl;
        return;
    }



    QString command ="";
    if (false){
        cout << "QProcess comming... \n";
        command = "" +  // botanist warns user as first thing if no drawgram in path
        NewickFileName + " " + // input file
        XBMFileName+ " " +QString::number(X,10)+ " " +QString::number(Y,10);   // output file
        QProcess comQt;
        QStringList clist;
        //comQt.execute("drawgram",clist);
        clist << NewickFileName << XBMFileName << QString::number(X,10) << QString::number(Y,10);
        comQt.start("drawgram",clist);
        if (!comQt.waitForStarted())
            cout << "ERROR 01 \n";

        comQt.write("Qt rocks! \n");
        comQt.closeWriteChannel();

        if (!comQt.waitForFinished())
            cout << "ERROR 02\n ";

        cout << "QProcess BYE... \n";
        //grep.waitForFinished(120000);
    }else{
        command = "drawgram " +  // botanist warns user as first thing if no drawgram in path
        NewickFileName + " " + // input file
        XBMFileName+ " " +QString::number(X,10)+ " " +QString::number(Y,10);   // output file
        system(command.toAscii().constData());

    }



    cout << command.toAscii().constData() <<"\n";
    command ="";
    command = "convert -flip " +  // botanist warns user as first thing if no drawgram in path (-rotate -180)
    XBMFileName + " " + // input file
    jpgName;   // output file
    system(command.toAscii().constData());
    cout << command.toAscii().constData() <<"\n";
    if(bv::debug) printf( " \n Living plotTreeUsingPhylip \n");

} // plotTreeUsingPhylip()

void MainWindow::showTreeJpg(QString volumeName, int bin, int k, int lambda, int X, int Y){
    if(bv::debug) printf( " \n showTreeJpg() \n");
    QStringList list = volumeName.split("/");
    QString base = list.at(list.size()-1).split(".").at(0);
    QString labImg=base+"_B"+QString::number(bin,10);
    labImg+="_K"+QString::number(k,10);
    labImg+="_L"+QString::number(lambda,10);
    labImg+="_"+QString::number(X,10);
    labImg+="x"+QString::number(Y,10);

    QString outdir="";
    for (int var=0 ; var < (list.size()-1); var++){
        outdir.append(list.at(var));
        outdir.append("/");
    }

    QString jpgTree = outdir + "output"+"/"+ base +"/"+labImg+ ".jpg";
    if(bv::debug) cout <<" jpgTree= "<<jpgTree.toAscii().constData() <<"\n";
    QImage image(jpgTree);
    if (image.isNull()) {
        //QMessageBox::information(this, tr("MATAV3"),
        //                         tr("Cannot load %1").arg(jpgTree));
        return;
    }
    if(bv::debug) cout <<"next step set   ui->labelImageTree->setPixmap(QPixmap::fromImage(image));\n";
    ui->labelImageTree->setPixmap(QPixmap::fromImage(image));
    if(bv::debug) cout <<" ui->labelImageTree->adjustSize() \n";
    ui->labelImageTree->adjustSize();
}
/**
This function show the jpg picture for a params_t/bin component tree. This
function is called with the user select a item in the list.
TODO_: show the tree wiht the same zoom factor that was used before. Rigth now starting show
the 800x600
**/

void MainWindow::displayTree(){
    if(bv::debug)  std::cout <<  " \n displayTree() \n";
    QList<QTreeWidgetItem *>  listT = ui->treeWidgetTree->selectedItems();
    QTreeWidgetItem  *item = listT.first();
    QString base = item->text(0);
    int bins = item->text(1).toInt();
    int lambda = item->text(2).toInt();
    int k = item->text(3).toInt();
    if(treeViewOption==0){
        if (item->text(1)!=""){
            if(bv::debug) cout <<" displayTree= "<< bins <<" "<< bins <<" "<< lambda <<" "<< k <<"\n";
            ui->scrollAreaTree->takeWidget();
            ui->scrollAreaTree->setWidget(ui->labelImageTree);
            showTreeJpg(mapVolumeName[base],bins,k,lambda,800,600);
        }
    }else{
        if (treeViewOption==1){
            m_view->setScene(m_tree->scene());
            m_view->setContextMenuPolicy(Qt::ActionsContextMenu);
            ui->scrollAreaTree->takeWidget();
            ui->scrollAreaTree->setWidget(m_view);
            ui->scrollAreaTree->setWindowTitle(tr("Tree Viewer"));
            createQplotTree(mapVolumeName[base],bins,k,lambda);
        }
    }
    updateVolumeList();
}

void  MainWindow::zoomPlus(){
    if (treeViewOption==0){
        if (scaleFactor<4){;
            scaleFactor=scaleFactor + 0.5;
            ui->actionZoomLess->setDisabled(false);
        }else{
            ui->actionZoomPlus->setDisabled(true);
        }
        QString NewickFileName = getJpgTree((int)800*scaleFactor,(int)600*scaleFactor);
    
        QFile Fout(NewickFileName);
        QListWidgetItem *item =ui->listWidgetVolume->currentItem();
        int i = ui->listWidgetVolume->row(item);
        int num_bins = maxBinsInVolumeVector.at(i);
    
        allmybins_t * allmybeans = allmyDescpt[acVolumeName];
        perbinsize_t* curBinP = (*allmybeans)[num_bins];
        int k = ui->lineEditK->text().toInt();
        int lambda = ui->lineEditLambda->text().toInt();
        params_t myparam = params_t(k,lambda);
        if(!Fout.exists()){
            plotTreeUsingPhylip(acVolumeName,num_bins,myparam, (int)800*scaleFactor,(int)600*scaleFactor);
            showTreeJpg(acVolumeName,num_bins,k,lambda,(int)800*scaleFactor,(int)600*scaleFactor);
        }else{
            showTreeJpg(acVolumeName,num_bins,k,lambda,(int)800*scaleFactor,(int)600*scaleFactor);
        }
    }else{
        m_view->zoomIn();
    }

}

void  MainWindow::zoomLess(){
    if (treeViewOption==0){
        if (scaleFactor>1){
            ui->actionZoomPlus->setDisabled(false);
            scaleFactor=scaleFactor - 0.5;
        }else{
            ui->actionZoomLess->setDisabled(true);
        }
        QString NewickFileName = getJpgTree((int)800*scaleFactor,(int)600*scaleFactor);

        QFile Fout(NewickFileName);
        QListWidgetItem *item =ui->listWidgetVolume->currentItem();
        int i = ui->listWidgetVolume->row(item);
        int num_bins = maxBinsInVolumeVector.at(i);

        allmybins_t * allmybeans = allmyDescpt[acVolumeName];
        perbinsize_t* curBinP = (*allmybeans)[num_bins];
        int k = ui->lineEditK->text().toInt();
        int lambda = ui->lineEditLambda->text().toInt();
        params_t myparam = params_t(k,lambda);
        if(!Fout.exists()){
            plotTreeUsingPhylip(acVolumeName,num_bins,myparam, (int)800*scaleFactor,(int)600*scaleFactor);
            showTreeJpg(acVolumeName,num_bins,k,lambda,(int)800*scaleFactor,(int)600*scaleFactor);
        }else{
            showTreeJpg(acVolumeName,num_bins,k,lambda,(int)800*scaleFactor,(int)600*scaleFactor);
        }
    }else{
        m_view->zoomOut();
    }
}

QString MainWindow::getJpgTree(int X, int Y){
    QListWidgetItem *item =ui->listWidgetVolume->currentItem();
    int i = ui->listWidgetVolume->row(item);
    int num_bins = maxBinsInVolumeVector.at(i);

    allmybins_t * allmybeans = allmyDescpt[acVolumeName];
    perbinsize_t* curBinP = (*allmybeans)[num_bins];
    int k = ui->lineEditK->text().toInt();
    int lambda = ui->lineEditLambda->text().toInt();

    QStringList list = acVolumeName.split("/");
    QString base = list.at(list.size()-1).split(".").at(0);
    QString labImg=base+"_B"+QString::number(num_bins,10);
    labImg+="_K"+QString::number(k,10);
    labImg+="_L"+QString::number(lambda,10);
    labImg+="_"+QString::number(X,10);
    labImg+="x"+QString::number(Y,10);
    QString outdir="";
    for (int var=0 ; var < (list.size()-1); var++){
        outdir.append(list.at(var));
        outdir.append("/");
    }
    QString NewickFileName = outdir + "output"+"/"+ base +"/"+labImg+ ".new";

    return NewickFileName;
}

bool  MainWindow::verifyTocomputeTree(){
    if(bv::debug) printf( " \n verifyTocomputeTree \n");
    if (ui->listWidgetVolume->count()==0){
        QMessageBox::warning(this,"MATAV3","To compute a tree you must OPEN a volume.",
                       "ok");
        return false;
    }

    if (!ui->listWidgetVolume->currentItem()->isSelected()){
        QMessageBox::warning(this,"MATAV3","To compute a tree you must SELECT a volume.",
                       "ok");
        return false;
    }
    if (ui->lineEditBins<0 && ui->lineEditK<0 && ui->lineEditLambda<0){
        QMessageBox::warning(this,"MATAV3","The number of bins, k and lambda should be >= 0",
                       "ok");
        return false;
    }
    //verify if the tree is already computed

    allmybins_t * allmybeans = allmyDescpt[acVolumeName];
    if (!allmybeans){
        return true;
    }else   {
        std::cerr << "1" << std::endl;
        int num_bins = ui->lineEditBins->text().toInt();
        if (ui->lineEditBins->text().toInt()==0){
           QListWidgetItem *item =ui->listWidgetVolume->currentItem();
           int i = ui->listWidgetVolume->row(item);
           num_bins= maxBinsInVolumeVector.at(i);
        }
        perbinsize_t* curBinP = (*allmybeans)[num_bins];
        if (!curBinP){
            std::cerr << "2" << std::endl;
            return true;//tree for a BIN never computed
        }
        std::cerr << "3" << std::endl;
        params_t selectedTuple(ui->lineEditK->text().toInt(),ui->lineEditLambda->text().toInt());
        edgies_t& edgies = curBinP->fg_edgies;
        if(edgies.find(selectedTuple) == edgies.end()) {
            std::cerr << "4" << std::endl;
            return true;
        }else{
            QMessageBox::warning(this,"MATAV3","The component tree form this parameter is already computed",
                           "ok");
            return false;
        }
    }

    return true;
}

bool  MainWindow::verifyToDisplayTree(){

    if(bv::debug) printf( " \n verifyTocomputeTree \n");
    if (ui->listWidgetVolume->count()==0){
        QMessageBox::warning(this,"MATAV3","To display a tree you must OPEN a volume.",
                       "ok");
        return false;
    }

    if (!ui->listWidgetVolume->currentItem()->isSelected()){
        QMessageBox::warning(this,"MATAV3","To display a tree you must SELECT a volume.",
                       "ok");
        return false;
    }

    if (!ui->treeWidgetTree->currentItem()->isSelected()){
        QMessageBox::warning(this,"MATAV3","To display a tree you must SELECT a tree.",
                       "ok");
        return false;
    }
    return true;
}

/**
This function will save a volume with the voxel of a node
*/
void MainWindow::saveVolume(){
    if(bv::debug) printf("saveVolumeSegmentation() \n");
    QTime myTimer;
    myTimer.start();
    QListWidgetItem *item =ui->listWidgetVolume->currentItem();
    int i = ui->listWidgetVolume->row(item);
    int bins = maxBinsInVolumeVector.at(i);

    QList<QTreeWidgetItem *>  listT = ui->treeWidgetTree->selectedItems();
    QTreeWidgetItem  *itemT = listT.first();
    int lambda = itemT->text(2).toInt();
    int k  =itemT->text(3).toInt();
    int num_bins = itemT->text(1).toInt();
    if (bins==num_bins)
        num_bins = 0;

    allmybins_t * acAllmybeansP = allmyDescpt[acVolumeName];
    perbinsize_t* curBinP = (*acAllmybeansP)[num_bins];

    std::cout << "The FHTS for selected volume is ALREADY computed. num_bins=" << num_bins<< std::endl;
    historytree * http;
    historytree * httpUnsimplied;
    if (lambda!=0 || k!=0){
        http = acAllmybeansP->get_historyTree(num_bins,k,lambda);
        httpUnsimplied = curBinP->htfgp;
    }else{
        http =curBinP->htfgp;
    }

    historynode* nodeThres = m_tree->getHistoryNodeSelected();
    if (!nodeThres){
        QMessageBox::warning(this,"MATAV3","You need select ONE node to create a volume.",
                       "ok");
        return;
    }


    /*htfgp->verifyCoordenatesTree();
    if (lambda!=0){
        htfgp->updateAllVoxesTreePruned(fcts_tree->hroot,nodeThres);
    }else{

    }
    htfgp->verifyCoordenatesTree();
    */

    //open the DIG file to be filled up
    openNewVolume=false;
    if (!openDIGFile()){
        return 0;
    }
    DIGBasis InBasis;
    digvorfil.GetBasis(&InBasis);
    if(InBasis != DIGBasis_VORONOI) {
        std::cerr << "Input image must have DIGBasis DIGBasis_VORONOI." << std::endl;
        digvorfil.Close();
        return 0;
    }

    DIGDimensions InDims;
    const DIGDimensions* InDimsP = &InDims;
    digvorfil.GetDimensions(&InDimsP);
    int N_Z, N_Y, N_X;
    N_X = InDimsP->x;
    N_Y = InDimsP->y;
    N_Z = InDimsP->z;
    DIGDimensions OutDims;
    OutDims.x = N_X;
    OutDims.y = N_Y;
    OutDims.z = N_Z;
    if(bv::debug) std::cout << "N_Z: " << N_Z << ", N_Y: " << N_Y << " , N_X:  " << N_X << std::endl;

    const DIGDimensions* OutDimsP = &OutDims;
    unsigned channels;
    digvorfil.GetChannels(&channels);
    if(channels != 1) {
    digvorfil.Close();
    std::cerr << "Input image must be single channel (channels==1)" << std::endl;
    return 0;
    }

    DIGValueType vt;
    digvorfil.GetValueType(&vt);
    if (vt != DIGValueType_REAL) {
    digvorfil.Close();
    std::cerr << "Input image must be real-valued (not complex)" << std::endl;
    return 0;
    }
    DIGDataType dt;
    digvorfil.GetDataType(&dt);
    DIGDataFormat df;
    digvorfil.GetDataFormat(&df);
    DIGGrid InGrid;
    digvorfil.GetGrid(&InGrid);
    if(InGrid != DIGGrid_SC) {
    std::cerr << "input Voronoi must be on SC." << std::endl;
    digvorfil.Close();
    return 0;
    }
    DIGUnit unit;
    digvorfil.GetUnit(&unit); // as a nicety
    struct DIGSampling InSampX, InSampY, InSampZ;
    const struct DIGSampling* InSampXp = &InSampX;

    digvorfil.GetSamplingX(&InSampXp);
    const struct DIGSampling* InSampYp = &InSampY;
    digvorfil.GetSamplingY(&InSampYp);
    const struct DIGSampling* InSampZp = &InSampZ; // just leave as is... no problem
    digvorfil.GetSamplingZ(&InSampZp);

    struct DIGSampling OutSampX = {1.0, 0,   0  }; // double x, y, z
    struct DIGSampling OutSampY = {0,   1.0, 0  }; // double x, y, z
    struct DIGSampling OutSampZ = {0,   0,   1.0}; // double x, y, z
    unsigned int bufsize;
    unsigned int NoOfItems;
    bool errstat = 0;
    errstat |= digvorfil.SelectArraySet(0);

    errstat |= digvorfil.SelectArray(0);
    errstat |= digvorfil.GetArrayBufferSize(&bufsize);
    errstat |= digvorfil.GetArrayNoOfItems(&NoOfItems);

    if(errstat) {
        std::cerr << "Error pulling some info out of input DIGFile" << std::endl;
        digvorfil.Close();
        return 0;
    }

    // Creat a Dir to put the calculed tree
    QStringList list = acVolumeName.split("/");

    QString outdir="";
    for (int var=0 ; var < (list.size()-1); var++){
        outdir.append(list.at(var));
        outdir.append("/");
    }


    // Creat a Dir to put the save the images
    QString outputbasedir_str2 = outdir+"volumes";
    QDir outputbasedir2(outputbasedir_str2);
    if(! outputbasedir2.exists ()) {
        QDir d;
        if ( d.mkdir(outputbasedir_str2) ) { // this feels weird
            std::cerr << "Folder 'volumes'' SUCCESSFULLY created." << outputbasedir_str2.toAscii().constData() << std::endl;
        } else {
            std::cerr << "FAILED to create the folder 'volumes'." << outputbasedir_str2.toAscii().constData()<< std::endl;
        }
    } else {
        std::cerr << "The folder 'volumes' alread exist." << outputbasedir_str2.toAscii().constData() << std::endl;
    }




    QString base = outdir+"volumes/volumeSeg_"+list.at(list.size()-1).split(".").at(0);
    QString labImg=base+"_B"+QString::number(num_bins,10);
    labImg+="_K"+QString::number(k,10);
    labImg+="_L"+QString::number(lambda,10);



    QString outputbasedir_str = base;
    QDir outputbasedir(outputbasedir_str);
    if(! outputbasedir.exists ()) {
            QDir d;
            if ( d.mkdir(outputbasedir_str) ) { // this feels weird
                std::cerr << "SUCCESSFULLY created " << outputbasedir_str.toAscii().constData() << std::endl;
            } else {
                std::cerr << "FAILED to create " << outputbasedir_str.toAscii().constData()<< std::endl;
            }
    } else {
            std::cerr << "did not need to create " << outputbasedir_str.toAscii().constData() << std::endl;
    }



    QString volume_seg= base+"/volume"+QString::number((int)nodeThres->gray)+"_0.vol";
    int addend=0;
    bool lokingFile=true;
    while(QFile::exists(volume_seg)){
        addend++;
        volume_seg=  base+"/volume"+QString::number((int)nodeThres->gray)+"_"+QString::number(addend)+".vol";
        QString qm= "ADDEND" +volume_seg +" ok";
        ui->statusBar->showMessage(qm,1000);
    }


    QString volume_vor=  base+"/volume"+QString::number((int)nodeThres->gray)+"_"+QString::number(addend)+".voronoi";
    std::string OutFileName = volume_vor.toAscii().constData();
    QString command = "../../historyTree/xmippToDIGFile/xmippToDIGFile 0 "+volume_vor+" "+volume_seg +" -1 0 0";
    QString command_rm = "rm "+volume_vor;

    //QString volume_map=  base+"/volume"+QString::number((int)nodeThres->gray)+"_"+QString::number(addend)+".mrc";
    //QString commandMap = "xmipp_image_convert -i "+volume_seg+" -o "+volume_map;

    if(bv::debug) printf("..opening DIGfile\n");
    DIGFile OutFile;
    if(0 != OutFile.Open(
               OutFileName.c_str(),  // const char* file name
               "GENAPP",             // const char* schema (maybe define own)
               "slices ",  // const char* title
               "v3d2rec output",   // const char* type
               1,                    // unsigned int channels
               // DIGEndian_LITTLE,     // DIGEndian endian
               DIGValueType_REAL,    // DIGValueType value_type
               DIGDataType_FLOAT,   // DIGDataType data_type
               DIGDataFormat_BINARY,// DIGDataFormat data_format
               DIGGrid_SC,           // DIGGrid grid
               DIGBasis_VORONOI,     // DIGBasis basis
               DIGUnit_UNSPECIFIED, // DIGUnit unit
               OutDimsP,  // &xyz     // const DIGDimensions* dimensions;
               &OutSampX, // &x-vec   // const DIGSampling* samplingX;
               &OutSampY, // &y-vec   // const DIGSampling* samplingY;
               &OutSampZ, // &z-vec   // const DIGSampling* samplingZ;
               "Contains a single 3D array",  // char* comm
               "NO_other_unit"       // char* other_unit;
               )) {
        std::cerr << "Unable to open output file: " << OutFileName << std::endl;
        digvorfil.Close();
        return 0;
    }


    //Setup to create the file from the HistoryTree
    if(0 != OutFile.AppendArraySet(
                 "CONV", // ??
                 OutFileName.c_str(),
                 "from 3d image",
                 "only array set for this reconstruction"
                 )
     )
    {
      std::cerr << "Error appending an array set to the file" << std::endl;
      return 0;
    }
    cout <<" .....DIGFile is created: [" << (myTimer.elapsed()) <<" ms]" << std::endl;
    int tim = myTimer.elapsed();
    float *outbufp;
    try {
        outbufp = new float[(N_X* N_Y*N_Z)];
    } catch(std::bad_alloc) {
        std::cerr << "Heap memory exhausted " << std::endl;
    }

    //-----------------------------------------
    //fill up the image with the small grayness
    //-----------------------------------------
    if(bv::debug) printf("..get the output file ready \n");
    float *outp = outbufp;
    float *outp2 = outbufp;
    for(int z = 0; z < N_Z; z++) {
        for(int y = 0; y < N_Y; y++) {
            for(int x = 0; x < N_X; x++) {
                //outbufp[z*(N_Y*N_Z)+(y*N_Y+x)] = 0.0;
                outp[z*((N_Y)*(N_X))+(y*(N_X)+x)] =(float) http->hroot->gray;
                //outbufp[z*((N_Y-1)*(N_Z-1))+(y*(N_Y-1)+x)] = nodeThres->gray;
            }
        }
    }

    // get the label unpruned
    std::queue<historynode*> s;
    dimIndexType min_z = - static_cast<dimIndexType>(N_Z/2);
    dimIndexType min_y = - static_cast<dimIndexType>(N_Y/2);
    dimIndexType min_x = - static_cast<dimIndexType>(N_X/2);
    dimIndexType max_z = static_cast<dimIndexType>(N_Z-1)/2;
    dimIndexType max_y = static_cast<dimIndexType>(N_Y-1)/2;
    dimIndexType max_x = static_cast<dimIndexType>(N_X-1)/2;
    dimIndexType z,y,x;

   if(bv::debug) printf(" Get the information from the FCTS  <%d>\n",lambda);

   // FOR SIMPLIFIED TREE
   if (lambda!=0){
        // SHOW INFORMATION OF UNSIMPLIFIED
        if (true){
            std::vector<coord> myCoords;
            s.push(nodeThres);
            int contador=0;
            int tirados =0;

            while(!(s.empty())) {
                historynode* np_; // node pointer
                np_ = s.front();
                s.pop();
                if(bv::debug) printf("np_->listCoord.size()=%d,  \n",np_->listCoord.size());
                for(int t = 0; t < np_->listCoord.size(); t++) {
                    myCoords.push_back(np_->listCoord[t]);
                    z=	np_->listCoord[t].zpos;
                    y=	np_->listCoord[t].ypos;
                    x=	np_->listCoord[t].xpos;

                    int min_z = static_cast<dimIndexType>(N_Z/2);
                    int min_y = static_cast<dimIndexType>(N_Y/2);
                    int min_x =  static_cast<dimIndexType>(N_X/2);

                    if (x>N_X) if(bv::debug) printf("\n+++++ %lld %d %d\n ",np_->gray,x+min_x,min_x);
                    if (y>N_Y) if(bv::debug) printf("\n+++++ %lld %d %d\n ",np_->gray,y+min_y,min_y);
                    if (z>N_Z) if(bv::debug) printf("\n+++++ %lld %d %d\n ",np_->gray,z+min_z,min_z);
                    z= z+min_z;
                    y= y+min_y;
                    x= x+min_x;

                    if ((z*(N_Y*N_X)+(y*N_X)+x)>(N_Y*N_X*N_Z)){
                        if(bv::debug) printf("\n np_->listCoord.size()=%d (%d,%d,%d)\n ",np_->listCoord.size(),x,y,z);
                    }else{
                        outbufp[z*(N_Y*N_X)+(y*N_X)+x]= (float) np_->gray;
                    }
                    contador++;
                    tirados++;
                }
                if (np_->hfirstchild!=0){
                    historynode* fc =np_->hfirstchild;
                    s.push(fc);
                    while(fc->hnext){
                        historynode* axu = fc->hnext;
                        s.push(axu);
                        fc = axu;
                    }
                }
            }
            if(bv::debug) printf("contador=%d geral=%d \n ",contador,(N_Y*N_X*N_Z));
            /*std::vector<coord> myCoords;
            if(bv::debug) printf("..getting the information from the Tree (%f)  [lambda!=0]\n",(float)nodeThres->gray);
            s.push(nodeThres);
            int contador=0;
            int tirados =0;

            while(!(s.empty())) {
                historynode* np_; // node pointer
                np_ = s.front();
                if(bv::debug) printf("+++>>>[%d,%d]\n",np_->gray, np_->listCoord.size());

                historynode* clone = http->cloneInUnsTree(np_,httpUnsimplied,k); // child pointer
                if (clone==np_)
                    if(bv::debug) printf("      DIDN'T GET THE CLONE IN THE ORIGINAL TREE\n");

                s.pop();

                // Since a node does not have the voxel of all their descendent, we need go over the denscendant to get the voxels.
                std::stack<historynode*> nodVxs;
                nodVxs.push(clone);
                while(!(nodVxs.empty())) {
                    historynode* np_nodVxs; // node pointer
                    historynode* chp_nodVxs; // child pointer
                    np_nodVxs = nodVxs.top();
                    nodVxs.pop();
                    //print the voxels for np_nodVxs
                    if(bv::debug) printf("\n (np_nodVxs=%d,%d)",np_nodVxs->gray,np_nodVxs->listCoord.size());
                    for(int t = 0; t < np_nodVxs->listCoord.size(); t++) {
                        if(bv::debug) printf("%d",t);
                        myCoords.push_back(np_nodVxs->listCoord[t]);
                        z=	np_nodVxs->listCoord[t].zpos;
                        y=	np_nodVxs->listCoord[t].ypos;
                        x=	np_nodVxs->listCoord[t].xpos;

                        int min_z = static_cast<dimIndexType>(N_Z/2);
                        int min_y = static_cast<dimIndexType>(N_Y/2);
                        int min_x =  static_cast<dimIndexType>(N_X/2);

                        z= z+min_z;
                        y= y+min_y;
                        x= x+min_x;

                        if ((z*(N_Y*N_X)+(y*N_X)+x)>(N_Y*N_X*N_Z)){
                            if(bv::debug) printf("\n np_nodVxs->listCoord.size()=%d (%d,%d,%d)\n ",np_nodVxs->listCoord.size(),x,y,z);
                        }else{
                            if(bv::debug) printf("\n (:-P) (%d,%d,%d)= %f ",x,y,z,(float)np_nodVxs->gray);
                            outbufp[z*(N_Y*N_X)+(y*N_X)+x]= (float) np_nodVxs->gray;
                        }
                        contador++;
                        tirados++;
                    }

                    if(bv::debug) printf("\n");
                    if (np_nodVxs->hfirstchild!=0){
                        chp_nodVxs = (np_nodVxs->hfirstchild);
                        if(bv::debug) printf(" {%d,%d} ",np_nodVxs->gray,np_nodVxs->listCoord.size());
                        if (chp_nodVxs->listCoord.size()>=k)
                            nodVxs.push(chp_nodVxs);
                        historynode* pnc_nodVxs = chp_nodVxs->hnext;
                        while(pnc_nodVxs!= 0) {
                            if(bv::debug) printf(" {%d,%d} ",pnc_nodVxs->gray,pnc_nodVxs->listCoord.size());
                            if (pnc_nodVxs->numberVoxel>=k)
                                nodVxs.push(pnc_nodVxs);
                            chp_nodVxs = pnc_nodVxs->hnext;
                            pnc_nodVxs = chp_nodVxs;
                        }
                    }
                }// finish the loop to get the descendent voxels.




                if(np_->sizeNewIPCD!=0) {
                    if(bv::debug) printf("[child=%lld  np_->sizeNewIPCD=%d] \n",np_->newIPCD[0]->gray,np_->newIPCD[0]->sizeNewIPCD);
                    s.push(np_->newIPCD[0]);
                    for (int i=1; i< np_->sizeNewIPCD;i++){
                        if(bv::debug) printf("[child=%lld  np_->sizeNewIPCD=%d]\n",np_->newIPCD[i]->gray,np_->newIPCD[i]->sizeNewIPCD);
                        s.push(np_->newIPCD[i]);
                    }
                }
            }*/
            if(bv::debug) printf("contador=%d geral=%d \n ",contador,(N_Y*N_X*N_Z));
        //SHOW SIMPLIFIED TREE
        }else{
            //htfgp->updateAllVoxesTreePruned(fcts_tree->hroot,nodeThres);
            //  s.push(http->hroot);
            historynode* r = http->hroot;
            historynode* m;
            if (http->isCriticalVertex(http->hroot)){
                m = r;
                if (lambda==0){
                    r = fcts_tree->immediateCriticalVertice(r);
                }else{
                    while (	(r->deltaLimit <= http->delta1) && (r->hfirstchild!=0) ){
                        r = r->successor;
                    }
                }
                //s.push(r);
                if(bv::debug) printf("\n The root is a CV. The m is= (%lld) \n",http->hroot->gray);
            }else{
                r = http->immediateCriticalVertice(http->hroot);
                if (lambda==0){
                    historynode* r = fcts_tree->immediateCriticalVertice(r);
                }else{
                    while (	(r->deltaLimit <= http->delta1) && (r->hfirstchild!=0) ){
                        r = r->successor;
                    }
                }
                //s.push(r);
                if(bv::debug) printf("\n The root is NOT a CV. The m is= (%lld) \n",(http->immediateCriticalVertice(http->hroot))->gray);
            }

            int contador=0;
            int tirados =0;
            float val= -1.0;
            s.push(nodeThres);
            while(!(s.empty())) {
                historynode* np_; // node pointer
                historynode* chp_; // child pointer
                np_ = s.front();
                if(bv::debug) printf(" sizeNewIPCD %d  gray=%lld (%x,%x,%lld) \n",np_->sizeNewIPCD,np_->gray,np_,nodeThres,nodeThres->gray);

                s.pop();
                //if (np_== nodeThres){
                    if(bv::debug)printf("THE GRAYNESS OF THE SELECTED NODE IS=%d - np_->listCoord.size()=%d <%d> \n",np_->gray,np_->listCoord.size(),nodeThres->listCoord.size());
                    for(int t = 0; t < np_->listCoord.size(); t++) {
                        z=	(N_Z/2)+ np_->listCoord[t].zpos;
                        y=	(N_Y/2)+ np_->listCoord[t].ypos;
                        x=	(N_X/2)+ np_->listCoord[t].xpos;
                        outbufp[z*(N_Y*N_Z)+(y*N_Y+x)]= (float) np_->gray;
                        contador++;
                        tirados++;
                        //outp[z*(N_Y*N_Z)+(y*N_Y+x)]= 255;
                    }
                //}
                if(bv::debug) printf("contador=%d contadorGeral=%d level=%f \n\n ",contador,tirados,(float) np_->gray);
                 contador=0;
                if (lambda!=0){
                    if(np_->sizeNewIPCD!=0) {
                        //if(bv::debug) printf("[child=%lld  np_->sizeNewIPCD=%d]",np_->newIPCD[0]->gray,np_->newIPCD[0]->sizeNewIPCD);
                        s.push(np_->newIPCD[0]);
                        for (int i=1; i< np_->sizeNewIPCD;i++){
                            //if(bv::debug) printf("[child=%lld  np_->sizeNewIPCD=%d]",np_->newIPCD[i]->gray,np_->newIPCD[i]->sizeNewIPCD);
                            s.push(np_->newIPCD[i]);
                        }
                    }
                }else{
                    if (np_->hfirstchild!=0){
                        historynode* fc =np_->hfirstchild;
                        s.push(fc);
                        while(fc->hnext){
                                historynode* axu = fc->hnext;
                                s.push(axu);
                                fc = axu;
                        }
                    }
                }
            }
        }
    //FOR UNSIMPLIFIED TREE
    }else{
        std::vector<coord> myCoords;
        if(bv::debug) printf("..getting the information from the Tree (%f)\n,",(float)nodeThres->gray);
        s.push(nodeThres);
        int contador=0;
        int tirados =0;

        while(!(s.empty())) {
            historynode* np_; // node pointer
            np_ = s.front();
            s.pop();
            printf("np_->listCoord.size()=%d, \n ",np_->listCoord.size());
            for(int t = 0; t < np_->listCoord.size(); t++) {
                myCoords.push_back(np_->listCoord[t]);
                z=	np_->listCoord[t].zpos;
                y=	np_->listCoord[t].ypos;
                x=	np_->listCoord[t].xpos;

                int min_z = static_cast<dimIndexType>(N_Z/2);
                int min_y = static_cast<dimIndexType>(N_Y/2);
                int min_x =  static_cast<dimIndexType>(N_X/2);

                if (x>N_X) if(bv::debug) printf("\n+++++ %lld %d %d\n ",np_->gray,x+min_x,min_x);
                if (y>N_Y) if(bv::debug) printf("\n+++++ %lld %d %d\n ",np_->gray,y+min_y,min_y);
                if (z>N_Z) if(bv::debug) printf("\n+++++ %lld %d %d\n ",np_->gray,z+min_z,min_z);
                z= z+min_z;
                y= y+min_y;
                x= x+min_x;

                if ((z*(N_Y*N_X)+(y*N_X)+x)>(N_Y*N_X*N_Z)){
                    if(bv::debug) printf("\n np_->listCoord.size()=%d (%d,%d,%d)\n ",np_->listCoord.size(),x,y,z);
                }else{
                    outbufp[z*(N_Y*N_X)+(y*N_X)+x]= (float) np_->gray;
                    if(bv::debug) printf("np_->listCoord.size()=%d (%d,%d,%d) %f\n ",np_->listCoord.size(),x,y,z,((float)np_->gray));
                }
                contador++;
                tirados++;
            }
            if (np_->hfirstchild!=0){
                historynode* fc =np_->hfirstchild;
                s.push(fc);
                while(fc->hnext){
                    historynode* axu = fc->hnext;
                    s.push(axu);
                    fc = axu;
                }
            }
        }
        printf("contador=%d geral=%d \n ",contador,(N_Y*N_X*N_Z));

    }

    if(bv::debug) printf("..putting the information in the output \n");

    if(0 != OutFile.AppendArray( 	// how much is determined by parameters to Open() above
                1, // enumeration number
                "z-slice through 3D image",
                static_cast<void*>(outbufp))
       ) {
      std::cerr << "Error appending the array to the file" << std::endl;
      std::cerr << "Quitting." << std::endl;
      digvorfil.Close();
      OutFile.Close();
      //delete [] outbufp;
      return 0;
    }else{
        std::cout << "Salve the file" << std::endl;
    }
    cout <<" .....Fit the DIGFile: [" << (myTimer.elapsed()-tim) <<" ms]" << std::endl;
    tim = myTimer.elapsed();

    std::cout << "Successful termination! "<< std::endl;
    OutFile.Close();
    //delete [] outbufp;
    cout << "command:: " << command.toAscii().constData() <<"\n";
    system(command.toAscii().constData()) ;
    cout <<" .....command Vol [" << (myTimer.elapsed()-tim) <<" ms]" << std::endl;
    tim = myTimer.elapsed();

    //cout << "commandMap:: "<< commandMap.toAscii().constData() <<"\n";
    //system(commandMap.toAscii().constData()) ;
    //cout <<" .....command Map [" << (myTimer.elapsed()-tim) <<" ms]" << std::endl;
    system(command_rm.toAscii().constData());
    cout <<" .....Final [" << (myTimer.elapsed()    ) <<" ms]" << std::endl;
    return 0;

}

void MainWindow::showDescendant(){

    historynode* nodeThres = m_tree->getHistoryNodeSelected();
    Node * nod = m_tree->getNodeSelected();
    if (!nodeThres){
        QMessageBox::warning(this,"MATAV3","You need select ONE node to create a volume.",
                       "ok");
        return;
    }

    QListWidgetItem *item =ui->listWidgetVolume->currentItem();
    int i = ui->listWidgetVolume->row(item);
    int bins = maxBinsInVolumeVector.at(i);

    QList<QTreeWidgetItem *>  listT = ui->treeWidgetTree->selectedItems();
    QTreeWidgetItem  *itemT = listT.first();
    int lambda = itemT->text(2).toInt();
    int k  =itemT->text(3).toInt();
    int num_bins = itemT->text(1).toInt();
    if (bins==num_bins)
        num_bins = 0;

    allmybins_t * acAllmybeansP = allmyDescpt[acVolumeName];
    perbinsize_t* curBinP = (*acAllmybeansP)[num_bins];

    std::cout << "The FHTS for selected volume is ALREADY computed. num_bins=" << num_bins<< std::endl;
    historytree * http;
    historytree * httpUnsimplied;
    if (lambda!=0 || k!=0){
        http = acAllmybeansP->get_historyTree(num_bins,k,lambda);
        httpUnsimplied = curBinP->htfgp;
    }else{
        http =curBinP->htfgp;
    }

    m_tree->plotDescendant(nodeThres, http, nod);

}

void MainWindow::computeOcuppance(){

    historynode* nodeThres = m_tree->getHistoryNodeSelected();
    if (!nodeThres){
        QMessageBox::warning(this,"MATAV3","You need select ONE node to create a volume.",
                       "ok");
        return;
    }

    QListWidgetItem *item =ui->listWidgetVolume->currentItem();
    int i = ui->listWidgetVolume->row(item);
    int bins = maxBinsInVolumeVector.at(i);

    QList<QTreeWidgetItem *>  listT = ui->treeWidgetTree->selectedItems();
    QTreeWidgetItem  *itemT = listT.first();
    int lambda = itemT->text(2).toInt();
    int k  =itemT->text(3).toInt();
    int num_bins = itemT->text(1).toInt();
    if (bins==num_bins)
        num_bins = 0;

    allmybins_t * acAllmybeansP = allmyDescpt[acVolumeName];
    perbinsize_t* curBinP = (*acAllmybeansP)[num_bins];

    std::cout << "The FHTS for selected volume is ALREADY computed. num_bins=" << num_bins<< std::endl;
    historytree * http;
    historytree * httpUnsimplied;
    if (lambda!=0){
        http = acAllmybeansP->get_historyTree(num_bins,k,lambda);
        httpUnsimplied = curBinP->htfgp;
    }else{
        http =curBinP->htfgp;
        float val = http->computeOccuppanceTree(nodeThres);
        val = val/nodeThres->numberVoxel;
        QString show = QString::number((double)val,'g',5);
        QMessageBox::warning(this,"MATAV3","The ocuppance is: "+show,"OK");
    }

    //m_tree->plotDescendant(nodeThres, http, nod);

}
