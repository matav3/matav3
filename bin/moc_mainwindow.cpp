/****************************************************************************
** Meta object code from reading C++ file 'mainwindow.h'
**
** Created: Sun Nov 11 02:08:48 2012
**      by: The Qt Meta Object Compiler version 63 (Qt 4.8.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../mainwindow.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'mainwindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_MainWindow[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
      17,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      12,   11,   11,   11, 0x08,
      33,   11,   28,   11, 0x08,
      47,   11,   11,   11, 0x08,
      74,   66,   11,   11, 0x08,
     103,   11,   11,   11, 0x08,
     117,   11,   11,   11, 0x08,
     131,   11,   11,   11, 0x08,
     142,   11,   11,   11, 0x08,
     153,   11,   11,   11, 0x08,
     173,   11,   11,   11, 0x08,
     195,   11,   11,   11, 0x08,
     217,   11,   28,   11, 0x08,
     239,   11,   28,   11, 0x08,
     261,   11,   11,   11, 0x08,
     278,   11,   11,   11, 0x08,
     293,   11,   11,   11, 0x08,
     311,   11,   11,   11, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_MainWindow[] = {
    "MainWindow\0\0numberChanged()\0bool\0"
    "openDIGFile()\0updateVolumeList()\0"
    "bin,par\0updateTreeList(int,params_t)\0"
    "computeTree()\0displayTree()\0zoomPlus()\0"
    "zoomLess()\0treeViewOptionJPG()\0"
    "treeViewOptionQPlot()\0treeViewOptionPSPDF()\0"
    "verifyTocomputeTree()\0verifyToDisplayTree()\0"
    "showQTreeLabel()\0showQTreeBox()\0"
    "showQTreeVoxels()\0showQTreeCV()\0"
};

void MainWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        MainWindow *_t = static_cast<MainWindow *>(_o);
        switch (_id) {
        case 0: _t->numberChanged(); break;
        case 1: { bool _r = _t->openDIGFile();
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = _r; }  break;
        case 2: _t->updateVolumeList(); break;
        case 3: _t->updateTreeList((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< params_t(*)>(_a[2]))); break;
        case 4: _t->computeTree(); break;
        case 5: _t->displayTree(); break;
        case 6: _t->zoomPlus(); break;
        case 7: _t->zoomLess(); break;
        case 8: _t->treeViewOptionJPG(); break;
        case 9: _t->treeViewOptionQPlot(); break;
        case 10: _t->treeViewOptionPSPDF(); break;
        case 11: { bool _r = _t->verifyTocomputeTree();
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = _r; }  break;
        case 12: { bool _r = _t->verifyToDisplayTree();
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = _r; }  break;
        case 13: _t->showQTreeLabel(); break;
        case 14: _t->showQTreeBox(); break;
        case 15: _t->showQTreeVoxels(); break;
        case 16: _t->showQTreeCV(); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData MainWindow::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject MainWindow::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_MainWindow,
      qt_meta_data_MainWindow, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &MainWindow::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *MainWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *MainWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_MainWindow))
        return static_cast<void*>(const_cast< MainWindow*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int MainWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 17)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 17;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
