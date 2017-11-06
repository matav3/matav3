/****************************************************************************
** Meta object code from reading C++ file 'TreeView.h'
**
** Created: Sat Oct 13 19:24:44 2012
**      by: The Qt Meta Object Compiler version 63 (Qt 4.8.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../src/treeViewer/TreeView.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'TreeView.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_TreeView[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       7,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: signature, parameters, type, tag, flags
      10,    9,    9,    9, 0x05,

 // slots: signature, parameters, type, tag, flags
      30,    9,   24,    9, 0x0a,
      37,    9,    9,    9, 0x0a,
      46,    9,    9,    9, 0x0a,
      56,    9,    9,    9, 0x0a,
      69,    9,    9,    9, 0x0a,
      85,    9,    9,    9, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_TreeView[] = {
    "TreeView\0\0zoomChanged()\0qreal\0zoom()\0"
    "zoomIn()\0zoomOut()\0zoomNormal()\0"
    "setHandCursor()\0setArrowCursor()\0"
};

void TreeView::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        TreeView *_t = static_cast<TreeView *>(_o);
        switch (_id) {
        case 0: _t->zoomChanged(); break;
        case 1: { qreal _r = _t->zoom();
            if (_a[0]) *reinterpret_cast< qreal*>(_a[0]) = _r; }  break;
        case 2: _t->zoomIn(); break;
        case 3: _t->zoomOut(); break;
        case 4: _t->zoomNormal(); break;
        case 5: _t->setHandCursor(); break;
        case 6: _t->setArrowCursor(); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData TreeView::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject TreeView::staticMetaObject = {
    { &QGraphicsView::staticMetaObject, qt_meta_stringdata_TreeView,
      qt_meta_data_TreeView, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &TreeView::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *TreeView::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *TreeView::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_TreeView))
        return static_cast<void*>(const_cast< TreeView*>(this));
    return QGraphicsView::qt_metacast(_clname);
}

int TreeView::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QGraphicsView::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 7)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 7;
    }
    return _id;
}

// SIGNAL 0
void TreeView::zoomChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 0, 0);
}
QT_END_MOC_NAMESPACE
