/*

  This file is part of DIG Library.
  DIG Library is a free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published
  by the Free Software Foundation; either version 2 of the License,
  or (at your option) any later version.
  DIG Library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with Foobar; if not, write to the Free Software Foundation,
  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

*/

// DIGFile.h,v 1.4 2007/06/07 20:48:35 jklukowska Exp


//defines added to enable large file support
//jklukowska, 06/07/2007
#ifndef _LARGEFILE_SOURCE
	#define _LARGEFILE_SOURCE     /* enable large file support  */
#endif
#ifndef _FILE_OFFSET_BITS
	#define _FILE_OFFSET_BITS 64     /* enable large file support  */
#endif


#ifndef DIG_FILE_H
#define DIG_FILE_H

#include <xercesc/dom/DOM.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/util/PlatformUtils.hpp>

//#include <DIG.h>
#include <DIGFile/DIG.h>
//#include "DOMTreeErrorReporter.h"
#include <DIGFile/DOMTreeErrorReporter.h>

#include <xercesc/framework/MemBufInputSource.hpp>

enum DIGEndian { 
  DIGEndian_LITTLE,
  DIGEndian_BIG
};

#define DIGEndianSize 2
extern const char* DIGEndianStr[DIGEndianSize];


enum DIGValueType {
  DIGValueType_REAL, 
  DIGValueType_COMPLEX
};

#define DIGValueTypeSize 2
extern const char* DIGValueTypeStr[DIGValueTypeSize];


enum DIGDataType {
  DIGDataType_BOOL, 
  DIGDataType_UCHAR, 
  DIGDataType_CHAR, 
  DIGDataType_USHORT, 
  DIGDataType_SHORT, 
  DIGDataType_UINT, 
  DIGDataType_INT, 
  DIGDataType_FLOAT, 
  DIGDataType_DOUBLE
};

#define DIGDataTypeSize 9
extern const char* DIGDataTypeStr[DIGDataTypeSize];


enum DIGDataFormat {
  DIGDataFormat_ASCII, 
  DIGDataFormat_BINARY
};

#define DIGDataFormatSize 2
extern const char* DIGDataFormatStr[DIGDataFormatSize];


enum DIGGrid {
  DIGGrid_SC, 
  DIGGrid_FCC1, 
  DIGGrid_FCC2, 
  DIGGrid_BCC1, 
  DIGGrid_BCC2, 
  DIGGrid_BCC3, 
  DIGGrid_BCC4, 
  DIGGrid_HEX1, 
  DIGGrid_HEX2
};

#define DIGGridSize 9
extern const char* DIGGridStr[DIGGridSize];


enum DIGBasis {
  DIGBasis_VORONOI, 
  DIGBasis_BLOB
};

#define DIGBasisSize 2
extern const char* DIGBasisStr[DIGBasisSize];


enum DIGUnit {
  DIGUnit_UNSPECIFIED, 
  DIGUnit_m, 
  DIGUnit_cm, 
  DIGUnit_mm, 
  DIGUnit_nm, 
  DIGUnit_micron, 
  DIGUnit_A, 
  DIGUnit_OTHER
};

#define DIGUnitSize 8
extern const char* DIGUnitStr[DIGUnitSize];


struct DIGDimensions {
	int x;
	int y;
	int z;
};

struct DIGSampling {
	double x;
	double y;
	double z;
};

///////////////////////////////////////////////////////////////////////////////
// DIGFile array set header
///////////////////////////////////////////////////////////////////////////////

class DIGFileMainHeader
{
public:
 	DIGString       Type;
 	DIGString       Title;
	unsigned int    Channels;
	//unsigned int    NoOfArrays;
	//unsigned int    NoOfSets;
	//DIGEndian       Endian;
	DIGValueType    ValueType;
	DIGDataType     DataType;
	DIGDataFormat   DataFormat;
	DIGGrid         Grid;
	DIGBasis        Basis;
	DIGUnit         Unit;
	DIGString       OtherUnit;
	DIGDimensions   Dimensions;
	DIGSampling     SamplingX;
	DIGSampling     SamplingY;
	DIGSampling     SamplingZ;
	DIGString       Comment;

public:
  DIGFileMainHeader()
  {
    // set defaults
	  Channels = 1;
	  //NoOfArrays = 0;
	  //NoOfSets = 0;
	  //Endian = DIGEndian_LITTLE;
	  ValueType = DIGValueType_REAL;
	  DataType = DIGDataType_DOUBLE;
	  DataFormat = DIGDataFormat_ASCII;
	  Grid = DIGGrid_SC;
	  Basis = DIGBasis_VORONOI;
	  Unit = DIGUnit_UNSPECIFIED;
	  Dimensions.x = 1;
	  Dimensions.y = 1;
	  Dimensions.z = 0;
	  SamplingX.x = 1;
	  SamplingX.y = 0;
	  SamplingX.z = 0;
	  SamplingY.x = 0;
	  SamplingY.y = 1;
	  SamplingY.z = 0;
	  SamplingZ.x = 0;
	  SamplingZ.y = 0;
	  SamplingZ.z = 1;
  };

  DIGFileMainHeader(
 	  const char*   pType,
 	  const char*   pTitle,
    unsigned int  pChannels, 
    //unsigned int  pArrays, 
    //unsigned int  pSets, 
    //DIGEndian     pEndian, 
    DIGValueType  pValueType, 
    DIGDataType   pDataType, 
    DIGDataFormat pDataFormat, 
    DIGGrid       pGrid, 
    DIGBasis      pBasis, 
    DIGUnit       pUnit, 
    DIGDimensions pDim, 
    DIGSampling   pSampX, 
    DIGSampling   pSampY, 
    DIGSampling   pSampZ, 
    const char*   pComment, 
    const char*   pOtherUnit
  )
  {
    Type.Set(pType);
    Title.Set(pTitle);
	  Channels = pChannels;
	  //NoOfArrays = pArrays;
	  //NoOfSets = pSets;
	  //Endian = pEndian;
	  ValueType = pValueType;
	  DataType = pDataType;
	  DataFormat = pDataFormat;
	  Grid = pGrid;
	  Basis = pBasis;
	  Unit = pUnit;
	  OtherUnit.Set(pOtherUnit);
	  Dimensions = pDim;
	  SamplingX = pSampX;
	  SamplingY = pSampY;
	  SamplingZ = pSampZ;
    Comment.Set(pComment);
  }; 
};

///////////////////////////////////////////////////////////////////////////////
// DIGFile array set header
///////////////////////////////////////////////////////////////////////////////

class DIGFileArraySetHeader
{
public:

	DIGString Type;
	DIGString Title;
	DIGString Parameters;
  DIGString Comment;
	
  DIGFileArraySetHeader() {};
  DIGFileArraySetHeader(const char* pType, const char* pTitle, const char* pParameters, const char* pComments = NULL)
  {
	  Type.Set(pType);
	  Title.Set(pTitle);
	  Parameters.Set(pParameters); 
    Comment.Set(pComments);
  };
};

///////////////////////////////////////////////////////////////////////////////
// DIGFile array header
///////////////////////////////////////////////////////////////////////////////

class DIGFileArrayHeader {

public:
	unsigned int EnumNo;
  DIGString    Comment;

  DIGFileArrayHeader() {};
  DIGFileArrayHeader(unsigned int pEnumNo, const char* pComment) {
    EnumNo = pEnumNo;
    Comment.Set(pComment);
  };
};



class ListElement 
{
public:
  void* Item;
  ListElement* NextListElement;

  inline ListElement(void* It)
  {
    Item = It;
    NextListElement = NULL;
  }
};

class List 
{
private:
  ListElement*  Head;

protected:
  unsigned int  Length;
  
public:
  List() {
    Head = NULL;
    Length = 0;
  }

  inline unsigned int getLength() {
    return Length;
  }

  inline void Append(ListElement* LElement) {
    ListElement** ListElementPtr = &Head;
    while((*ListElementPtr) != NULL) {
      ListElementPtr = &((*ListElementPtr)->NextListElement);
    }
    *ListElementPtr = LElement;
    Length++;
  }

  inline int getListElement(ListElement** ListElementPtr, unsigned int index) {

    *ListElementPtr = Head;
    if((*ListElementPtr) == NULL) {
      return -1;  // error empty list
    }

    for(unsigned int i = 0; i < index; i++) {
      *ListElementPtr = (*ListElementPtr)->NextListElement;
      if((*ListElementPtr) == NULL) {
        return -1;  // error
      }
    }
    return 0;
  }

  inline int removeListElement(ListElement** ListElementPtr, unsigned int index) {

    *ListElementPtr = Head;
    if((*ListElementPtr) == NULL) {
      return -1;  // error empty list
    }

    // first element of the list
    ListElement** prevPtr = &Head;

    for(unsigned int i = 0; i < index; i++) {

      prevPtr = &((*ListElementPtr)->NextListElement);
      *ListElementPtr = (*ListElementPtr)->NextListElement;
      if((*ListElementPtr) == NULL) {
        return -1;  // error
      }
    }

    *prevPtr = (*ListElementPtr)->NextListElement;
    (*ListElementPtr)->NextListElement = NULL;
    Length--;
    return 0;
  }
};

class PointerList: private List
{
  
public:
/*
  List() {
    Head = NULL;
    Length = 0;
  }
  */

  inline unsigned int getLength() {
    return Length;
  }

  inline void Append(void* Pointer) {
    ListElement* LI = new ListElement(Pointer);
    List::Append(LI);
  }

  inline int getListElement(void** Pointer, unsigned int index) {

    ListElement* LI;
    if(List::getListElement(&LI, index) != 0) {
      return -1; // error
    }
    *Pointer = LI->Item;
    return 0;
  }

  inline int removeListElement(void** Pointer, unsigned int index) {

    ListElement* LI;

    if(List::removeListElement(&LI, index) != 0) {
      return -1; // error
    }

    *Pointer = LI->Item;
    delete LI;
    return 0;
  }
};


class ArraySet 
{
private:
  PointerList            ArrayHeaderList;  // list of array heders in the set

public:
  unsigned int           FirstArrayOffset; // offset of first array in the set
  DIGFileArraySetHeader* ArraySetHeader;

  
  inline unsigned int getLength() {
    return ArrayHeaderList.getLength();
  }

  inline void Append(DIGFileArrayHeader* ArrayHeader) {
    ArrayHeaderList.Append(ArrayHeader);
  }

  inline int getArrayHeader(DIGFileArrayHeader** ArrayHeader, unsigned int index) {
    return ArrayHeaderList.getListElement((void**) ArrayHeader, index);
  }

  inline int removeArrayHeader(DIGFileArrayHeader** ArrayHeader, unsigned int index) {
    return ArrayHeaderList.removeListElement((void**) ArrayHeader, index);
  }

  inline int Release() {
    DIGFileArrayHeader* ArrayHeader;

    int Lenght = getLength();
    for(int i = 0; i < Lenght; i++) {
      // remove first item from the list
      if(removeArrayHeader(&ArrayHeader, 0) != 0) {
        return -1; // error
      }
      delete ArrayHeader;
    }
    return 0;
  }
};


class ArraySetList 
{
  PointerList            ASList;
  
public:

  inline unsigned int getLength() {
    return ASList.getLength();
  }

  inline void Append(ArraySet* ASet) {
    ASList.Append(ASet);
  }

  inline int getArraySet(ArraySet** ASet, unsigned int index) {
    return ASList.getListElement((void**) ASet, index);
  }

  inline int removeArraySet(ArraySet** ASet, unsigned int index) {
    return ASList.removeListElement((void**) ASet, index);
  }

  inline int Release() {
    ArraySet* ASet;

    int Lenght = getLength();
    for(int i = 0; i < Lenght; i++) {
      // remove first item from the list
      if(removeArraySet(&ASet, 0) != 0) {
        return -1; // error
      }
      ASet->Release();
    }
    return 0;
  }
};

class DIGFile
{
protected:

  // File Variables
  DIGString           SchemaName;         // XML schema name
  DIGFileMainHeader   MainHeader;         // main header

  unsigned int        NoOfSets;
  unsigned int        NoOfArrays;
  DIGEndian           Endian;             // File Endian

  FILE*               File;
  FILE*               TmpDataFile;

  int                   XMLHeaderLength;
  bool                  MainHeaderWritten;
  bool                  Dirty;

  ArraySetList          ArraySets;        // array sets list

  unsigned int          NoOfArrayItems;   // number of items in array
  unsigned int          ArrayBufferSize;  // array buffer size

  // Curent Array Set Variables
  ArraySet*             CurrenArraySet;

  // Cureny Array Variables
  DIGFileArrayHeader*   CurrentArrayHeader;
  int          CurrenArrayIndex; // index of curent array


  // XML Parsing and Writing

  bool                  ArraySetOpen;
  //int                   ArrayCounter;     // arrays counter

  // DOM 

  XercesDOMParser*      parser;
  DOMTreeErrorReporter* errReporter;

  DOMImplementation* impl;

  DOMDocument* Document;
  DOMElement*  RootElement;
  DOMElement*  MainHeaderElement;
  DOMElement*  ApplicationHeaderElement;
  DOMElement*  ArraySetHeaderElement;     // active array set header
  DOMElement*  ArrayHeaderElement;        // active array header

  DOMNodeList* ArraySetHeaderList;
  DOMNodeList* ArrayHeaderList;


  // XML Strings
  static bool XMLStringsInitialized;


  static const XMLCh* RootXMLStr;

  static const XMLCh* MainHeaderXMLStr;

  static const XMLCh* TitleXMLStr;
  static const XMLCh* TypeXMLStr;
  static const XMLCh* ChannelsXMLStr;
  static const XMLCh* NumberOfArraysXMLStr;
  static const XMLCh* NumberOfArraySetsXMLStr;
  static const XMLCh* EndianXMLStr;
  static const XMLCh* ValueTypeXMLStr;
  static const XMLCh* DataTypeXMLStr;
  static const XMLCh* DataFormatXMLStr;
  static const XMLCh* GridTypeXMLStr;
  static const XMLCh* BasisFunctionXMLStr;
  static const XMLCh* UnitXMLStr;
  static const XMLCh* OtherUnitXMLStr;
  static const XMLCh* CommentsXMLStr;
  static const XMLCh* DimensionsXMLStr;

  static const XMLCh* DimmensionXXMLStr;
  static const XMLCh* DimmensionYXMLStr;      
  static const XMLCh* DimmensionZXMLStr;

  static const XMLCh* SamplingRateXXMLStr;
  static const XMLCh* SamplingRateYXMLStr;
  static const XMLCh* SamplingRateZXMLStr;

  static const XMLCh* ArraySetHeaderXMLStr;

  //static const XMLCh* TypeXMLStr;
  //static const XMLCh* TitleXMLStr;
  static const XMLCh* ParametersXMLStr;

  static const XMLCh* ArrayHeaderXMLStr;
  static const XMLCh* EnumerationNumberXMLStr;

public:
  DIGFile();
  virtual ~DIGFile() { };

  int Open(const char* FileName);
  int Open(
    const char*           pFileName, 
    const char*           pSchema, 
 	  const char*           pTitle,
 	  const char*           pType,
	  unsigned int          pChannels,
	  //DIGEndian             pEndian,
	  DIGValueType          pValueType,
	  DIGDataType           pDataType,
	  DIGDataFormat         pDataFormat,
	  DIGGrid               pGrid,
	  DIGBasis              pBasis,
	  DIGUnit               pUnit,
	  const DIGDimensions*  pDimensions,
	  const DIGSampling*    pSamplingX,
	  const DIGSampling*    pSamplingY,
	  const DIGSampling*    pSamplingZ,
	  const char*           pComment,
	  const char*           pOtherUnit
  );

  int Close();

  int GetArrayBufferSize(unsigned int* pArrayBufferSize);
  int GetArrayNoOfItems(unsigned int* pArrayNoOfItems);
  int GetNoOfArraySets(unsigned int* pNoOfArraySets);
  int GetNoOfArrays(unsigned int* pNoOfArrays);

  // main header
  int SetComment(const char* pComment);
  int SetTitle(const char* pTitle);

  int GetTitle(const char** pTitle);
  int GetType(const char** pType);
  int GetChannels(unsigned int* pChannels);
  int GetEndian(DIGEndian* pEndian);
  int GetValueType(DIGValueType* pValueType);
  int GetDataType(DIGDataType* pDataType);
  int GetDataFormat(DIGDataFormat* pDataFormat);
  int GetGrid(DIGGrid* pGrid);
  int GetBasis(DIGBasis* pBasis);
  int GetUnit(DIGUnit* pUnit);
  int GetOtherUnit(const char** pOtherUnit);
  int GetDimensions(const DIGDimensions** pDimensions);
  int GetSamplingX(const DIGSampling** pSampling);
  int GetSamplingY(const DIGSampling** pSampling);
  int GetSamplingZ(const DIGSampling** pSampling);
  int GetComment(const char** pComment);


  //int AppendArraySet(DIGFileArraySetHeader* SetHeader);
  int AppendArraySet(
    const char* pType,
    const char* pTitle,
    const char* pParameters,
    const char* pComments
  );

   //int AppendArray(DIGFileArrayHeader* ArrayHeader, void* Buffer); 
  int AppendArray(
    unsigned int pEnumerationNumber,
    const char* pComment,
    const void* pBuffer
  );

  //int GetArraySetHeader(DIGFileArraySetHeader** SetHeader);
  int GetArraySetType(const char** pType);
  int GetArraySetTitle(const char** pTitle);
  int GetArraySetParameters(const char** pParameters);
  int GetArraySetComment(const char** pComments);

  //int GetArrayHeader(DIGFileArrayHeader** ArrayHeader);
  int GetArrayEnumNo(unsigned int* pEnumNo);
  int GetArrayComment(const char** pComment);


  int SelectArraySet(unsigned int pNo);
  int SelectArray(unsigned int pNo);

  int GetArrayData(void* pArrayData);

protected:

  unsigned int ComputeArrayNoOfItems();
  unsigned int ComputeArrayBufferSize();

  int SymbolToCode(char* Symbol, const char** SymbolTable, int NumberOfSymbols);

  int WriteArrayData(FILE* f, const void* ArrayData);

  int FindSrting(char* String, FILE* File);
  int GetXMLHeader(char** XMLBuff);

  // Writing Methods
  int WriteXMLHeader();
  void CreateRootElement();
  virtual int CreateApplicationHeader();
  int CreateMainHeader(DIGFileMainHeader* header);
  int AppendArraySetHeader(DIGFileArraySetHeader* SetHeader);
  int AppendArrayHeader(DIGFileArrayHeader* ArrayHeader);


  // Reading Methods
  int InitParser2(const MemBufInputSource& memBufIS);
  int ResetParser();
  int GetIntAttributeValue(DOMElement* Element, const XMLCh* AttributeName, int* value);
  int GetUintAttributeValue(DOMElement* Element, const XMLCh* AttributeName, unsigned int* value);
  int GetDoubleAttributeValue(DOMElement* Element, const XMLCh* AttributeName, double* value);
  int GetTxtAttributeValue(DOMElement* Element, const XMLCh* AttributeName, char** Text);
  int GetElementText(DOMElement* Element, char** Text);

  int GetMainHeader(DIGFileMainHeader* MainHeader);
  virtual int ParseApplicationHeader();
  int GetArraySetHeader(DIGFileArraySetHeader* ArraySetHeader, unsigned int index);
  int GetArrayHeader(DIGFileArrayHeader* ArrayHeader, unsigned int index);

  int ParseXML(const char* MemBuf);

  int CreateMainHeader();
};

#endif
