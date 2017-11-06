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

// DIG.h,v 1.2 2005/07/27 16:39:05 jdubowy Exp

#ifndef DIG_H
#define DIG_H

#include <cstring>

//#include <DIGEndian.h>
#include <DIGFile/DIGEndian.h>

#if defined _MSC_VER

typedef signed char    DIGint8;
typedef signed __int16 DIGint16;
typedef signed __int32 DIGint32;
typedef signed __int64 DIGint64;

typedef unsigned char    DIGuint8;
typedef unsigned __int16 DIGuint16;
typedef unsigned __int32 DIGuint32;
typedef unsigned __int64 DIGuint64;

typedef float        DIGreal32;
typedef double       DIGreal64;
typedef long double  DIGreal80;

#if defined _M_ALPHA //Defined for DEC ALPHA platforms.

#elif defined _M_IX86 //Defined for x86 processors.

#elif defined _M_MPPC //Defined for Power Macintosh platforms.

#elif defined _M_MRX000 //Defined for MIPS platforms.

#elif defined _M_PPC //Defined for PowerPC platforms.

#endif

#elif defined __GNUC__

typedef signed char      DIGint8;
typedef signed short     DIGint16;
typedef signed int       DIGint32;
typedef signed long long DIGint64;

typedef unsigned char      DIGuint8;
typedef unsigned short     DIGuint16;
typedef unsigned int       DIGuint32;
typedef unsigned long long DIGuint64;

typedef float        DIGreal32;
typedef double       DIGreal64;
typedef long double  DIGreal80;

#if defined __i386__ //Defined for x86 processors.


#endif

#endif


class DIGString
{
protected:

  char*   String;

public:

  DIGString()
  {
    String = NULL;
  }

  int Set(const char* pString)
  {
    if(String != NULL) {         // if not empty 
      delete String;             // destroy old buffer
    }

    if(pString == NULL) {        // if NULL as parameter
      String = NULL;             // stet to NULL
      return 0;
    }

    int len = strlen(pString);   // get length of string
    String = new char[len + 1];  // create new buffer
    strcpy(String, pString);     // copy string to buffer

    return 0;
  }

  int Get(const char** pString)
  {
    *pString = String;
    return 0;
  }

  const char* Get()
  {
    return String;
  }

  ~DIGString()
  {
    if(String != NULL) {
      delete String;
    }
  }
};

#endif // DIG_H
