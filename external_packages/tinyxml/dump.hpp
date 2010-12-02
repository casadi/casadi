/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#ifndef DUMP_HPP
#define DUMP_HPP

#include "tinyxml.h"

const char * getIndent( unsigned int numIndents );

/** \brief  same as getIndent but no "+" at the end */
const char * getIndentAlt( unsigned int numIndents );

int dump_attribs_to_stdout(TiXmlElement* pElement, unsigned int indent);

void dump_to_stdout( TiXmlNode* pParent, unsigned int indent = 0 );

/** \brief  load the named file and dump its structure to STDOUT */
void dump_to_stdout(const char* pFilename);


#endif