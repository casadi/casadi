/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#include "vertcat.hpp"
#include "../stl_vector_tools.hpp"
#include <cassert>
#include <iterator>


#include "slicer.hpp"

using namespace std;

namespace CasADi{

//Slicer::Slicer(char s_) { std::string s;s=s_;p=new SlicerPrimitiveString(s);}
Slicer::Slicer(const char *s_) 			{ p=new SlicerPrimitiveString(std::string(s_));}
Slicer::Slicer(const std::string &s_)		{ p=new SlicerPrimitiveString(s_);}
Slicer::Slicer(std::vector<int> &i_) 		{ p=new SlicerPrimitiveList(i_);}
Slicer::Slicer(int i)				{ p=new SlicerPrimitiveSingle(i);}
Slicer::Slicer(const SlicerPrimitiveAll &p_)	{ p=new SlicerPrimitiveAll();}
Slicer::Slicer(const SlicerPrimitiveFromTo &p_)	{ p=new SlicerPrimitiveFromTo(p_);}
Slicer::Slicer(const Slicer &s)	{
	if (s.p->type==SINGLE) {
		p=new SlicerPrimitiveSingle(*static_cast<SlicerPrimitiveSingle*>(s.p));
	} else if (s.p->type==LIST) {
		p=new SlicerPrimitiveList(*static_cast<SlicerPrimitiveList*>(s.p));
	} else if (s.p->type==ALL) {
		p=new SlicerPrimitiveAll(*static_cast<SlicerPrimitiveAll*>(s.p));
 	} else if (s.p->type==FROMTO) {
		p=new SlicerPrimitiveFromTo(*static_cast<SlicerPrimitiveFromTo*>(s.p));
	} else if (s.p->type==STRING) {
		p=new SlicerPrimitiveString(*static_cast<SlicerPrimitiveString*>(s.p));
	} else if (s.p->type==UNSPECIFIED) {
		cout << "TROUBLE" << endl;
	}
}

Slicer::~Slicer()	{
	delete p;
}

vector<int>::iterator Slicer::begin() {return p->begin();}
vector<int>::iterator Slicer::end()  {return p->end();}
void Slicer::initialize(int end) {p->initialize(end);}



SlicerPrimitive::SlicerPrimitive(SlicerPrimitiveType type_): type(type_) {}
void SlicerPrimitive::initialize(int end) {}
vector<int>::iterator SlicerPrimitive::begin() {return ind.begin();}
vector<int>::iterator SlicerPrimitive::end() {return ind.end();}

SlicerPrimitiveList::SlicerPrimitiveList(std::vector<int> &i_): SlicerPrimitive(LIST) {ind = i_;}
void SlicerPrimitiveList::initialize(int end_) {}

SlicerPrimitiveFromTo::SlicerPrimitiveFromTo(int from_, int to_, int step_): SlicerPrimitive(FROMTO), from(from_), to(to_), step(step_) {}

void SlicerPrimitiveFromTo::initialize(int end) {
	if (end==0) {
		throw CasadiException("Attempt to slice with list size of 0.");
	}
	while (to<0) {to+=end+1;}
	while (from<0) {from+=end+1;}
	for (int i=from;i<to;i+=step) {
		ind.push_back(i);
	}
}

SlicerPrimitiveAll::SlicerPrimitiveAll(): SlicerPrimitiveFromTo(0,-1) {}

SlicerPrimitiveString::SlicerPrimitiveString(const std::string &s_): SlicerPrimitive(STRING), s(s_) {}

void SlicerPrimitiveString::initialize(int end) {
	if (s.compare(":")==0) {
		SlicerPrimitiveFromTo sl(0,end);
		sl.initialize(end);
		ind=sl.ind;
	} else {
		throw CasadiException("SlicePrimitiveString only accepts ':'.");
	}
}


SlicerPrimitiveSingle::SlicerPrimitiveSingle(int i): SlicerPrimitiveFromTo(i,i+1) {}




} // namespace CasADi
