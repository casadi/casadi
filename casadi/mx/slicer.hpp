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

#ifndef SLICER_HPP
#define SLICER_HPP


/*
	These classes have some pretty strange design.

	The following decisions & facts led to the current structure:
	- we want to have one type 'Slicer' that can accept a broad range of types via implicit type conversion with the help of constructors.
	- each flavour of accepted type, uses its own algorithm. we want to use specialization: inheritance from a baseclass 'SlicerPrimitive'.
	- it appears that C++, when doing implicit type conversion, does not care to look to subclasses
	- so Slicer itself needs to reimplement all constructors.
	- Slicer should somehow store the 'SlicerPrimitive' objects
	- Storing them as 'SlicerPrimitive' destroys polymorphism, and storing 'SlicerPrimitive &' didn't work
	- So 'SlicerPrimitive' has to work with pointers
	- But weird things happened because of implicit cloning
	- Had to resort to ugly tricks: add type to 'SlicerPrimitive', have typechecking clone function in Slicer.


	todo: exploit sparsity if available
*/

namespace CasADi{

enum SlicerPrimitiveType {UNSPECIFIED, SINGLE, LIST, ALL, STRING, FROMTO};

/**
 \brief Slicer primitive is the base class for various slicer algorithms.
*/
class SlicerPrimitive {
public:

//! Default constructor
SlicerPrimitive(SlicerPrimitiveType type=UNSPECIFIED);
/**
	\brief Performs initialization.
	Not done in constructor because the value of 'end' may be needed.
*/
virtual void initialize(int end=-1);
//!	\brief Let SlicerPrimitive return iterators
vector<int>::iterator begin();
//!	\brief Let SlicerPrimitive return iterators
vector<int>::iterator end();
//! Return resultant size of slice
int size();

//! index of end

int endind;
//! \brief vector of indices
vector<int> ind;

//! Return slice index at i
int operator()(int i);

enum SlicerPrimitiveType type;

};

/**
 \brief Slicer primitive that takes a vector of integers
*/

class SlicerPrimitiveList : public SlicerPrimitive {
public:
	/** \brief  Constructor */
	SlicerPrimitiveList(std::vector<int> &i);

	void initialize(int end=-1);

private:


};

/**
 \brief Slicer primitive that implements basic (from, to, step)

 If to is negative, (length + 1) is added.
 In this way, SlicerPrimitiveFromTo(0,-1) behaves the same as SlicerPrimitiveAll()
*/

class SlicerPrimitiveFromTo : public SlicerPrimitive {
public:
	/** \brief  Constructor */
	SlicerPrimitiveFromTo(int from, int to, int step=1);
	void initialize(int end=-1);
private:
	int from;
	int to;
	int step;
	
};

/**
 \brief Slicer primitive that selects all. (like ':' in matlab)
*/
class SlicerPrimitiveAll : public SlicerPrimitiveFromTo {
public:
	/** \brief  Constructor */
	SlicerPrimitiveAll();

private:

};

/**
 \brief Slicer primitive that selects all. (like ':' in matlab)
*/
class SlicerPrimitiveString : public SlicerPrimitive {
public:
	/** \brief  Constructor */
	SlicerPrimitiveString(const std::string &s);
	void initialize(int end=-1);

private:
	const std::string &s;
};


/**
 \brief Slicer primitive that takes a single interger
*/

class SlicerPrimitiveSingle : public SlicerPrimitiveFromTo {
public:
	/** \brief  Constructor */
	SlicerPrimitiveSingle(int i);

private:

};



/** \brief Generates indices iterators for array slices.

  Use cases:

  SlicerPrimitive understands:
	SlicerPrimitiveAll()  (equivalent for ":" in matlab)
	SlicerPrimitiveFromTo()
	":"
	single index (int)
	std::vector of indices

  \author Joris Gillis 
  \date 2010
*/
class Slicer {


public:

//Slicer(char s_);
Slicer(const char *s_);
Slicer(const std::string &s);
Slicer(std::vector<int> &i);
Slicer(int i);
Slicer(const SlicerPrimitiveAll &p_);
Slicer(const SlicerPrimitiveFromTo &p_);
//! Copy constructor
Slicer(const Slicer &s);

void initialize(int end=-1);
vector<int>::iterator begin();
vector<int>::iterator end();

//! Return resultant size of slice
int size();

//! Return slice index at i
int operator()(int i);

~Slicer();

protected:
	SlicerPrimitive *p;
};


} // namespace CasADi

#endif // SLICER_HPP
