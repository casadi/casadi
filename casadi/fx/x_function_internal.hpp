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

#ifndef X_FUNCTION_INTERNAL_HPP
#define X_FUNCTION_INTERNAL_HPP

#include "x_function.hpp"
#include <map>
#include <stack>
#include "fx_internal.hpp"

namespace CasADi{

/** \brief  Internal node class for XFunction
    \author Joel Andersson 
    \date 2011
*/
class XFunctionInternal : public FXInternal{
  friend class XFunction;
  public:
    
    /** \brief  Constructor  */
    XFunctionInternal();
    
    /** \brief  Destructor */
    virtual ~XFunctionInternal();
  
    /** \brief  Topological sorting of the nodes based on Depth-First Search (DFS) */
    template<typename Node>
    static void sort_depth_first(std::stack<Node*>& s, std::vector<Node*>& algnodes);

    /** \brief  Topological (re)sorting of the nodes based on Bredth-First Search (BFS) (Kahn 1962) */
    template<typename Node>
    static void resort_bredth_first(std::vector<Node*>& algnodes);

    /** \brief  Topological (re)sorting of the nodes with the purpose of postponing every calculation as much as possible, as long as it does not influence a dependent node */
    template<typename Node>
    static void resort_postpone(std::vector<Node*>& algnodes, std::vector<int>& lind);
  
};

// Template implementations

template<typename Node>
void XFunctionInternal::sort_depth_first(std::stack<Node*>& s, std::vector<Node*>& algnodes){

    while(!s.empty()){
      // Get the topmost element
      Node* t = s.top();
      
      // If the last element on the stack has not yet been added
      if (!t->temp){
        
        if(t->hasDep()){
          // If a dependent node was added to the stack
          bool added_dep = false;
          
          // Add dependent nodes if not already added
          for(int i=0; i<t->ndep(); ++i){
            if(t->dep(i)->hasDep() && static_cast<Node*>(t->dep(i).get())->temp == 0) {
              // if the first child has not yet been added
              s.push(static_cast<Node*>(t->dep(i).get()));
              added_dep=true;
              break;
            }
          }
          
          if(!added_dep){    // if both children have already been added
            // Add to algorithm
            algnodes.push_back(t);

            // Mark the node as found
            t->temp = 1;

            // Remove from stack
            s.pop();
            
          }
        } else {  // If the element is a constant or symbolic node

          // Remove from stack
          s.pop();
        }
      } else {
        // If the last element on the stack has already been added
        s.pop();
      }
    }

  // Reset the node counters
  for(typename std::vector<Node*>::iterator it=algnodes.begin(); it!=algnodes.end(); ++it){
    (**it).temp = 0;
  }
}

template<typename Node>
void XFunctionInternal::resort_postpone(std::vector<Node*>& algnodes, std::vector<int>& lind){

  // Number of levels
  int nlevels = lind.size()-1;

  // Set the counter to be the corresponding place in the algorithm
  for(int i=0; i<algnodes.size(); ++i)
    algnodes[i]->temp = i;

  // Save the level of each element
  std::vector<int> level(algnodes.size());
  for(int i=0; i<nlevels; ++i)
    for(int j=lind[i]; j<lind[i+1]; ++j)
      level[j] = i;

  // Count the number of times each node is referenced inside the algorithm
  std::vector<int> numref(algnodes.size(),0);
  for(int i=0; i<algnodes.size(); ++i){
    for(int c=0; c<algnodes[i]->ndep(); ++c){ // for both children
      Node* child = static_cast<Node*>(algnodes[i]->dep(c).get());
      if(child->hasDep())
        numref[child->temp]++;
    }
  }

  // Stacks of additional nodes at the current and previous level
  std::stack<int> extra[2];

  // Loop over the levels in reverse order
  for(int i=nlevels-1; i>=0; --i){

    // The stack for the current level (we are removing elements from this stack)
    std::stack<int>& extra_this = extra[i%2]; // i odd -> use extra[1]

    // The stack for the previous level (we are adding elements to this stack)
    std::stack<int>& extra_prev = extra[1-i%2]; // i odd -> use extra[0]

    // Loop over the nodes of the level
    for(int j=lind[i]; j<lind[i+1]; ++j){
      // element to be treated
      int el = j;

      // elements in the stack have priority
      if(!extra_this.empty()){
        // Replace the element with one from the stack
        el = extra_this.top();
        extra_this.pop();
        --j; // redo the loop
      }

      // Skip the element if belongs to a higher level (i.e. was already treated)
      if(level[el] > i) continue;

      // for both children
      for(int c=0; c<algnodes[el]->ndep(); ++c){

        Node* child = static_cast<Node*>(algnodes[el]->dep(c).get());

        if(child->hasDep()){
          // Decrease the reference count of the children
          numref[child->temp]--;

          // If this was the last time the child was referenced ...
          // ... and it is not the previous level...
          if(numref[child->temp]==0 && level[child->temp] != i-1){

            // ... then assign a new level ...
            level[child->temp] = i-1;

            // ... and add to stack
            extra_prev.push(child->temp);

          } // if no more references
        } // if binary
      } // for c = ...
    } // for j
  } // for i

  // Count the number of elements on each level
  for(std::vector<int>::iterator it=lind.begin(); it!=lind.end(); ++it)
    *it = 0;
  for(std::vector<int>::const_iterator it=level.begin(); it!=level.end(); ++it)
    lind[*it + 1]++;

  // Cumsum to get the index corresponding to the first element of each level
  for(int i=0; i<nlevels; ++i)
    lind[i+1] += lind[i];

  // New index for each element
  std::vector<int> runind = lind; // running index for each level
  std::vector<int> newind(algnodes.size());
  for(int i=0; i<algnodes.size(); ++i)
    newind[i] = runind[level[algnodes[i]->temp]]++;

  // Resort the algorithm and reset the temporary
  std::vector<Node*> oldalgnodes = algnodes;
  for(int i=0; i<algnodes.size(); ++i){
    algnodes[newind[i]] = oldalgnodes[i];
    oldalgnodes[i]->temp = 0;
  }

}

template<typename Node>
void XFunctionInternal::resort_bredth_first(std::vector<Node*>& algnodes){

  // We shall assign a "level" to each element of the algorithm. A node which does not depend on other binary nodes are assigned level 0 and for nodes that depend on other nodes of the algorithm, the level will be the maximum level of any of the children plus 1. Note that all nodes of a level can be evaluated in parallel. The level will be saved in the temporary variable

  // Total number of levels
  int nlevels = 0;  

  // Get the earliest posible level
  for(typename std::vector<Node*>::iterator it=algnodes.begin(); it!=algnodes.end(); ++it){
    // maximum level of any of the children
    int maxlevel = -1;
    for(int c=0; c<(*it)->ndep(); ++c){    // Loop over the children
      Node* child = static_cast<Node*>((*it)->dep(c).get());
      if(child->hasDep() && child->temp > maxlevel)
        maxlevel = child->temp;
    }

    // Save the level of this element
    (*it)->temp = 1 + maxlevel;

    // Save if new maximum reached
    if(1 + maxlevel > nlevels)
      nlevels = 1 + maxlevel;
  }
  nlevels++;

  // Index of the first node on each level
  std::vector<int> lind;

  // Count the number of elements on each level
  lind.resize(nlevels+1,0); // all zeros to start with
  for(int i=0; i<algnodes.size(); ++i)
    lind[algnodes[i]->temp+1]++;

  // Cumsum to get the index of the first node on each level
  for(int i=0; i<nlevels; ++i)
    lind[i+1] += lind[i];

  // Get a new index for each element of the algorithm
  std::vector<int> runind = lind; // running index for each level
  std::vector<int> newind(algnodes.size());
  for(int i=0; i<algnodes.size(); ++i)
    newind[i] = runind[algnodes[i]->temp]++;

  // Resort the algorithm accordingly and reset the temporary
  std::vector<Node*> oldalgnodes = algnodes;
  for(int i=0; i<algnodes.size(); ++i){
    algnodes[newind[i]] = oldalgnodes[i];
    oldalgnodes[i]->temp = 0;
  }

#if 0

 int maxl=-1;
  for(int i=0; i<lind.size()-1; ++i){
    int l = (lind[i+1] - lind[i]);
//if(l>10)    cout << "#level " << i << ": " << l << endl;
  cout << l << ",";
    if(l>maxl) maxl= l;
  }
    cout << endl << "maxl = " << maxl << endl;

  for(int i=0; i<algnodes.size(); ++i){
    algnodes[i]->temp = i;
  }


  maxl=-1;
  for(int i=0; i<lind.size()-1; ++i){
    int l = (lind[i+1] - lind[i]);
    cout << endl << "#level " << i << ": " << l << endl;

int ii = 0;

    for(int j=lind[i]; j<lind[i+1]; ++j){

  std::vector<Node*>::const_iterator it = algnodes.begin() + j;

cout << "  "<< ii++ << ": ";

    int op = (*it)->op;
    stringstream s,s0,s1;
    s << "i_" << (*it)->temp;

    int i0 = (*it)->child[0].get()->temp;
    int i1 = (*it)->child[1].get()->temp;

    if((*it)->child[0]->hasDep())  s0 << "i_" << i0;
    else                             s0 << (*it)->child[0];
    if((*it)->child[1]->hasDep())  s1 << "i_" << i1;
    else                             s1 << (*it)->child[1];

    cout << s.str() << " = ";
    print_c[op](cout,s0.str(),s1.str());
    cout << ";" << endl;




    }

  cout << l << ",";
    if(l>maxl) maxl= l;
  }
    cout << endl << "maxl (before) = " << maxl << endl;


  for(int i=0; i<algnodes.size(); ++i){
    algnodes[i]->temp = 0;
  }


#endif

  // Resort in order to postpone all calculations as much as possible, thus saving cache
 resort_postpone(algnodes,lind);


#if 0

  for(int i=0; i<algnodes.size(); ++i){
    algnodes[i]->temp = i;
  }



  maxl=-1;
  for(int i=0; i<lind.size()-1; ++i){
    int l = (lind[i+1] - lind[i]);
    cout << endl << "#level " << i << ": " << l << endl;

int ii = 0;

    for(int j=lind[i]; j<lind[i+1]; ++j){

  std::vector<Node*>::const_iterator it = algnodes.begin() + j;

cout << "  "<< ii++ << ": ";

    int op = (*it)->op;
    stringstream s,s0,s1;
    s << "i_" << (*it)->temp;

    int i0 = (*it)->child[0].get()->temp;
    int i1 = (*it)->child[1].get()->temp;

    if((*it)->child[0]->hasDep())  s0 << "i_" << i0;
    else                             s0 << (*it)->child[0];
    if((*it)->child[1]->hasDep())  s1 << "i_" << i1;
    else                             s1 << (*it)->child[1];

    cout << s.str() << " = ";
    print_c[op](cout,s0.str(),s1.str());
    cout << ";" << endl;




    }

  cout << l << ",";
    if(l>maxl) maxl= l;
  }
    cout << endl << "maxl = " << maxl << endl;


//  return;




  for(int i=0; i<algnodes.size(); ++i){
    algnodes[i]->temp = 0;
  }



/*assert(0);*/
#endif

}

} // namespace CasADi

#endif // X_FUNCTION_INTERNAL_HPP
