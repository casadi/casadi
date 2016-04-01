//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

// /**
// \file Option.cpp
// \brief Implement the Option class
// \author Ashutosh Mahajan, Argonne National Laboratory
// */

#include <iostream>
#include <iomanip>
#include <ios>
#include <sstream>

#include "MinotaurConfig.h"
#include "Operations.h"
#include "Option.h"

using namespace Minotaur;

template <class T> Option<T>::Option(const std::string &name, 
    const std::string &desc, bool is_known, T val)
  : everUsed_(false),
    isKnown_(is_known)
{
  name_ = std::string(name);
  desc_ = std::string(desc);
  val_ = val;
  toLowerCase(name_);
}


template <class T> Option<T>::~Option()
{
}


template <class T> const std::string & Option<T>::getName()
{
  return name_;
}


template <class T> const std::string & Option<T>::getDesc()
{
  return desc_;
}


template <class T> void Option<T>::write(std::ostream &out) const
{
  std::stringstream name_stream;
  name_stream << name_ << " " << val_;
  out.flags(std::ios::left);
  out.width(30);
  out << name_stream.str();

  //out.width(4);
  //out << isKnown_;

  out << " # " << desc_;
}


// leave these lines here. They should come after the definition of all
// Option.functions(). See, for example, 
// http://www.parashift.com/c++-faq-lite/templates.html#faq-35.15
template class Option<int>;
template class Option<bool>;
template class Option<double>;
template class Option<std::string>;


//------------------------------------------------------------------------ // 
//------------------------------------------------------------------------ // 

OptionDB::OptionDB()
{
}


OptionDB::~OptionDB()
{
  bool_ops_.clear();
  int_ops_.clear();
  double_ops_.clear();
  string_ops_.clear();
}


void OptionDB::insert(BoolOptionPtr option, bool is_flag)
{
  if (is_flag) {
    flag_ops_.insert(option);
  } else {
    bool_ops_.insert(option);
  }
}


void OptionDB::insert(IntOptionPtr option)
{
  int_ops_.insert(option);
}


void OptionDB::insert(DoubleOptionPtr option)
{
  double_ops_.insert(option);
}


void OptionDB::insert(StringOptionPtr option)
{
  string_ops_.insert(option);
}


BoolOptionPtr OptionDB::findBool(const std::string &cname)
{
  BoolOptionPtr option = BoolOptionPtr(); //NULL
  std::string name(cname);
  toLowerCase(name);
  for (BoolOptionSetIter iter=bool_ops_.begin(); iter!=bool_ops_.end(); 
      ++iter) {
    if ((*iter)->getName()==name) {
      option = (*iter);
      break;
    }
  }
  return option;
}


IntOptionPtr OptionDB::findInt(const std::string &cname)
{
  IntOptionPtr option = IntOptionPtr(); //NULL
  std::string name(cname);
  toLowerCase(name);
  for (IntOptionSetIter iter=int_ops_.begin(); iter!=int_ops_.end(); 
      ++iter) {
    if ((*iter)->getName()==name) {
      option = (*iter);
      break;
    }
  }
  return option;
}


DoubleOptionPtr OptionDB::findDouble(const std::string &cname)
{
  DoubleOptionPtr option = DoubleOptionPtr(); //NULL
  std::string name(cname);
  toLowerCase(name);
  for (DoubleOptionSetIter iter=double_ops_.begin(); iter!=double_ops_.end(); 
      ++iter) {
    if ((*iter)->getName()==name) {
      option = (*iter);
      break;
    }
  }
  return option;
}


StringOptionPtr OptionDB::findString(const std::string &cname)
{
  StringOptionPtr option = StringOptionPtr(); //NULL
  std::string name(cname);
  toLowerCase(name);
  for (StringOptionSetIter iter=string_ops_.begin(); iter!=string_ops_.end(); 
      ++iter) {
    if ((*iter)->getName()==name) {
      option = (*iter);
      break;
    }
  }
  return option;
}


FlagOptionPtr OptionDB::findFlag(const std::string &cname)
{
  FlagOptionPtr option = FlagOptionPtr(); //NULL
  std::string name(cname);
  toLowerCase(name);
  for (FlagOptionSetIter iter=flag_ops_.begin(); iter!=flag_ops_.end(); 
      ++iter) {
    if ((*iter)->getName()==name) {
      option = (*iter);
      break;
    }
  }
  return option;
}


BoolOptionSetIter OptionDB::boolBegin()
{
  return bool_ops_.begin();
}

BoolOptionSetIter OptionDB::boolEnd()
{
  return bool_ops_.end();
}

IntOptionSetIter OptionDB::intBegin()
{
  return int_ops_.begin();
}

IntOptionSetIter OptionDB::intEnd()
{
  return int_ops_.end();
}

DoubleOptionSetIter OptionDB::dblBegin()
{
  return double_ops_.begin();
}

DoubleOptionSetIter OptionDB::dblEnd()
{
  return double_ops_.end();
}

StringOptionSetIter OptionDB::strBegin()
{
  return string_ops_.begin();
}

StringOptionSetIter OptionDB::strEnd()
{
  return string_ops_.end();
}

FlagOptionSetIter OptionDB::flagBegin()
{
  return flag_ops_.begin();
}

FlagOptionSetIter OptionDB::flagEnd()
{
  return flag_ops_.end();
}


void OptionDB::write(std::ostream &out) const
{
  out << "## boolean options:" << std::endl;
  //out.flush();
  for (BoolOptionSetIter iter=bool_ops_.begin(); iter!=bool_ops_.end(); 
      ++iter) {
    (*iter)->write(out);
    out << std::endl;
  }
  out << std::endl;

  out << "## int options:" << std::endl;
  //out.flush();
  for (IntOptionSetIter iter=int_ops_.begin(); iter!=int_ops_.end(); 
      ++iter) {
    (*iter)->write(out);
    out << std::endl;
  }
  out << std::endl;

  out << "## double options:" << std::endl;
  //out.flush();
  for (DoubleOptionSetIter iter=double_ops_.begin(); iter!=double_ops_.end(); 
      ++iter) {
    (*iter)->write(out);
    out << std::endl;
  }
  out << std::endl;

  out << "## String options:" << std::endl;
  //out.flush();
  for (StringOptionSetIter iter=string_ops_.begin(); iter!=string_ops_.end(); 
      ++iter) {
    (*iter)->write(out);
    out << std::endl;
  }
  out << std::endl;

  out << "## Flags:" << std::endl;
  for (FlagOptionSetIter iter=flag_ops_.begin(); iter!=flag_ops_.end(); 
      ++iter) {
    (*iter)->write(out);
    out << std::endl;
  }
  out << std::endl;
}


// Local Variables: 
// mode: c++ 
// eval: (c-set-style "k&r") 
// eval: (c-set-offset 'innamespace 0) 
// eval: (setq c-basic-offset 2) 
// eval: (setq fill-column 78) 
// eval: (auto-fill-mode 1) 
// eval: (setq column-number-mode 1) 
// eval: (setq indent-tabs-mode nil) 
// End:
