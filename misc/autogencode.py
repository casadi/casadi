#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             K.U. Leuven. All rights reserved.
#     Copyright (C) 2011-2014 Greg Horn
#
#     CasADi is free software; you can redistribute it and/or
#     modify it under the terms of the GNU Lesser General Public
#     License as published by the Free Software Foundation; either
#     version 3 of the License, or (at your option) any later version.
#
#     CasADi is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     Lesser General Public License for more details.
#
#     You should have received a copy of the GNU Lesser General Public
#     License along with CasADi; if not, write to the Free Software
#     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
#
"""

This script generates parts of CasADi source code that are highly repetitive and easily automated

"""
import os, fnmatch
import re

# Walk into directories in filesystem
# Ripped from os module and slightly modified
# for alphabetical sorting
#
def sortedWalk(top, topdown=True, onerror=None):
  from os.path import join, isdir, islink

  names = os.listdir(top)
  names.sort()
  dirs, nondirs = [], []
  for name in names:
    if isdir(os.path.join(top, name)):
      dirs.append(name)
    else:
      nondirs.append(name)
  if topdown:
    yield top, dirs, nondirs
  for name in dirs:
    path = join(top, name)
    if not os.path.islink(path):
      for x in sortedWalk(path, topdown, onerror):
        yield x
  if not topdown:
    yield top, dirs, nondirs

def locate(pattern, root=os.curdir):
    """
    Locate all files matching supplied filename pattern in and below
    supplied root directory.
    """
    for path, dirs, files in sortedWalk(os.path.abspath(root)):
        dirs.sort()

        for filename in fnmatch.filter(files, pattern):
            yield os.path.join(path, filename)


# Part 1:    helper functions for inputs/outputs


# Some regexes
re_introIn = re.compile("\s*/// (Input arguments .*) \[(\w+)\]")
re_introOut = re.compile("\s*/// (Output arguments .*) \[(\w+)\]")
re_introStruct = re.compile("\s*/// (Structure specification of .*) \[(\w+)\]")
re_doc = re.compile("\s*///\s*(.*) \[([\w\d]+)\]\s*")
re_docraw = re.compile("\s*///\s*(.*)")
re_enumheader = re.compile("\s*enum (\w+)\s*{")
re_enum = re.compile("\s*(\w+),\s*")
re_end = re.compile("\s*}\s*;")

class Enum:
  def __init__(self,name, doc):
    self.name    = name
    self.doc     = doc + "\n"
    self.entries = []

  def addArg(self, name, doc):
    self.entries.append([name,doc])

  def getDoc(self):
    return self.doc.replace("\\n","")

  def addDoc(self,doc):
    self.doc+=doc + "\n"

  def addEnum(self, enum):
    if len(self.entries[-1])>2:
      raise Exception("Enum %s already given" % enum)
    self.entries[-1].append(enum)

  def addEnumheader(self,header):
    self.enum = header

  def __str__(self):
    s = self.name + "(" + self.doc + ")" + "[" + self.enum  + "]" + ": " + str(self.entries)
    return s

  def checkconsistency(self):
    for e in self.entries:
      if not(len(e))==3:
         raise Exception("Consistency"+str(e))

    print self.enum, self.__class__.__name__
    assert(self.enum.endswith(self.__class__.__name__))
    prefix = self.enum[:-len(self.__class__.__name__)]

    for name, doc, enum in self.entries:
      assert(enum.startswith(enum))
      assert(not("_NUM" in enum))

  def num(self):
    #return self.enum[:-len(self.__class__.__name__)] +  "_NUM_" + self.PRE
    return str(len(self.entries))

  def hppcode(self):
    s= "/// \cond INTERNAL\n/// Helper function for '" + self.enum + "'\n"
    s+="""
template<class M>
class CASADI_CORE_EXPORT %sIOSchemeVector : public IOSchemeVector<M> {
  public:
    explicit %sIOSchemeVector(const std::vector<M>& t)
      : IOSchemeVector<M>(t, SCHEME_%s) {}
};
/// \endcond
""" % (self.enum,self.enum,self.enum)
    #if self.enum.endswith("Struct"):
    #  s+="typedef %sIOSchemeVector<Sparsity> %sure;\n" % (self.enum,self.enum)
    s+= "\n".join(map(lambda x: ("/// " + x).rstrip(),self.doc.split("\n")))+"\n"
    s+= "/// \\copydoc scheme_" + self.enum +  "\n"
    s+= "template<class M>" + "\n"
    s+= self.enum + "IOSchemeVector<M> " + self.name + "("
    for i, (name, doc, enum) in enumerate(self.entries):
      s+="""\n    const std::string &arg_s%d ="", const M &arg_m%d =M()""" % (i,i) + ","
    s=s[:-1] + ") {\n"
    s+= "  std::vector<M> ret(%d);\n" % len(self.entries)


    s+= "  std::map<std::string, M> arg;\n"
    for i,_ in enumerate(self.entries):
      s+="""  if (arg_s%d != "") arg.insert(make_pair(arg_s%d, arg_m%d));\n""" % (i,i,i)
    s+="""  typedef typename std::map<std::string, M>::const_iterator it_type;\n"""
    s+="""  for (it_type it = arg.begin(); it != arg.end(); it++) {
    int n = getSchemeEntryEnum(SCHEME_%s, it->first);
    if (n==-1)
      casadi_error("Keyword error in %s: '" << it->first
        << "' is not recognized. Available keywords are: "
        "%s");  // NOLINT(whitespace/line_length)
    ret[n] = it->second;
  }
"""  % (self.enum,self.enum,", ".join([name for name, doc, enum in self.entries]))
    s+= "  return %sIOSchemeVector<M>(ret);" % self.enum
    s+="\n}\n"
    s+= "template<class M>" + "\n"
    s+= "std::vector<M> " + self.name + "(const std::vector<M>& args,"
    for i, (name, doc, enum) in enumerate(self.entries):
      s+='\n    const std::string &arg_s%d=""' % i + ","
    s=s[:-1] + ") {\n"
    s+= "  std::vector<M> ret;\n"
    for i,_ in enumerate(self.entries):
      s+="""  if (arg_s%d != "") ret.push_back(args.at(getSchemeEntryEnum(SCHEME_%s, arg_s%d))); // NOLINT(whitespace/line_length)\n""" % (i,self.enum,i)
    s+= "  return ret;\n"
    s+="\n}\n"

    return s
    
  def cppcode(self):

    return ""

  def swigcode(self):
    s="namespace casadi {\n"
    #s+= "%warnfilter(302) " + self.name+ ";\n" -- does not seem to work
    if self.enum.endswith('VecStruct'):
      s+="%template(" + self.name + ") " + self.name + "< std::vector<casadi::Sparsity> >;\n"
    elif self.enum.endswith('Struct'):
      s+="%template(" + self.name + ") " + self.name + "<casadi::Sparsity>;\n"
    else:
      s+="%template(" + self.name + ") " + self.name + "<casadi::SX>;\n"
      s+="%template(" + self.name + ") " + self.name + "<casadi::MX>;\n"
      s+="%template(" + self.name + ") " + self.name + "<casadi::Matrix<double> >;\n"
      s+="%template(" + self.name + ") " + self.name + "<casadi::Sparsity>;\n"
      s+="%template(" +  "IOSchemeVector" + self.enum + "SX) " + self.enum + "IOSchemeVector<SX>;\n"
      s+="%template(" +  "IOSchemeVector" + self.enum + "MX) " + self.enum + "IOSchemeVector<MX>;\n"
      s+="%template(" +  "IOSchemeVector" + self.enum + "D) " + self.enum + "IOSchemeVector< Matrix<double> >;\n"
      s+="%template(" +  "IOSchemeVector" + self.enum + "Sparsity) " + self.enum + "IOSchemeVector<Sparsity>;\n"
      #s+="%rename(" + self.name + ") " + self.name + "<casadi::SX>;\n"
      #s+="%rename(" + self.name + ") " + self.name + "<casadi::MX>;\n"
      #s+="%rename(" + self.name + ") " + self.name + "<casadi::Sparsity>;\n"
      s+="%rename(" + "IOSchemeVector" + self.enum + ") " + "IOSchemeVector" + self.enum + "SX;\n"
      s+="%rename(" + "IOSchemeVector" + self.enum + ") " + "IOSchemeVector" + self.enum + "MX;\n"
      s+="%rename(" + "IOSchemeVector" + self.enum + ") " + "IOSchemeVector" + self.enum + "D;\n"
      s+="%rename(" + "IOSchemeVector" + self.enum + ") " + "IOSchemeVector" + self.enum+ "Sparsity;\n"
    s+="}\n"
    return s

  def pureswigcode(self):
    s="namespace casadi {\n"
    if self.enum.endswith('VecStruct'):
      s+="%template(" + self.enum + "ure) " + self.enum + "IOSchemeVector< std::vector<Sparsity> >;\n"
    elif self.enum.endswith('Struct'):
      s+="%template(" + self.enum + "ure) " + self.enum + "IOSchemeVector<Sparsity>;\n"
    s+="}\n"
    return s

  def pycode(self):
    s="def " + self.name + "(*dummy,**kwargs):\n"
    s+='  """\n'
    s+= "  Helper function for '" + self.enum + "'\n\n"
    s+= "  Two use cases:\n"
    s+= "     a) arg = %s(%s)\n" % (self.name , ", ".join(["%s=my_%s" % (name,name) for name, doc, enum in self.entries]))
    s+= "          all arguments optional\n"
    s+= "     b) %s = %s(arg,%s)\n" % ( ", ".join([name for name, doc, enum in self.entries]), self.name , ", ".join(['"' + name+'"' for name, doc, enum in self.entries]))
    s+= "          all arguments after the first optional\n"
    s+= "\n".join(map(lambda x: "  " + x.rstrip(),self.getDoc().split("\n"))) + "\n"
    s+= "  Keyword arguments::\n\n"
    maxlenname = max([len(name) for name, doc, enum in self.entries])
    for name, doc, enum in self.entries:
      s+="    " + name + (" "*(maxlenname-len(name))) +  " -- " +  doc + " [" + enum + "]\n"
    s+='  """\n'
    s+="""  if (len(dummy)>0 and len(kwargs)>0): raise Exception("Cannot mix two use cases of %s. Either use keywords or non-keywords ")\n""" % self.name
    s+="""  if len(dummy)>0: return [ dummy[0][getSchemeEntryEnum(SCHEME_%s,n)] for n in dummy[1:]]\n""" % self.enum
    for name, doc, enum in self.entries:
      s+="  %s = %s\n  if '%s' in kwargs:\n    %s = kwargs['%s']\n" % (name,"Sparsity()" if self.enum.endswith("Struct") and not self.enum.endswith("VecStruct") else "[]",name,name,name)
    s+="""  for k in kwargs.keys():\n    if not(k in [%s]):\n      raise Exception("Keyword error in %s: '%%s' is not recognized. Available keywords are: %s" %% k )\n""" % (",".join(["'%s'" % name for name, doc, enum in self.entries]),self.name,", ".join([name for name, doc, enum in self.entries]))

    if (self.enum.endswith("Struct")):
      s+="  return %sure([" % self.enum
      for name, doc, enum in self.entries:
        s+=name+","
      s=s[:-1] + "])\n"
    else:
      s+="  return IOSchemeVector(["
      for name, doc, enum in self.entries:
        s+=name+","
      s=s[:-1] + "], IOScheme(SCHEME_%s))\n" % self.enum
    return s

class Input(Enum):
  PRE = "IN"
  tp = "Input"
  pass

class Output(Enum):
  PRE = "OUT"
  tp = "Output"
  pass

class Struct(Enum):
  PRE = "STRUCT"
  tp = "Structure"
  pass

# This has a write method which keepds adding lines, but on .close()
# it will only write out to a file if the file doesn't exist, or if there
# are changes to the file.
# This prevents a rebuild from happening if nothing changed
class LazyFile(object):
  def __init__(self, path):
    self._path = path
    self._lines = []

  def write(self, txt):
    self._lines.append(txt)

  def close(self):
    new_content = ''.join(self._lines)

    if os.path.isfile(self._path):
      with open(self._path,'r') as content_file:
        old_content = content_file.read()
      if old_content == new_content:
        print '(unchanged) "' + self._path + '"'
      else:
        print '(changed)   "' + self._path + '"'
        self._reallyWrite(new_content)

    else:
      print '(new)       "' + self._path + '"'
      self._reallyWrite(new_content)


  def _reallyWrite(self,contents):
    f = file(self._path,'w')
    f.write(''.join(self._lines))
    f.close()
    
  def __del__(self):
    self.close()

autogenmetadatahpp = LazyFile(os.path.join(os.curdir,"..","casadi","core","function","schemes_metadata.hpp"))
autogenhelpershpp = LazyFile(os.path.join(os.curdir,"..","casadi","core","function","schemes_helpers.hpp"))
autogenpy = LazyFile(os.path.join(os.curdir,"..","swig","autogenerated.i"))
autogencpp = LazyFile(os.path.join(os.curdir,"..","casadi","core","function","schemes_metadata.cpp"))


autogenmetadatahpp.write(file('license_header.txt','r').read())
autogenmetadatahpp.write("/** All edits to this file will be lost - autogenerated by misc/autogencode.py */\n")
autogenmetadatahpp.write("""#ifndef SCHEMES_METADATA_HPP\n#define SCHEMES_METADATA_HPP\n#include <vector>\n#include <string>\n#include <utility>\n#include <map>\n#include "../casadi_exception.hpp"\nnamespace casadi {\ntemplate <class T>
class IOSchemeVector;class Sparsity;\n""")
autogenhelpershpp.write(file('license_header.txt','r').read())
autogenhelpershpp.write("/** All edits to this file will be lost - autogenerated by misc/autogencode.py */\n")
autogenhelpershpp.write("""#ifndef SCHEMES_HELPERS_HPP\n#define SCHEMES_HELPERS_HPP\n#include <vector>\n#include <string>\n#include <utility>\n#include <map>\n#include "io_scheme_vector.hpp"\nnamespace casadi {\n\n""")
autogencpp.write(file('license_header.txt','r').read())
autogencpp.write("/** All edits to this file will be lost - autogenerated by misc/autogencode.py */\n")
autogencpp.write("""#include "schemes_metadata.hpp"\n#include <string>\nnamespace casadi {\n""")

autogenpy.write(file('license_header.txt','r').read())
autogenpy.write("/** All edits to this file will be lost - autogenerated by misc/autogencode.py */\n")
autogenpy.write("#ifndef AUTOGENERATED_I\n")
autogenpy.write("#define AUTOGENERATED_I\n")
autogenpy.write("""%include "casadi/core/function/schemes_metadata.hpp"\n%include "casadi/core/function/schemes_helpers.hpp"\n""")

schemes = []

for h in locate("*.hpp",os.path.join(os.curdir,"../casadi/core")):
  f =  file(h,'r')
  while 1:
    line = f.readline()
    if not(line):
      break
    p = None
    m = re.match(re_introIn,line)
    if m:
      p = Input(m.group(2),m.group(1))
    m = re.match(re_introOut,line)
    if m:
      p = Output(m.group(2),m.group(1))
    m = re.match(re_introStruct,line)
    if m:
      p = Struct(m.group(2),m.group(1))
    if p:
      line = f.readline()
      while re.search(re_docraw, line):
        p.addDoc(re.search(re_docraw,line).group(1))
        line = f.readline()

      while not(re.search(re_end, line)):
        mm = re.search(re_doc, line )
        if mm:
          p.addArg(mm.group(2),mm.group(1))
        mm = re.match(re_enum, line )
        if mm:
          p.addEnum(mm.group(1))
        mm = re.match(re_enumheader, line )
        if mm:
          p.addEnumheader(mm.group(1))
        line = f.readline()
      p.checkconsistency()
      schemes.append(p)

autogenmetadatahpp.write("enum InputOutputScheme {\n  %s };\n\n" % ",\n  ".join(["SCHEME_"+p.enum for p in schemes]) )

autogenmetadatahpp.write("CASADI_CORE_EXPORT std::string getSchemeEntryName(InputOutputScheme scheme, int i);\n")
autogenmetadatahpp.write("CASADI_CORE_EXPORT std::string getSchemeEntryDoc(InputOutputScheme scheme, int i);\n")
autogenmetadatahpp.write("CASADI_CORE_EXPORT std::string getSchemeEntryEnumName(InputOutputScheme scheme, int i);\n")
autogenmetadatahpp.write("CASADI_CORE_EXPORT int getSchemeEntryEnum(InputOutputScheme scheme, const std::string &name);\n")
autogenmetadatahpp.write("CASADI_CORE_EXPORT int getSchemeSize(InputOutputScheme scheme);\n")
autogenmetadatahpp.write("CASADI_CORE_EXPORT std::string getSchemeName(InputOutputScheme scheme);\n")
autogenmetadatahpp.write("CASADI_CORE_EXPORT std::string getSchemeEntryNames(InputOutputScheme scheme);\n")

autogenpy.write("#ifdef SWIGPYTHON\n")
autogenpy.write("%pythoncode %{\n")
autogenpy.write("""
def IOSchemeVector(arg,io_scheme):
  try:
    return IOSchemeVectorD(arg,io_scheme)
  except:
    pass
  try:
    return IOSchemeVectorSX(arg,io_scheme)
  except:
    pass
  try:
    return IOSchemeVectorMX(arg,io_scheme)
  except:
    pass
  try:
    arg = map(lambda x: sp_dense(0,0) if isinstance(x,list) and len(x)==0 else x,arg)
    return IOSchemeVectorSparsity(arg,io_scheme)
  except:
    pass
  raise TypeError("IOSchemeVector called with faulty arguments. Individual values must be SX, MX or Sparsity.")
""")
autogenpy.write("%}\n")
autogenpy.write("#endif //SWIGPYTHON\n")
for p in schemes:
  print p.name
  autogencpp.write(p.cppcode())
  autogenhelpershpp.write(p.hppcode())
  autogenpy.write("#ifdef SWIGPYTHON\n")
  autogenpy.write("%pythoncode %{\n")
  autogenpy.write(p.pycode())
  autogenpy.write("%}\n")
  autogenpy.write("#endif //SWIGPYTHON\n")
  autogenpy.write("#ifndef SWIGPYTHON\n")
  autogenpy.write(p.swigcode())
  autogenpy.write("#endif //SWIGPYTHON\n")
  autogenpy.write(p.pureswigcode())
  
autogenhelpershpp.write("#define INSTANTIATE_IOSCHEME_HELPERS(T) \\\n")
for p in schemes:
  autogenhelpershpp.write("template class %sIOSchemeVector<T>;\\\n" % p.enum)
autogenhelpershpp.write("\n")

autogencpp.write("std::string getSchemeName(InputOutputScheme scheme) {\n  switch (scheme) {\n")
for p in schemes:
  autogencpp.write("""    case SCHEME_%s: return "%s";\n""" % (p.enum,p.enum))
autogencpp.write("""  default: casadi_error("getSchemeName: Scheme '" << scheme <<  "' does not exist.");\n""")
autogencpp.write("  }\n}\n")

autogencpp.write("std::string getSchemeEntryNames(InputOutputScheme scheme) {\n  switch (scheme) {\n")
for p in schemes:
  autogencpp.write("""    case SCHEME_%s:\n      return "%s";\n""" % (p.enum,", ".join([name for name, doc, enum in p.entries])))
autogencpp.write("""  default: casadi_error("getSchemeName: Scheme '" << scheme <<  "' does not exist.");\n""")
autogencpp.write("  }\n}\n")

autogencpp.write("std::string getSchemeEntryName(InputOutputScheme scheme, int i) {\n  switch (scheme) {\n")
for p in schemes:
  autogencpp.write("    case SCHEME_%s:\n" % p.enum)
  for i, (name, doc, enum) in enumerate(p.entries):
    autogencpp.write("""      if (i==%d) return "%s";\n""" % (i,name))
  autogencpp.write("      break;\n")
autogencpp.write("  }\n")
autogencpp.write("""\
  casadi_error("getSchemeEntryName: supplied number is out of range. Scheme '"
               << getSchemeName(scheme) << "' has only " << getSchemeSize(scheme)
               << " entries: " << getSchemeEntryNames(scheme) << ".");
""")
autogencpp.write("}\n")

autogencpp.write("std::string getSchemeEntryDoc(InputOutputScheme scheme, int i) {\n  switch (scheme) {\n")
for p in schemes:
  autogencpp.write("    case SCHEME_%s:\n" % p.enum)
  for i, (name, doc, enum) in enumerate(p.entries):
    autogencpp.write("""      if (i==%d) return "%s";  // NOLINT(whitespace/line_length)\n""" % (i,doc))
  autogencpp.write("      break;\n")
autogencpp.write("  }\n")
autogencpp.write("""\
  casadi_error("getSchemeEntryDoc: supplied number is out of range. Scheme '"
               << getSchemeName(scheme) << "' has only " << getSchemeSize(scheme)
               << " entries: " << getSchemeEntryNames(scheme) << ".");\n""")
autogencpp.write("}\n")

autogencpp.write("std::string getSchemeEntryEnumName(InputOutputScheme scheme, int i) {\n  switch (scheme) {\n")
for p in schemes:
  autogencpp.write("    case SCHEME_%s:\n" % p.enum)
  for i, (name, doc, enum) in enumerate(p.entries):
    autogencpp.write("""      if (i==%d) return "%s";\n""" % (i,enum))
  autogencpp.write("      break;\n")
autogencpp.write("  }\n")
autogencpp.write("""\
  casadi_error("getSchemeEntryEnumName: supplied number is out of range. Scheme '"
               << getSchemeName(scheme) << "' has only "
               << getSchemeSize(scheme) << " entries: "
                << getSchemeEntryNames(scheme) << ".");\n""")
autogencpp.write("}\n")

autogencpp.write("int getSchemeSize(InputOutputScheme scheme) {\n  switch (scheme) {\n")
for p in schemes:
  autogencpp.write("    case SCHEME_%s:\n" % p.enum)
  autogencpp.write("""      return %d;\n""" % len(p.entries))
  autogencpp.write("      break;\n")
autogencpp.write("""  default: casadi_error("getSchemeSize: Scheme '" << scheme <<  "' does not exist.");\n""")
autogencpp.write("  }\n")
autogencpp.write("}\n")

autogencpp.write("int getSchemeEntryEnum(InputOutputScheme scheme, const std::string &name) {\n  switch (scheme) {\n")
for p in schemes:
  autogencpp.write("    case SCHEME_%s:\n" % p.enum)
  for i, (name, doc, enum) in enumerate(p.entries):
    autogencpp.write("""      if (name=="%s") return %d;\n""" % (name,i))
  autogencpp.write("      break;\n")
autogencpp.write("  }\n")
autogencpp.write("""\
  casadi_error("getSchemeEntryEnum: Scheme '" << getSchemeName(scheme)
               <<  "' has no entry named '" << name
               <<  "'. Available entries are: "
                << getSchemeEntryNames(scheme) << ".");\n""")
autogencpp.write("}\n")
autogenmetadatahpp.write("} // namespace casadi\n#endif //SCHEMES_METADATA_HPP")
autogenhelpershpp.write("} // namespace casadi\n#endif //SCHEMES_HELPERS_HPP")
autogencpp.write("} // namespace casadi")

autogenpy.write("#endif //AUTOGENERATED_I\n")

# Endline at the end of files
autogenmetadatahpp.write("\n\n")
autogenhelpershpp.write("\n\n")
autogencpp.write("\n\n")
autogenpy.write("\n")
