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
    s= "/// Helper function for '" + self.enum + "'\n"
    s+="""
template<class M>
class %sIOSchemeVector : public IOSchemeVector<M> {
  public:
    explicit %sIOSchemeVector(const std::vector<M>& t) : IOSchemeVector<M>(t,SCHEME_%s){} 
};
""" % (self.enum,self.enum,self.enum)
    #if self.enum.endswith("Struct"):
    #  s+="typedef %sIOSchemeVector<CRSSparsity> %sure;\n" % (self.enum,self.enum)
    s+= "\n".join(map(lambda x: "/// " + x.rstrip(),self.doc.split("\n")))+"\n"
    s+= "/// \\copydoc scheme_" + self.enum +  "\n"
    s+= "template<class M>" + "\n"
    s+= self.enum + "IOSchemeVector<M> " + self.name + "("
    for i, (name, doc, enum) in enumerate(self.entries):
      s+="""const std::string arg_s%d="",M arg_m%d=M()""" % (i,i) + ","
    s=s[:-1] + "){" + "\n"
    s+= "  std::vector<M> ret(%d);\n" % len(self.entries)
    
    
    s+= "  std::map<std::string,M> arg;\n"
    for i,_ in enumerate(self.entries):
      s+="""  if (arg_s%d!="") arg.insert(make_pair(arg_s%d,arg_m%d));\n""" % (i,i,i)
    s+="""  typedef typename std::map<std::string,M>::const_iterator it_type;\n"""
    s+="""  for(it_type it = arg.begin(); it != arg.end(); it++) {
    int n = getSchemeEntryEnum(SCHEME_%s,it->first);
    if (n==-1)
      casadi_error("Keyword error in %s: '" << it->first << "' is not recognized. Available keywords are: %s");
    ret[n] = it->second;
  }
"""  % (self.enum,self.enum,", ".join([name for name, doc, enum in self.entries]))
    s+= "  return %sIOSchemeVector<M>(ret);" % self.enum
    s+="\n}\n"
    s+= "template<class M>" + "\n"
    s+= "std::vector<M> " + self.name + "(const std::vector<M>& args,"
    for i, (name, doc, enum) in enumerate(self.entries):
      s+='const std::string arg_s%d=""' % i + ","
    s=s[:-1] + "){" + "\n"
    s+= "  std::vector<M> ret;\n"
    for i,_ in enumerate(self.entries):
      s+="""  if (arg_s%d!="") ret.push_back(args.at(getSchemeEntryEnum(SCHEME_%s,arg_s%d)));\n""" % (i,self.enum,i)
    s+= "  return ret;\n"
    s+="\n}\n"

    return s

  def cppcode(self):

    return ""
    
  def swigcode(self):
    s="namespace CasADi {\n"
    if self.enum.endswith('Struct'):
      s+="%template(" + self.name + ") " + self.name + "<CRSSparsity>;\n"
    else:
      s+="%template(" + self.name + ") " + self.name + "<SXMatrix>;\n"
      s+="%template(" + self.name + ") " + self.name + "<MX>;\n"
      s+="%template(" + self.name + ") " + self.name + "<CRSSparsity>;\n"
      s+="%template(" +  "IOSchemeVector" + self.enum + ") " + self.enum + "IOSchemeVector<SXMatrix>;\n"
      s+="%template(" +  "IOSchemeVector" + self.enum + ") " + self.enum + "IOSchemeVector<MX>;\n"
      s+="%template(" +  "IOSchemeVector" + self.enum + ") " + self.enum + "IOSchemeVector<CRSSparsity>;\n"
    s+="}\n"
    return s

  def pureswigcode(self):
    s="namespace CasADi {\n"
    if self.enum.endswith('Struct'):
      s+="%template(" + self.enum + "ure) " + self.enum + "IOSchemeVector<CRSSparsity>;\n"
    s+="}\n"
    return s
    
  def pycode(self):
    s="def " + self.name + "(*dummy,**kwargs):\n"
    s+='  """\n'
    s+= "  Helper function for '" + self.enum + "'\n\n"
    s+= "  Two use cases:\n"
    s+= "     a) arg = %s(%s) \n" % (self.name , ", ".join(["%s=my_%s" % (name,name) for name, doc, enum in self.entries]))
    s+= "          all arguments optional\n"
    s+= "     b) %s = %s(arg,%s) \n" % ( ", ".join([name for name, doc, enum in self.entries]), self.name , ", ".join(['"' + name+'"' for name, doc, enum in self.entries]))
    s+= "          all arguments after the first optional\n"
    s+= "\n".join(map(lambda x: "  " + x.rstrip(),self.getDoc().split("\n"))) + "\n"
    s+= "  Keyword arguments:\n"
    maxlenname = max([len(name) for name, doc, enum in self.entries])
    for name, doc, enum in self.entries:
      s+="    " + name + (" "*(maxlenname-len(name))) +  " -- " +  doc + " [" + enum + "]\n"
    s+='  """\n'
    s+="""  if(len(dummy)>0 and len(kwargs)>0): raise Exception("Cannot mix two use cases of %s. Either use keywords or non-keywords ")\n""" % self.name
    s+="""  if len(dummy)>0: return [ dummy[0][getSchemeEntryEnum(SCHEME_%s,n)] for n in dummy[1:]]\n""" % self.enum
    for name, doc, enum in self.entries:
      s+="  %s = []\n  if '%s' in kwargs:\n    %s = kwargs['%s']\n" % (name,name,name,name)
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

autogenmetadatahpp = file(os.path.join(os.curdir,"..","symbolic","fx","schemes_metadata.hpp"),"w")
autogenhelpershpp = file(os.path.join(os.curdir,"..","symbolic","fx","schemes_helpers.hpp"),"w")
autogenpy = file(os.path.join(os.curdir,"..","swig","autogenerated.i"),"w")
autogencpp = file(os.path.join(os.curdir,"..","symbolic","fx","schemes_metadata.cpp"),"w")


autogenmetadatahpp.write(file('license_header.txt','r').read())
autogenmetadatahpp.write("/** All edits to this file will be lost - autogenerated by misc/autogencode.py */\n")
autogenmetadatahpp.write("""#ifndef SCHEMES_METADATA_HPP\n#define SCHEMES_METADATA_HPP\n#include <vector>\n#include <string>\n#include <utility>\n#include <map>\n#include "../casadi_exception.hpp"\nnamespace CasADi{ \ntemplate <class T>
class IOSchemeVector;class CRSSparsity;\n""")
autogenhelpershpp.write(file('license_header.txt','r').read())
autogenhelpershpp.write("/** All edits to this file will be lost - autogenerated by misc/autogencode.py */\n")
autogenhelpershpp.write("""#ifndef SCHEMES_HELPERS_HPP\n#define SCHEMES_HELPERS_HPP\n#include <vector>\n#include <string>\n#include <utility>\n#include <map>\n#include "io_scheme_vector.hpp"\nnamespace CasADi{ \n\n""")
autogencpp.write(file('license_header.txt','r').read())
autogencpp.write("/** All edits to this file will be lost - autogenerated by misc/autogencode.py */\n")
autogencpp.write("""#include "schemes_metadata.hpp"\n#include <string>\nnamespace CasADi{\n""")

autogenpy.write(file('license_header.txt','r').read())
autogenpy.write("/** All edits to this file will be lost - autogenerated by misc/autogencode.py */\n")
autogenpy.write("#ifndef AUTOGENERATED_I\n")
autogenpy.write("#define AUTOGENERATED_I\n")
autogenpy.write("""%include "symbolic/fx/schemes_metadata.hpp"\n%include "symbolic/fx/schemes_helpers.hpp"\n""")

schemes = []

for h in locate("*.hpp",os.path.join(os.curdir,"..")):
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
      
autogenmetadatahpp.write("enum InputOutputScheme { %s };\n" % ", ".join(["SCHEME_"+p.enum for p in schemes]) )

autogenmetadatahpp.write("std::string getSchemeEntryName(InputOutputScheme scheme, int i);\n")
autogenmetadatahpp.write("std::string getSchemeEntryDoc(InputOutputScheme scheme, int i);\n")
autogenmetadatahpp.write("std::string getSchemeEntryEnumName(InputOutputScheme scheme, int i);\n")
autogenmetadatahpp.write("int getSchemeEntryEnum(InputOutputScheme scheme, const std::string &name);\n")
autogenmetadatahpp.write("int getSchemeSize(InputOutputScheme scheme);\n")
autogenmetadatahpp.write("std::string getSchemeName(InputOutputScheme scheme);\n")
autogenmetadatahpp.write("std::string getSchemeEntryNames(InputOutputScheme scheme);\n")

autogenpy.write("#ifdef SWIGPYTHON\n")
autogenpy.write("%pythoncode %{\n")
autogenpy.write("""
def IOSchemeVector(arg,io_scheme):
  try:
    return IOSchemeVectorSXMatrix(arg,io_scheme)
  except:
    pass
  try:
    return IOSchemeVectorMX(arg,io_scheme)
  except:
    pass
  try:
    arg = map(lambda x: sp_dense(0,0) if isinstance(x,list) and len(x)==0 else x,arg)
    return IOSchemeVectorCRSSparsity(arg,io_scheme)
  except:
    pass
    
def customIO(**kwargs):
  items = kwargs.items()
  
  return IOSchemeVector(zip(*items)[1], IOScheme(zip(*items)[0]))
  
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

autogenhelpershpp.write("""
/// Helper function for 'customIO'
template<class M>
IOSchemeVector<M> customIO(""" + ",".join(['const std::string arg_s%d="",M arg_m%d=M()' % (i,i) for i in range(20)]) +"""){
  std::vector<std::string> k;
  std::vector<M> v;
"""+
  "\n".join(['  if (arg_s%d!="") { k.push_back(arg_s%d);  v.push_back(arg_m%d); }' % (i,i,i) for i in range(20) ])
+"""
  return IOSchemeVector<M>(v,IOScheme(k));
}
""")
  
autogencpp.write("std::string getSchemeName(InputOutputScheme scheme) {\n  switch (scheme) {\n")
for p in schemes:
  autogencpp.write("""    case SCHEME_%s: return "%s";\n""" % (p.enum,p.enum))
autogencpp.write("  }\n}\n")

autogencpp.write("std::string getSchemeEntryNames(InputOutputScheme scheme) {\n  switch (scheme) {\n")
for p in schemes:
  autogencpp.write("""    case SCHEME_%s: return "%s";\n""" % (p.enum,", ".join([name for name, doc, enum in p.entries])))
autogencpp.write("  }\n}\n")

autogencpp.write("std::string getSchemeEntryName(InputOutputScheme scheme, int i) {\n  switch (scheme) {\n")
for p in schemes:
  autogencpp.write("    case SCHEME_%s: \n" % p.enum)
  for i, (name, doc, enum) in enumerate(p.entries):
    autogencpp.write("""      if(i==%d) return "%s";\n""" % (i,name))
  autogencpp.write("      break;\n")
autogencpp.write("  }\n")
autogencpp.write("""  casadi_error("getSchemeEntryName: supplied number is out of range. Scheme '" << getSchemeName(scheme) << "' has only " << getSchemeSize(scheme) << " entries: " << getSchemeEntryNames(scheme) << ".");\n""")
autogencpp.write("}\n")

autogencpp.write("std::string getSchemeEntryDoc(InputOutputScheme scheme, int i) {\n  switch (scheme) {\n")
for p in schemes:
  autogencpp.write("    case SCHEME_%s: \n" % p.enum)
  for i, (name, doc, enum) in enumerate(p.entries):
    autogencpp.write("""      if(i==%d) return "%s";\n""" % (i,doc))
  autogencpp.write("      break;\n")
autogencpp.write("  }\n")
autogencpp.write("""  casadi_error("getSchemeEntryDoc: supplied number is out of range. Scheme '" << getSchemeName(scheme) << "' has only " << getSchemeSize(scheme) << " entries: " << getSchemeEntryNames(scheme) << ".");\n""")
autogencpp.write("}\n")

autogencpp.write("std::string getSchemeEntryEnumName(InputOutputScheme scheme, int i) {\n  switch (scheme) {\n")
for p in schemes:
  autogencpp.write("    case SCHEME_%s: \n" % p.enum)
  for i, (name, doc, enum) in enumerate(p.entries):
    autogencpp.write("""      if(i==%d) return "%s";\n""" % (i,enum))
  autogencpp.write("      break;\n")
autogencpp.write("  }\n")
autogencpp.write("""  casadi_error("getSchemeEntryEnumName: supplied number is out of range. Scheme '" << getSchemeName(scheme) << "' has only " << getSchemeSize(scheme) << " entries: " << getSchemeEntryNames(scheme) << ".");\n""")
autogencpp.write("}\n")

autogencpp.write("int getSchemeSize(InputOutputScheme scheme) {\n  switch (scheme) {\n")
for p in schemes:
  autogencpp.write("    case SCHEME_%s: \n" % p.enum)
  autogencpp.write("""      return %d;\n""" % len(p.entries))
  autogencpp.write("      break;\n")
autogencpp.write("  }\n")
autogencpp.write("}\n")

autogencpp.write("int getSchemeEntryEnum(InputOutputScheme scheme, const std::string &name) {\n  switch (scheme) {\n")
for p in schemes:
  autogencpp.write("    case SCHEME_%s: \n" % p.enum)
  for i, (name, doc, enum) in enumerate(p.entries):
    autogencpp.write("""      if(name=="%s") return %d;\n""" % (name,i))
  autogencpp.write("      break;\n")
autogencpp.write("  }\n")
autogencpp.write("""  casadi_error("getSchemeEntryEnum: Scheme '" << getSchemeName(scheme) <<  "' has no entry named '" << name <<  "'. Available entries are: " << getSchemeEntryNames(scheme) << ".");\n""")
autogencpp.write("}\n")
autogenmetadatahpp.write("}\n#endif //SCHEMES_METADATA_HPP")
autogenhelpershpp.write("}\n#endif //SCHEMES_HELPERS_HPP")
autogencpp.write("}")

autogenpy.write("#endif //AUTOGENERATED_I\n")

# Endline at the end of files
autogenmetadatahpp.write("\n\n")
autogenhelpershpp.write("\n\n")
autogencpp.write("\n\n")
autogenpy.write("\n")
