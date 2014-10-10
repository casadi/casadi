/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <fstream>

#include "casadi/core/profiling.hpp"
#include "casadi/core/casadi_math.hpp"
#include "casadi/core/casadi_calculus.hpp"

/**float(m.group(1))*1e-6, float(m.group(2))*1e-3, m.group(3), int(m.group(4)), m.group(5) , m.group(6)
linereg = re.compile("^([\d\.]+) ns \| ([\d\.]+) ms \| (0x[\w\d]+:[\w_]+):(\d+)\|(0x[\w\d]+:[\w_]+)?\|(.*)$")

void parseline(const std::string& line,double& ns,double &ms, std::string &id, std::string &name, int & line_number, )
*/

using namespace casadi;


struct linestat {
  double total_time;
  std::string code;
  int count;
  int opcode;
  long dependency;
};

typedef ProfilingData_IO iostat;

struct functionstat {
  std::string name;
  std::vector<linestat> lines;
  double total_time;
  double external_time; // Part of time spent in calls
  double overhead_time; // Part of time not spent in evaluating lines or external calls
  int count;
  int algorithm_size;
  ProfilingData_FunctionType type;
  std::vector<iostat> inputs;
  std::vector<iostat> outputs;
};

typedef std::map<long,functionstat> Stats;

int main(int argc, char* argv[])
{
    Stats data;
    if (argc < 1) {
        std::cerr << "Usage: " << argv[0] << "prof.log" << std::endl;
        return 1;
    }
    
  std::ifstream myfile (argv[1]);
  if (myfile.is_open())
  {
    while (!myfile.eof()) {
    ProfilingHeader hd;
    myfile.read(reinterpret_cast<char*>(&hd), sizeof(hd));
    switch (hd.type) {
     case (ProfilingData_Type_TIMELINE) : {
      ProfilingData_TIMELINE s;
      myfile.read(reinterpret_cast<char*>(&s), sizeof(s));
      Stats::iterator it = data.find(s.thisp);
      std::vector<linestat> & v = it->second.lines;
      v[s.line_number].count+=1;
      v[s.line_number].total_time+=s.local;
      //std::cout << s.thisp << ":" << s.line_number << "|" << s.local << "," << s.total << std::endl;
     }; break;
     case (ProfilingData_Type_SOURCE) : {
      ProfilingData_SOURCE s;
      myfile.read(reinterpret_cast<char*>(&s), sizeof(s));
      char sourceline[s.length+1];
      myfile.read(reinterpret_cast<char*>(&sourceline), s.length);
      Stats::iterator it = data.find(s.thisp);
      std::vector<linestat> & v = it->second.lines;
      linestat L;
      L.total_time = 0;
      L.count = 0;
      L.opcode = 0;
      L.dependency = 0;
      if (s.line_number-int(v.size())+1>0) {
        v.insert(v.end(),s.line_number-int(v.size())+1,L);
      }
      //std::cout << v.size() << s.line_number << std::endl;
      sourceline[s.length] = 0;
      v[s.line_number].code = sourceline;
      v[s.line_number].opcode = s.opcode;
      v[s.line_number].dependency = s.dependency;
      //std::cout << s.thisp << ":" << s.line_number << ": " << sourceline << std::endl;
     }; break;
     case (ProfilingData_Type_NAME) : {
      ProfilingData_NAME s;
      myfile.read(reinterpret_cast<char*>(&s), sizeof(s));
      char name[s.length+1];
      myfile.read(name, sizeof(char)*s.length);
      Stats::iterator it = data.find(s.thisp);
      if (it==data.end()) {
        functionstat f;
        f.count = 0;
        f.total_time = 0;
        f.external_time = 0;
        f.overhead_time = 0;
        data[s.thisp] = f;
        it = data.find(s.thisp);
      }
      name[s.length] = 0; // Null-terminated
      
      functionstat &f = it->second;
      f.name = name;
      f.algorithm_size = s.algorithm_size;
      f.type = s.type;
      
      f.inputs.resize(s.numin);
      //std::cout << name << std::endl;
      //std::cout << "n" << s.numin << "," << s.numout << std::endl;
      for (int i=0;i<s.numin;++i) {
        myfile.read(reinterpret_cast<char*>(&f.inputs[i]), sizeof(iostat));
      }
      f.outputs.resize(s.numout);
      for (int i=0;i<s.numout;++i) {
        myfile.read(reinterpret_cast<char*>(&f.outputs[i]), sizeof(iostat));
      }
     }; break;
     case (ProfilingData_Type_ENTRY) : {
      ProfilingData_ENTRY s;
      myfile.read(reinterpret_cast<char*>(&s), sizeof(s));
      Stats::iterator it = data.find(s.thisp);      
      it->second.count+=1;
      //std::cout << "Entry " << s.thisp << std::endl;
     }; break;
     case (ProfilingData_Type_EXIT) : {
      ProfilingData_EXIT s;
      myfile.read(reinterpret_cast<char*>(&s), sizeof(s));
      Stats::iterator it = data.find(s.thisp);
      it->second.total_time +=s.total;
      //std::cout << "Exit " << s.thisp << ": " << s.total << std::endl;
     }; break;
     default:
       std::cerr << "Unknown type in profile header: " << hd.type << std::endl;
    }
    }
  } else {
    std::cerr << "Unable to open file" << std::endl;
  }
  std::cout << "== Here comes the stats ==" << std::endl;
  
  for (Stats::const_iterator it = data.begin();it!=data.end();++it) {
    const functionstat & f = it->second;
    if (f.count ==0 ) continue;
    std::cout << f.name << ":" << f.count << ":" << f.algorithm_size << ":" << bool(f.type==ProfilingData_FunctionType_MXFunction) << std::endl;
    /**
    for (int i=0;i<f.lines.size();++i) {
     const linestat &L = f.lines[i];
     std::cout << " " << i << ":" << L.count << ":" << L.total_time << ":" << L.code;
     if (L.dependency!=0) std::cout << "| call " << L.dependency << std::endl;
    }*/
  }
  
  // Calculate external_time and overhead_time
  for (Stats::iterator it = data.begin();it!=data.end();++it) {
    functionstat & f = it->second;
    if (f.type==ProfilingData_FunctionType_MXFunction) {
      f.overhead_time = f.total_time;
      for (std::vector<linestat>::const_iterator it2 = f.lines.begin();it2!=f.lines.end();++it2) {
        const linestat & l = *it2;
        f.overhead_time-= l.total_time;
        if (l.dependency) {
          f.external_time+= l.total_time;
        }
      }
    }
  }
  
  std::ofstream report ("prof.html");
  
  report << "<!DOCTYPE html PUBLIC '-//W3C//DTD XHTML 1.0 Strict//EN' 'http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd'>\n<html xmlns='http://www.w3.org/1999/xhtml'>\n<head><script src='http://d3js.org/d3.v3.min.js'></script><script src='http://cpettitt.github.io/project/dagre-d3/v0.1.5/dagre-d3.min.js'></script>";
  
  report << "<style>"
"svg {"
"    overflow: hidden;"
"}"
""
".node rect {"
"    stroke: #333;"
"    stroke-width: 1.5px;"
"    fill: #fff;"
"}"
""
".edgeLabel rect {"
"    fill: #fff;"
"}"
""
".edgePath {"
"    stroke: #333;"
"    stroke-width: 1.5px;"
"    fill: none;"
"}"
".outer {"
"    width: 1024px;"
"    height: 960px;"
"    overflow: auto;"
"}"
".inner {"
"    width: 8000px;"
"    height: 6000px;"
"}"
"svg {"
"    display: block;"
"    width: 100%;"
"    height: 100%;"
"}"
  "</style>";
  
  report << "</head><body>\n";
  report << "<div id='chart3'></div>" << std::endl;
  


  report << "<table><thead><tr><th>Id</th><th>#calls</th><th>Internal (s)</th><th>External (s)</th><th>Overhead (s)</th><th>Algorithm size</th><th>#inputs</th><th>#input nz</th><th>#outputs</th><th>#output nz</th></tr></thead>" << std::endl;
  

  
  for (Stats::const_iterator it = data.begin();it!=data.end();++it) {
    const functionstat & f = it->second;
    if (f.count ==0 ) continue;
    int nz_in = 0;
    int nz_out = 0;
    for (int i=0;i<f.inputs.size();++i) {
      nz_in+= f.inputs[i].ndata;
    }
    for (int i=0;i<f.outputs.size();++i) {
      nz_out+= f.outputs[i].ndata;
    }
    
    report << "<tr><td><a href='#" << it->first << "'>" << f.name << "</a></td><td>" << f.count << "</td><td>"  << f.total_time-f.external_time-f.overhead_time << "</td><td>"  << f.external_time << "</td><td>"  << f.overhead_time << "</td><td>" << f.algorithm_size <<  "</td><td>" << f.inputs.size() <<  "</td><td>" << nz_in <<  "</td><td>" << f.outputs.size() <<  "</td><td>" << nz_out <<  "</td></tr>\n";
  }
  
  report << "</table>";
  
  // data
  report << "<script>var functions={";
  for (Stats::const_iterator it = data.begin();it!=data.end();++it) {
    const functionstat & f = it->second;
    
    if (f.type== ProfilingData_FunctionType_Other ) continue;

    report << it->first  << ": {name:" << '"' << f.name << '"' << "}," << std::endl;
  }
  report << "};" << std::endl;
  
  report << "var functioncalls=[";
  for (Stats::const_iterator it = data.begin();it!=data.end();++it) {
    const functionstat & f = it->second;
    if (f.type== ProfilingData_FunctionType_Other ) continue;
    
    std::map<long,int> callmap;
    for (int i=0;i<f.lines.size();++i) {
      if (f.lines[i].dependency!=0) {
      
        Stats::iterator itd = data.find(f.lines[i].dependency);
        
        if (itd->second.type== ProfilingData_FunctionType_Other ) continue;
        
        std::map<long,int>::iterator it = callmap.find(f.lines[i].dependency);
        if (it==callmap.end()) {
          callmap[f.lines[i].dependency] = f.lines[i].count;
        } else {
          it->second += f.lines[i].count;
        }
      }
    }
    
    for (std::map<long,int>::iterator it2=callmap.begin();it2!=callmap.end();++it2) {
      report << " {source:" << it->first << ", target:" << it2->first << "," << "ncalls:" << it2->second << "},";
    }
    
    
  }
  report << "];\n</script>\n";
  
  
  
  //Binning
  std::vector<int> type_binning_ncalls(casadi::NUM_BUILT_IN_OPS);
  std::vector<double> type_binning_time(casadi::NUM_BUILT_IN_OPS);
  
  for (Stats::const_iterator it = data.begin();it!=data.end();++it) {
    const functionstat & f = it->second;
    if (f.type==ProfilingData_FunctionType_MXFunction) {
      for (int i=0;i<f.lines.size();++i) {
        const linestat &L = f.lines[i];
        type_binning_ncalls[L.opcode]+=1;
        type_binning_time[L.opcode]+=L.total_time;
      }
    }
  }
  

  report << "<div class='outer'><div class='inner'><svg>"
"    <g transform='translate(20,20)'/>"
"</svg></div></div>";
  
  report << "<script>"
"var g = new dagreD3.Digraph();"
"Object.keys(functions).forEach(function(key) {"
"g.addNode(   key.toString(), { label: functions[key].name });"
"});"
"functioncalls.forEach(function(link) {"
"g.addEdge(null, link.source.toString(),   link.target.toString(), { label: link.ncalls.toString() });"
"});"
"g=g.filterNodes(function(u) { return g.neighbors(u).length >0 });"
"var layout = dagreD3.layout()"
"                    .nodeSep(20)"
"                    .rankDir('LR');"
"var renderer = new dagreD3.Renderer();"
"renderer.layout(layout).run(g, d3.select('svg g'));"
"</script>";
  
  report << "<table><thead><tr><th>Operation code</th><th>Operation</th><th>total (s)</th><th>ncalls</th></tr></thead>\n";
  for (int i=0;i<casadi::NUM_BUILT_IN_OPS;++i) {
    report << "<tr><td>" << i << "</td><td>";
    casadi::casadi_math<double>::printPre(i,report);
    casadi::casadi_math<double>::printSep(i,report);
    casadi::casadi_math<double>::printPost(i,report);
    report << "</td><td>" << type_binning_time[i] << "</td><td>" << type_binning_ncalls[i] <<  "</td><td>";
  }
  report << "</table>";
  
  
  for (Stats::const_iterator it = data.begin();it!=data.end();++it) {
    const functionstat & f = it->second;
    if (f.count ==0 ) continue;
    report << "<a name='" << it->first << "'><h2>" << f.name << "</h2></a><dl><dt>#calls</dt><dd>" << f.count << "</dd><dt>Algorithm size</dt><dd>" << f.algorithm_size << "</dd><dt>Total (s)</dt><dd>" << f.total_time << "</dd>" << std::endl;
    report << "<dt>Inputs</dt><dd><table><thead><tr><td>i</td><td>rows</td><td>cols</td><td>nonzeros</td></tr></thead>" << std::endl;
    for (int i=0;i<f.inputs.size();++i) {
      report << "<tr><td>" << i << "</td><td>" << f.inputs[i].nrow << "</td><td>" << f.inputs[i].ncol <<  "</td><td>" << f.inputs[i].ndata << "</td></tr>" << std::endl;
    }
    report << "</table></dd>" << std::endl;
    report << "<dt>Outputs</dt><dd><table><thead><tr><td>i</td><td>rows</td><td>cols</td><td>nonzeros</td></tr></thead>" << std::endl;
    for (int i=0;i<f.outputs.size();++i) {
      report << "<tr><td>" << i << "</td><td>" << f.outputs[i].nrow << "</td><td>" << f.outputs[i].ncol <<  "</td><td>" << f.outputs[i].ndata << "</td></tr>" << std::endl;
    }
    report << "</table></dd>" << std::endl;
    
    report << "</dl>\n";

    if (f.type==ProfilingData_FunctionType_MXFunction || f.type==ProfilingData_FunctionType_Other) {
      report << "<table><thead><tr><th>Codeline</th><th>total (ms)</th><th>ncalls</th><th>souce</th></tr></thead>\n";
      int total = f.lines.size();
      if (total>10000) total = 10000;
      for (int i=0;i<total;++i) {
       const linestat &L = f.lines[i];
       //if (L.dependency!=0) std::cout << "| call " << L.dependency << std::endl;
       report << "<tr><td>" << i << "</td><td>" << L.total_time*1000 << "</td><td>" << L.count <<  "</td><td>";
       if (L.dependency==0) {
         report << L.code << "</td></tr>" << std::endl;
       } else {
         report << "<a href='#"<< L.dependency << "'>" << L.code << "</a></td></tr>" << std::endl;
       }
       
        
      }
      report << "</table>";
      if (f.lines.size()>10000) {
        report << "<p>Clipped " << f.lines.size() << " lines</p>" << std::endl;
      }
    }

  }
  
  report << "</body></html>";
  
}
