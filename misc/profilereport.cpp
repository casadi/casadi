#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <fstream>

#include "symbolic/profiling.hpp"
#include "symbolic/casadi_math.hpp"
#include "symbolic/casadi_calculus.hpp"

/**float(m.group(1))*1e-6, float(m.group(2))*1e-3, m.group(3), int(m.group(4)), m.group(5) , m.group(6)
linereg = re.compile("^([\d\.]+) ns \| ([\d\.]+) ms \| (0x[\w\d]+:[\w_]+):(\d+)\|(0x[\w\d]+:[\w_]+)?\|(.*)$")

void parseline(const std::string& line,double& ns,double &ms, std::string &id, std::string &name, int & line_number, )
*/



struct linestat {
  double total_time;
  std::string code;
  int count;
  int opcode;
  long dependency;
};

struct functionstat {
  std::string name;
  std::vector<linestat> lines;
  double total_time;
  int count;
  int algorithm_size;
  ProfilingData_FXType type;
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
        data[s.thisp] = f;
        it = data.find(s.thisp);
      }
      name[s.length] = 0; // Null-terminated
      it->second.name = name;
      it->second.algorithm_size = s.algorithm_size;
      it->second.type = s.type;
      //std::cout << s.thisp << ": " << s.length << std::endl;
      //std::cout << s.thisp << ": " << name << std::endl;
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
    std::cout << f.name << ":" << f.count << ":" << f.algorithm_size << ":" << bool(f.type==ProfilingData_FXType_MXFunction) << std::endl;
    /**
    for (int i=0;i<f.lines.size();++i) {
     const linestat &L = f.lines[i];
     std::cout << " " << i << ":" << L.count << ":" << L.total_time << ":" << L.code;
     if (L.dependency!=0) std::cout << "| call " << L.dependency << std::endl;
    }*/
  }
  
  std::ofstream report ("prof.html");
  
  report << "<!DOCTYPE html PUBLIC ""-//W3C//DTD XHTML 1.0 Strict//EN"" ""http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd"">\n<html xmlns=""http://www.w3.org/1999/xhtml"">\n<body>\n<img src=""callgraph.png""/>";

  report << "<table><thead><tr><th>Id</th><th>#calls</th><th>Algorithm size</th><th>Total (s)</th></tr></thead>" << std::endl;
  
  for (Stats::const_iterator it = data.begin();it!=data.end();++it) {
    const functionstat & f = it->second;
    if (f.count ==0 ) continue;
    report << "<tr><td><a href='#" << it->first << "'>" << f.name << "</a></td><td>" << f.count << "</td><td>" << f.algorithm_size <<  "</td><td>" << f.total_time <<  "</td></tr>\n";
  }
  
  report << "</table>";
  
  //Binning
  std::vector<int> type_binning_ncalls(CasADi::NUM_BUILT_IN_OPS);
  std::vector<double> type_binning_time(CasADi::NUM_BUILT_IN_OPS);
  
  for (Stats::const_iterator it = data.begin();it!=data.end();++it) {
    const functionstat & f = it->second;
    if (f.type==ProfilingData_FXType_MXFunction) {
      for (int i=0;i<f.lines.size();++i) {
        const linestat &L = f.lines[i];
        type_binning_ncalls[L.opcode]+=1;
        type_binning_time[L.opcode]+=L.total_time;
      }
    }
  }
  
  report << "<table><thead><tr><th>Operation code</th><th>Operation</th><th>total (s)</th><th>ncalls</th></tr></thead>\n";
  for (int i=0;i<CasADi::NUM_BUILT_IN_OPS;++i) {
    report << "<tr><td>" << i << "</td><td>";
    CasADi::casadi_math<double>::printPre(i,report);
    CasADi::casadi_math<double>::printSep(i,report);
    CasADi::casadi_math<double>::printPost(i,report);
    report << "</td><td>" << type_binning_time[i] << "</td><td>" << type_binning_ncalls[i] <<  "</td><td>";
  }
  report << "</table>";
  
  
  for (Stats::const_iterator it = data.begin();it!=data.end();++it) {
    const functionstat & f = it->second;
    if (f.count ==0 ) continue;
    report << "<a name='" << it->first << "'><h2>" << f.name << "</h2></a><dl><dt>#calls</dt><dd>" << f.count << "</dd><dt>Algorithm size</dt><dd>" << f.algorithm_size << "</dd><dt>Total (s)</dt><dd>" << f.total_time << "</dd></dl>\n";

    if (f.type==ProfilingData_FXType_MXFunction || f.type==ProfilingData_FXType_Other) {
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
