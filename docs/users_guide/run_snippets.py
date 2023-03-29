import glob
import subprocess
import os
import sys
import re

langs = ["python","octave"]

if len(sys.argv)>1:
  langs = sys.argv[1:]

success = True

for lang in langs:
    print("%s snippets" % lang)
    snippets = glob.glob('snippets/*.%s.in' % lang)
    for f in snippets:
       out_file = f.replace("in","out")
       if lang == "python":
         args = ["python",f]
       if lang == "octave":
         args = ["octave","--no-gui","--no-window-system","--eval","try,run('%s'),catch e,e.message,exit(1),end" % f]
     
       p = subprocess.run(args,stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
       with open(f.replace("in","out"), 'w') as myoutput:
       
         out = p.stdout
         # Filter out Ipopt copyright notice
         out = re.sub(r'\*{5,}.*Ipopt.*\*{5,}','',out,flags=re.MULTILINE|re.DOTALL)
         # Filter out qpOASES copyright notice
         out = re.sub(r"""qpOASES -- An Implementation of the Online Active Set Strategy.
Copyright \(C\) 2007-2015 by Hans Joachim Ferreau, Andreas Potschka,
Christian Kirches et al. All rights reserved.

qpOASES is distributed under the terms of the 
GNU Lesser General Public License 2.1 in the hope that it will be 
useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
See the GNU Lesser General Public License for more details.""",'',out,flags=re.MULTILINE|re.DOTALL)

         # remove leading empty lines
         out_split = out.split("\n")
         for i,line in enumerate(out_split):
            line = line.rstrip()
            if len(line)!=0:
                out = "\n".join(out_split[i:])
                break
         
         myoutput.write(out)
       if p.returncode!=0:
         print(("="*5) + f)
         print(open(f,'r').read())
         print("="*50)
         print(open(out_file,'r').read())
         success = False
    print("Processed %d snippets" % len(snippets))

if not success:
  exit(1)
