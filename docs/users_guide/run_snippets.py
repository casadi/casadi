import glob
import subprocess
import os
import sys

langs = ["python","octave"]

if len(sys.argv)>1:
  langs = sys.argv[1:]

success = True

for lang in langs:
    print("%s snippets" % lang)
    snippets = glob.glob('snippets/*.%s.in' % lang)
    for f in snippets:
       out_file = f.replace("in","out")
       with open(f.replace("in","out"), 'w') as myoutput:
         if lang == "python":
           args = ["python",f]
         if lang == "octave":
           args = ["octave-cli","--eval","try,run('%s'),catch e,e.message,exit(1),end" % f]
         
         p = subprocess.Popen(args,stdout=myoutput, stderr=myoutput)
         p.wait()
       if p.returncode!=0:
         print(("="*5) + f)
         print(open(f,'r').read())
         print("="*50)
         print(open(out_file,'r').read())
         success = False
    print("Processed %d snippets" % len(snippets))

if not success:
  os.exit(1)
