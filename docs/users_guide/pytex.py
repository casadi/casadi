import glob
import subprocess
import os

# Process all snippets in the current directory
for f in glob.glob("pytex_*.py"):
  logfilename = f[:-3]+'.log'
  logfile = file(logfilename,'w')
  
  # Execute snippet
  p = subprocess.Popen(['python',f],stdout=logfile)
  p.wait()
  
  logfile = file(logfilename,'r')
  
  outfile = False
  contents = False
  for l in logfile.readlines():
    if l.startswith('pytex snippet::'):
      num=l[len('pytex snippet::'):].rstrip()
      
      # Remove logfiles that are empty
      #if outfile and not(contents):
      #  outfile.close()
      #  os.remove(outfile.name)
    
      outfile = file(f[:-3]+'_' + num + '.log','w')
      contents = False
    else:
      if outfile:
        contents = True
        outfile.write(l)
  os.remove(logfilename)
  # Remove logfiles that are empty
  #if outfile and not(contents):
  #  outfile.close()
  #  os.remove(outfile.name)
