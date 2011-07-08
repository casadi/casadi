import os
import sys
from subprocess import *
import time
import re



def is_exe(root,name):
  fpath = os.path.join(root,name)
  return os.path.exists(fpath) and os.access(fpath, os.X_OK)
    
class TestSuite:
  def __init__(self,suffix=None,dirname=None,preRun=None,postRun=None,command=None,skipdirs=[],skipfiles=[],inputs={},workingdir = lambda x: x,allowable_returncodes=[],args=[]):
    """
    
    dirname: The directory that should be crawled for test problems.
    
    command: A mapping of (dirname,filename) -> Popen command list.
             e.g  lambda dir,fn:  ['./'+fn]  for executable files
             
    skipdirs: A list of directories that should be skipped during recursive traversal of dirname.
    
    skipfiles: A list of matched files that should be skipped during recursive traversal of dirname.
    
    suffix:  Only the files that end with this suffix will be tested.
             Omit argument or set to None to test any file.
             
    inputs:  Can be a dictionary of filename -> string to pass numbers to STDIN of the process
    
    workingdir: Specify the working directory for the process. The path may be relative to the /trunk/test path
    
    args: a list of command line options:
       -skipfiles="file1 file2"   get's added to skipfiles
       
    
    """
    
    # Don't buffer
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
    sys.stderr = os.fdopen(sys.stderr.fileno(), 'w', 0)

    self.stats={'numtests':0,'numfails':0}
    self.suffix = suffix
    self.dirname = dirname
    self.preRun = preRun
    self.postRun = postRun
    self.command = command
    self.skipdirs = skipdirs
    self.skipfiles = skipfiles
    self.inputs = inputs
    self.allowable_returncodes = allowable_returncodes 
    self.workingdir = workingdir
    self.args=args
    for arg in args:
      okay = False
      m = re.search('-skipfiles=(.*)', arg)
      if m:
        print "foo:", m.group(1)
        self.skipfiles+=m.group(1).split(' ')
        okay = True
        
      if not(okay):
        print "Unknown argument: ", arg

  def run(self):
    print "Running test in " + self.dirname
    for root, dirs, files in os.walk(self.dirname):
      if not(self.preRun is None):
        self.preRun(root)
      for name in files:
        for skip in self.skipdirs:
          if skip in dirs:
            dirs.remove(skip)
        if (self.suffix is None and is_exe(root,name)) or (not(self.suffix is None) and name.endswith('.'+self.suffix)):
          if not(name in self.skipfiles):
            self.test(root,name,self.command)
      if not(self.postRun is None):
        self.postRun(root)
        
    print "Ran %d tests, %d fails." % (self.stats['numtests'],self.stats['numfails'])

    if self.stats['numfails']>0:
      sys.exit(1)


  def test(self,dir,fn,commands):
    self.stats['numtests']+=1
    print ("%02d. " % self.stats['numtests']) + fn
    t0 = time.clock()
    p=Popen(self.command(dir,fn),cwd=self.workingdir(dir),stdout=PIPE, stderr=PIPE, stdin=PIPE)
    inp = None
    if fn in self.inputs:
      inp = self.inputs[fn]
    stdoutdata, stderrdata = p.communicate(inp)
    t = time.clock() - t0
    if (p.returncode==0 or p.returncode in self.allowable_returncodes):
      pass
      #print "  > Succes: %0.2f [ms]" % (t*1000)
    else :
      self.stats['numfails']+=1
      print "In %s, %s failed:" % (dir,fn)
      print "="*30
      print p.returncode
      print stdoutdata
      print stderrdata
      print "="*30
  
