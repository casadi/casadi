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
from __future__ import print_function
import os
import sys
from subprocess import *
import time
import re
import tempfile

class TimeoutEvent(Exception):
    pass

def alarm_handler(signum, frame):
    raise TimeoutEvent


def is_exe(root,name):
  fpath = os.path.join(root,name)
  return os.path.exists(fpath) and os.access(fpath, os.X_OK)

from os import kill
from signal import signal


try:
    from signal import alarm, SIGALRM, SIGKILL
    alarm_available = True
except:
    def alarm(amount):
        pass
    alarm_available = False

from subprocess import PIPE, Popen

# Snippet from http://stackoverflow.com/questions/1191374/subprocess-with-timeout
def run(args, input=None, cwd = None, shell = False, kill_tree = True, timeout = -1, env = None):
    '''
    Run a command with a timeout after which it will be forcibly
    killed.
    '''
    class Alarm(Exception):
        pass
    def alarm_handler(signum, frame):
        raise Alarm
    p = Popen(args, shell = shell, cwd = cwd, stdout = PIPE, stderr = PIPE, env = env)
    if timeout != -1:
        if alarm_available:
            signal(SIGALRM, alarm_handler)
        alarm(timeout)
    try:
        t0 = time.time()
        stdout, stderr = p.communicate()
        print("Ran for", time.time()-t0, "seconds.")
        if timeout != -1:
            alarm(0)
    except Alarm:
        pids = [p.pid]
        if kill_tree:
            pids.extend(get_process_children(p.pid))
        for pid in pids:
            # process might have died before getting to this line
            # so wrap to avoid OSError: no such process
            try:
                kill(pid, SIGKILL)
            except OSError:
                pass
        return -9, '', ''
    return p.returncode, stdout, stderr

def killProcess(pid):
  pids = [pid]
  if kill_tree:
      pids.extend(get_process_children(p.pid))
  for pid in pids:
      # process might have died before getting to this line
      # so wrap to avoid OSError: no such process
      try:
          kill(pid, SIGKILL)
      except OSError:
          pass

def get_process_children(pid):
    p = Popen('ps --no-headers -o pid --ppid %d' % pid, shell = True,
              stdout = PIPE, stderr = PIPE)
    stdout, stderr = p.communicate()
    return [int(p) for p in stdout.split()]

if __name__ == '__main__':
    print(run('find /', shell = True, timeout = 3))
    print(run('find', shell = True))

deprecated = re.compile(r"\b[dD]epr[ei]c[ie]?at[ei]")
warning = re.compile("warning")

class TestSuite:
  def __init__(self,suffix=None,dirname=None,preRun=None,postRun=None,command=None,skipdirs=[],skipfiles=[],inputs={},workingdir = lambda x: x,allowable_returncodes=[],args=[],stderr_trigger=[],stdout_trigger=[],check_depreciation=True,check_warning=False,custom_stdout=None,default_fail=False):
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


    stderr_trigger: list of strings/regexes. If any string is found in std_err, it is considered an error, regardless of the return_code. A message may be attached to a trigger by packing the string/regex in a 2-tuple, with the second element the message.

    stdout_trigger: list of strings/regexes. If any string is found in std_out, it is considered an error, regardless of the return_code. A message may be attached to a trigger by packing the string/regex in a 2-tuple, with the second element the message.

    check_depreciation: raise an error if the output contains a depreciation warning

    check_warning: raise an error if the output contains a warning

    args: a list of command line options:
       -skipfiles="file1 file2"   get's added to skipfiles
       -memcheck           Include a check for memory leaks
       -passoptions="option1 option2"   get's passed onto the command constructor as third argument


    """
    if alarm_available:
        signal(SIGALRM, alarm_handler)

    # Don't buffer
    if sys.version_info < (3, 0):
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
    self.memcheck = False
    self.passoptions = []
    self.stderr_trigger = stderr_trigger
    self.stdout_trigger = stdout_trigger
    self.custom_stdout = custom_stdout
    self.default_fail=default_fail
    if check_depreciation:
      self.stderr_trigger.append((deprecated,"deprecated"))
      self.stdout_trigger.append((deprecated,"deprecated"))
    if check_warning:
      self.stderr_trigger.append((warning,"warning"))
      self.stdout_trigger.append((warning,"warning"))
    self.args=args
    for arg in args:
      okay = False
      m = re.search('-skipfiles=(.*)', arg)
      if m:
        self.skipfiles+=m.group(1).split(' ')
        okay = True
      m = re.search('-skipdirs=(.*)', arg)
      if m:
        self.skipdirs+=m.group(1).split(' ')
        okay = True
      m = re.search('-memcheck', arg)
      if m:
        self.memcheck = True
        okay = True
      m = re.search('-passoptions=(.*)', arg)
      if m:
        self.passoptions+=m.group(1).split(' ')
        okay = True
      if not(okay):
        print("Unknown argument: ", arg)

  def run(self):
    print("Running test in " + self.dirname)
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

    print("Ran %d tests, %d fails." % (self.stats['numtests'],self.stats['numfails']))

    if self.stats['numfails']>0:
      sys.exit(1)


  def test(self,dir,fn,commands):
    self.stats['numtests']+=1
    print(("%02d. " % self.stats['numtests']) + fn)
    t0 = time.time()

    p=Popen(self.command(dir,fn,self.passoptions),cwd=self.workingdir(dir),stdout=PIPE, stderr=PIPE, stdin=PIPE)

    inp = None
    if callable(self.inputs):
      inputs = self.inputs(dir,fn)
    else:
      inputs = self.inputs
    if fn in inputs:
      inp = inputs[fn]


    alarm(60*60) # 1 hour
    try:
      stdoutdata, stderrdata = p.communicate(inp)
    except TimeoutEvent:
      killProcess(p.pid)
      raise Exception("Timout.")
    try:
      stdoutdata = stdoutdata.decode("ascii")
      stderrdata = stderrdata.decode("ascii")
    except:
      pass
    alarm(0) # Remove alarm
    t = time.time() - t0
    if self.custom_stdout is not None:
      stdoutdata = self.custom_stdout(dir,fn)
    
    print("Ran for",t, "seconds")

    stderr_trigger = False
    for trigger in self.stderr_trigger:
      trigger_message = str(trigger)
      if isinstance(trigger,tuple):
        trigger_message = trigger[1]
        trigger = trigger[0]
      trigger_raised = re.search(trigger, stderrdata)
      if trigger_raised:
        print("stderr_trigger '%s' was raised." % trigger_message)
        stderr_trigger = True
    stdout_trigger = False
    for trigger in self.stdout_trigger:
      trigger_message = str(trigger)
      if isinstance(trigger,tuple):
        trigger_message = trigger[1]
        trigger = trigger[0]
      trigger_raised = re.search(trigger, stdoutdata)
      if trigger_raised:
        print("stdout_trigger '%s' was raised." % trigger_message)
        stdout_trigger = True
        
    if self.default_fail:
      passed = stdout_trigger or stderr_trigger
    else:
      passed = (not(stderr_trigger) and (p.returncode==0 or (p.returncode in self.allowable_returncodes)))
    
    if passed:
      pass
      #print "  > Succes: %0.2f [ms]" % (t*1000)
    else :
      self.stats['numfails']+=1
      print("In %s, %s failed:" % (dir,fn))
      print("="*30)
      print("returncode: ", p.returncode)
      print(stdoutdata)
      print(stderrdata)
      print("="*30)
      return

    if self.memcheck:
      # --suppressions=../internal/valgrind-python.supp
      suppressions = ["internal/valgrind-python.supp","internal/valgrind-casadi.supp"]
      supps = ['--suppressions='+os.path.join(os.getcwd(),s) for s in suppressions]
      stdoutfile = tempfile.TemporaryFile()
      p=Popen(['valgrind','--leak-check=full']+supps+['--show-possibly-lost=no','--error-limit=no',"--gen-suppressions=all"]+self.command(dir,fn,self.passoptions),cwd=self.workingdir(dir),stdout=stdoutfile, stderr=PIPE, stdin=PIPE)
      f=Popen(['grep','-E','-C','50', "definitely lost|leaks|ERROR SUMMARY|Invalid write|casadi"],stdin=p.stderr,stdout=PIPE)
      p.stderr.close()
      t0 = time.time()
      alarm(60*60) # 1 hour
      try:
        stdoutdata, stderrdata = f.communicate()
        stdoutdata = stdoutfile.read() + "\n"+ stdoutdata
      except TimeoutEvent:
        killProcess(p.pid)
        killProcess(f.pid)
        raise Exception("Timeout.")
      alarm(0) # Remove alarm
      print("Ran for", time.time()-t0, "seconds")
      m = re.search('definitely lost: (.*) bytes', stdoutdata)
      lost = "0"
      if m:
        lost = m.group(1)
      else:
        m = re.search("no leaks are possible",stdoutdata)
        if not(m):
          if p.returncode==0:
            return
          debug = str(stdoutdata) + " -- " + str(stderrdata) + " -- " + str(f.returncode)
          raise Exception("valgrind output is not like expected: %s" % debug)

      m = re.search('ERROR SUMMARY: (.*) errors', stdoutdata)
      errors = "0"
      if m:
        errors = m.group(1)
      else:
        print(stdoutdata)
        raise Exception("valgrind output is not like expected: %s")

      error_casadi = False
      # Filter out stuff after address header
      diagnose_lines = []
      error_log = True
      for l in stdoutdata.split("\n"):
        if l.startswith("=="):
          if "  Address" in l:
            error_log = False
          elif "  at" in l or "  by" in l:
            pass
          else:
            error_log = True

        if error_log:
          diagnose_lines.append(l)
          if l.startswith("==") and re.search('casadi', l):
            error_casadi = True

      diagnosis = "\n".join(diagnose_lines)

      errors = "0"  # disabling valgrind error-checking for now: samples are flooded with errors
      if not(lost=="0" and errors=="0") or error_casadi:
        if not(lost=="0"):
          print("Memory leak: lost %s bytes" % (lost))
        if not(errors=="0"):
          print("Valgrind errors: %s" % (errors))
        print("="*30)
        print(stdoutdata)
        print("="*30)
        self.stats['numfails']+=1

