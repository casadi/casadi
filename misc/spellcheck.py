from pyparsing import *
import hunspell
hobj = hunspell.HunSpell('/usr/share/hunspell/en_US.dic', '/usr/share/hunspell/en_US.aff')
import subprocess
from multiprocessing import Process, Queue, Lock

lt_exclusions = [
  'WHITESPACE_RULE',
  'EN_A_VS_AN',
  'COMMA_PARENTHESIS_WHITESPACE',
  'GENERAL_XX',
  'SENTENCE_WHITESPACE',
  'EN_QUOTES',
  'DT_DT',
  'A_PLURAL',
  'DOUBLE_PUNCTUATION',
  'ALL_OF_THE',
  'ENGLISH_WORD_REPEAT_BEGINNING_RULE',
  'UPPERCASE_SENTENCE_START',
  'ENGLISH_WORD_REPEAT_RULE',
  'IN_A_X_MANNER',
  'HE_VERB_AGR',
  'PHRASE_REPETITION',
  'KIND_OF_A',
  'EN_UNPAIRED_BRACKETS',
  'MAY_BE'
]

import sys
import re
import os

dir = sys.argv[1]
files = sys.argv[2:]

if 'symbolic' in dir:
  sys.exit(0)

success = True

def checkfile(f,lock,queue):
  localsuccess = True
  text = ""
  spells = ""
  interface = None
  m = re.search("interfaces/(\w+)/",f)
  if m:
    interface = re.compile(r"\b"+m.group(1)+r"\b",re.I)
  for (m,s,e) in cppStyleComment.scanString(file(f,'r').read()):
    t = m.asList()[0]
    if t.startswith("/*") and not t.startswith("/**"):
      continue
    #if "//" in t[3:-3]:
    #  continue
    t = re.sub('\[\w+\]','',t)
    t = re.sub('\(\w+\s+x\s+\w+\)','',t)
    t = re.sub('[\\\\@](copydoc|a|e|p|param|defgroup)\s+\w+','',t)
    t = re.sub(r'\b[A-Z_]{2,}\w+','',t)
    t = re.sub(r'(?<!\w)-+(?!\w)','',t)
    t = re.sub('\\\\verbatim(.*?)\\\\endverbatim','',t,flags=re.DOTALL)
    t = re.sub('\\\\f\$(.*?)\\\\f\$','',t,flags=re.DOTALL)
    t = re.sub('\\\\f\[(.*?)\\\\f\]\$','',t,flags=re.DOTALL)
    
    t = re.sub(r'\b\w+[A-Z]\w+\b','',t) # camelcase
    t = re.sub(r'\b\w+_\w+\b','',t) # camelcase
    t = re.sub('[\\\\@][\w{}]+','',t)
    t = re.sub('<tt>.*?</tt>','',t)
    t = re.sub('</?\w+>','',t)
    t = re.sub('#\w+','',t)
    t = re.sub('\*','',t)
    t = re.sub('\b[\w\d]+(-by-|x)[\w+\d]\b','',t)
    t = re.sub('C\d+','',t) # Visual studio warnings
    t = re.sub("'[^\s]*?'",'',t)
    t = re.sub(r'\b[ntdpcx][a-zA-Z0-9]\b',lambda e: e.group(0) if hobj.spell(e.group()) else '',t) # should be escaped with \e
    t = re.sub("-\d+",'',t)
    t = re.sub("\b\w+\d+\b",'',t)
    if interface is not None:
      t = interface.sub('',t)
    text+= t+"\n"
    #for w in t.split(" "):
    #  if not hobj.spell(w):
    #    success = False
    #    spells+="Spelling error: %s (%s)\n" % (w,str(hobj.suggest(w)))
  p = subprocess.Popen(['java','-jar','/home/jg/programs/LanguageTool-2.5/languagetool-commandline.jar','-l','en','-d',",".join(lt_exclusions),'-'],stdin = subprocess.PIPE, stdout = subprocess.PIPE)
  out, err = p.communicate(text)
  if len(out.split("\n")[2:-2]) > 0:
    localsuccess = False
  ph = subprocess.Popen(['hunspell','-l','-p','casadi.dic'],stdin = subprocess.PIPE, stdout = subprocess.PIPE)
  outh, errh = ph.communicate(text)
  
  if len(outh) > 0:
    localsuccess = False
  if not localsuccess:
    lock.acquire()
    print "In file ", f
    print "\n".join(out.split("\n")[2:-2])
    print outh
    lock.release()
    return False
  

lock = Lock()
queue = Queue()
for f in files:
  localsuccess = True
  ff = os.path.join(dir,f)
  if f.endswith("cpp"):
    continue
  p = Process(target=checkfile, args=(ff,lock,queue))
  p.start()
  
success = True
    
if success:
  sys.exit(0)
else:
  sys.exit(1)
