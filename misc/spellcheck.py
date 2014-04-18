from pyparsing import *
import hunspell
hobj = hunspell.HunSpell('/usr/share/hunspell/en_US.dic', '/usr/share/hunspell/en_US.aff')
import subprocess

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

success = True

for f in files:
  localsuccess = True
  ff = os.path.join(dir,f)
  if f.endswith("cpp"):
    continue
  text = ""
  spells = ""
  for (m,s,e) in cppStyleComment.scanString(file(ff,'r').read()):
    t = m.asList()[0]
    #if "//" in t[3:-3]:
    #  continue
    t = re.sub('\[\w+\]','',t)
    t = re.sub('\(\w+\s+x\s+\w+\)','',t)
    t = re.sub('[\\\\@](copydoc|a|e|param) \w+','',t)
    t = re.sub(r'\b[A-Z_]{2,}\w+','',t)
    t = re.sub(r'(?<!\w)-+(?!\w)','',t)
    t = re.sub('\\\\verbatim(.*?)\\\\endverbatim','',t,flags=re.DOTALL)
    t = re.sub('\\\\f\$(.*?)\\\\f\$','',t,flags=re.DOTALL)
    t = re.sub('\\\\f\[(.*?)\\\\f\]\$','',t,flags=re.DOTALL)
    
    t = re.sub(r'\b\w+[A-Z]\w+\b','',t) # camelcase
    t = re.sub(r'\b\w+_\w+\b','',t) # camelcase
    t = re.sub('[\\\\@][\w{}]+','',t)
    t = re.sub('<tt>.*?</tt>','',t)
    t = re.sub(r'\b[ntdpcx][a-zA-Z0-9]\b',lambda e: e.group(0) if hobj.spell(e.group()) else '',t) # should be escaped with \e
    t = re.sub('</?\w+>','',t)
    t = re.sub('#\w+','',t)
    t = re.sub('\*','',t)
    t = re.sub('[\w\d]+-by-[\w+\d]','',t)
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
    print "In file ", ff
    print "\n".join(out.split("\n")[2:-2])
    print outh
    success = False
    
if success:
  sys.exit(0)
else:
  sys.exit(1)
