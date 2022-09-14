import os
import re
import baseconv # pip install python-baseconv

next_identifier = 0

try:
    with open("next_identifier.txt","r") as inp:
        next_identifier = int(inp.read())
except:
    pass
    
print("next_identifier",next_identifier)
    
def id2str(a):
    return baseconv.base36.encode(a)
def str2id(a):
    return int(baseconv.base36.decode(a))

re_comments = re.compile(r"^ */\*\* (.*?)\*/$", re.MULTILINE | re.DOTALL)
identifiers = set()
location = {}


for root, dirs, files in os.walk("../casadi"):
   for name in files:
      if "runtime" in root: continue
      if "core" not in root: continue
      if name.endswith(".hpp"):
          file = os.path.join(root, name)
          with open(file,"r") as inp:
            data = inp.read()
            orig = data
            def process(match):
                global next_identifier
                if "\\brief" in match.group(0):
                 if "\\identifier" in match.group(0):
                    identifier = str2id(re.search(r"\\identifier\{(.*?)\}",match.group(0)).group(1))
                    line = data[:match.start()].count("\n")+1
                    loc = file + ":" + str(line) + "\n" + match.group(0)
                    if identifier in identifiers:
                        raise Exception(f"Duplicate identifiers:\n{loc}\n{location[identifier]}")
                    identifiers.add(identifier)
                    location[identifier] = loc
                 else:
                    m = match.group(0).rstrip()
                    ident = m.find("\\brief")*" "
                    next_identifier+=1
                    return m[:-2].rstrip() + "\n" + ident+ "\\identifier{" + id2str(next_identifier) + "}" + " */"
                return match.group(0)
            data = re.sub(re_comments, process, data)
          if data!=orig:
            print("Updating", file)
            with open(file,"w") as outp:
              outp.write(data)
            #print(data)

print("next_identifier",next_identifier)
with open("next_identifier.txt","w") as outp:
    outp.write(str(next_identifier))
