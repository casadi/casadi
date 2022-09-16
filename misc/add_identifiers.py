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

def is_blank(line):
    line = line.strip()
    return line=="*" or line==""


for root, dirs, files in os.walk("../casadi"):
   for name in files:
      if "runtime" in root: continue
      # "core" not in root and  "core" not in root : continue
      if name.endswith(".hpp"):
          file = os.path.join(root, name)
          with open(file,"r") as inp:
            data = inp.read()
            orig = data
            def process(match):
                global next_identifier
                if ("core" in root and ("\\brief" in match.group(0))) or "\\defgroup main_" in match.group(0) or "\\defgroup plugin_" in match.group(0):
                 if "\\identifier" in match.group(0):
                    identifier = str2id(re.search(r"\\identifier\{(.*?)\}",match.group(0)).group(1))
                    line = data[:match.start()].count("\n")+1
                    loc = file + ":" + str(line) + "\n" + match.group(0)
                    if identifier in identifiers:
                        raise Exception(f"Duplicate identifiers:\n{loc}\n{location[identifier]}")
                    identifiers.add(identifier)
                    location[identifier] = loc
                    
                    # Make sure there is a blank line after \brief when multi-line
                    # If there is no such newline, doxygen xml lumps brief and subsequent lines into a single para
                    if "\\brief" in match.group(0):
                        if "\n" in match.group(0):
                            split = match.group(0).split("\n")
                            if not is_blank(split[1]):
                                split = [split[0]] + [""]+split[1:]
                                return "\n".join(split)
                    # Make sure a defgroup starts with \par (needed for doxygeb 1.8.17)
                    if "defgroup main_" in match.group(0) or "defgroup plugin_" in match.group(0):
                        split = match.group(0).split("\n")
                        if "\\par" not in split[1]:
                            split = [split[0]] + [(" "*match.group(0).index("\\defgroup"))+"\\par"]+split[1:]
                            return "\n".join(split)
                        
                    # Make sure there is a blank line before \identifier
                    split = match.group(0).split("\n")
                    identifier_index = [i for i,e in enumerate(split) if e.strip().startswith("\\identifier")][0]
                    if not is_blank(split[identifier_index-1]):
                        split = split[:identifier_index] + [""] + split[identifier_index:]
                        return "\n".join(split)
                 else:
                    m = match.group(0).rstrip()
                    if "\\brief" in m:
                        ident = m.find("\\brief")*" "
                    else:
                        ident = m.find("\\defgroup")*" "
                    next_identifier+=1
                    return m[:-2].rstrip() + "\n\n" + ident+ "\\identifier{" + id2str(next_identifier) + "}" + " */"
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
