import os
import sys
import struct
import shutil
import hashlib
import subprocess

# Authenticode-sign the unsigned Windows binaries in `dir` with Azure Artifact Signing
# The access token is read from the AZURE_SIGNING_TOKEN environment variable
# If `like` is given, only binaries that also occur in that directory are signed

dir = sys.argv[1]
jsign = sys.argv[2]
endpoint = sys.argv[3]
alias = sys.argv[4]
like = sys.argv[5] if len(sys.argv)>5 else None

def is_signed(f):
    f.seek(0x3c)
    pe = struct.unpack("<I",f.read(4))[0]
    f.seek(pe+24)
    magic = struct.unpack("<H",f.read(2))[0]
    # Certificate table entry of the data directories
    f.seek(pe+24+(112 if magic==0x20b else 96)+4*8)
    return struct.unpack("<II",f.read(8))[1]>0

# Unsigned PE files (.dll, .exe, .pyd, .mexw64, .mex, .oct), grouped by content
def signables(dir):
    ret = {}
    for root, dirs, files in os.walk(dir):
        for file in files:
            path = os.path.join(root, file)
            if os.path.islink(path): continue
            with open(path, "rb") as f:
                if f.read(2)!=b"MZ" or is_signed(f): continue
                f.seek(0)
                ret.setdefault(hashlib.sha256(f.read()).hexdigest(),[]).append(path)
    return ret

groups = signables(dir)
if like is not None:
    keep = signables(like)
    groups = dict((k,v) for k,v in groups.items() if k in keep)
if len(groups)==0:
    print("Nothing to sign")
    sys.exit(0)

# Keep modification times: the build tree is reused
mtimes = dict((path,os.stat(path).st_mtime) for paths in groups.values() for path in paths)

# Single invocation: every invocation costs an extra request
filelist = "sign_windows.txt"
with open(filelist,"w") as f:
    for paths in groups.values():
        print("Signing",paths[0])
        f.write(paths[0]+"\n")

assert subprocess.call(["java","-jar",jsign,"--storetype","TRUSTEDSIGNING","--keystore",endpoint,"--storepass","env:AZURE_SIGNING_TOKEN","--alias",alias,"@"+filelist])==0

os.remove(filelist)

for paths in groups.values():
    for path in paths[1:]:
        print("Copying to",path)
        shutil.copyfile(paths[0],path)
    for path in paths:
        os.utime(path,(mtimes[path],mtimes[path]))
