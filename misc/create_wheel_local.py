import os
import base64
from wheel.archive import archive_wheelfile
import hashlib
import csv
from email.message import Message
from email.generator import Generator
import shutil
import re
import sys

def open_for_csv(name, mode):
    if sys.version_info[0] < 3:
        kwargs = {}
        mode += 'b'
    else:
        kwargs = {'newline': '', 'encoding': 'utf-8'}

    return open(name, mode, **kwargs)

version   = sys.argv[1]
pyversion = sys.argv[2]
os_name   = sys.argv[3]
bitness   = sys.argv[4]
arch      = sys.argv[5]
dir_name  = sys.argv[6]


if "+" in version:
  [pre,post] = version.split("+")
  post = re.sub(r"[^a-zA-Z0-9]",".",post)
  version=pre+"+"+post

if version.startswith("v"):
  version = version[1:]
if os_name=="linux":
  if arch=="manylinux1-x86":
    arch = "manylinux1_i686"
  elif arch=="manylinux1-x64":
    arch = "manylinux2010_x86_64"
  elif arch=="manylinux2014-x64":
    arch = "manylinux2014_x86_64"
  elif arch=="manylinux2014-x86":
    arch = "manylinux2014_i686"
  tag = "cp%s-none-%s" % (pyversion,arch.replace("-","_"))
elif os_name=="osx":
  if arch=="osx":
    tag = ["cp%s-none-macosx_10_13_x86_64" % (pyversion),
         "cp%s-none-macosx_10_13_intel" % (pyversion)]
  elif arch=="osx-m1":
    tag = "cp%s-none-macosx_11_0_arm64" % (pyversion)
elif os_name=="windows":
  if bitness=="64":
    tag = "cp%s-none-win_amd64" % pyversion
  else:
    tag = "cp%s-none-win32" % pyversion
else:
  raise Exception()

def write_record(bdist_dir, distinfo_dir):        

      record_path = os.path.join(distinfo_dir, 'RECORD')
      record_relpath = os.path.relpath(record_path, bdist_dir)

      def walk():
          for dir, dirs, files in os.walk(bdist_dir):
              dirs.sort()
              for f in sorted(files):
                  yield os.path.join(dir, f)

      def skip(path):
          """Wheel hashes every possible file."""
          return (path == record_relpath)

      with open_for_csv(record_path, 'w+') as record_file:
          writer = csv.writer(record_file)
          for path in walk():
              relpath = os.path.relpath(path, bdist_dir)
              if skip(relpath):
                  hash = ''
                  size = ''
              else:
                  with open(path, 'rb') as f:
                      data = f.read()
                  digest = hashlib.sha256(data).digest()
                  hash = 'sha256=' + base64.urlsafe_b64encode(digest).decode("ascii").strip("=")
                  size = len(data)
              record_path = os.path.relpath(
                  path, bdist_dir).replace(os.path.sep, '/')
              writer.writerow((record_path, hash, size))
             
wheel_dist_name = "casadi" + "-" + version
bdist_dir = dir_name

distinfo_dir = os.path.join(bdist_dir,'%s.dist-info' % wheel_dist_name)
if not os.path.exists(distinfo_dir):
  os.mkdir(distinfo_dir)

for d in os.listdir(dir_name):
  if not d.startswith("casadi") and os.path.isdir(os.path.join(dir_name, d)):
    shutil.copytree(os.path.join(dir_name, d),os.path.join(bdist_dir,"casadi",d),dirs_exist_ok=True)
    shutil.rmtree(os.path.join(dir_name, d))

msg = Message()
msg['Wheel-Version'] = '1.0'  # of the spec
#msg['Generator'] = generator
msg['Root-Is-Purelib'] = "false"
if isinstance(tag,list):
  for t in tag: msg["Tag"] = t
else:
  msg["Tag"] = tag



wheelfile_path = os.path.join(distinfo_dir, 'WHEEL')
with open(wheelfile_path, 'w') as f:
    Generator(f, maxheaderlen=0).flatten(msg)

metadata_path = os.path.join(distinfo_dir, 'METADATA')
with open(metadata_path, 'w') as f:
    f.write("""Metadata-Version: 2.0
Name: casadi
Version: %s
Summary: CasADi -- framework for algorithmic differentiation and numeric optimization
Home-page: http://casadi.org
Author: Joel Andersson, Joris Gillis, Greg Horn
Author-email:  developer_first_name@casadi.org
License: GNU Lesser General Public License v3 or later (LGPLv3+)
Download-URL: http://install.casadi.org
Project-URL: https://github.com/casadi/casadi
Platform: Windows
Platform: Linux
Platform: Mac OS-X
Classifier: Development Status :: 5 - Production/Stable
Classifier: Intended Audience :: Science/Research
Classifier: License :: OSI Approved
Classifier: Programming Language :: C++
Classifier: Programming Language :: Python
Classifier: Programming Language :: Python :: 2
Classifier: Programming Language :: Python :: 2.7
Classifier: Programming Language :: Python :: 3
Classifier: Programming Language :: Python :: 3.4
Classifier: Programming Language :: Python :: 3.5
Classifier: Programming Language :: Python :: 3.6
Classifier: Programming Language :: Python :: 3.6
Classifier: Programming Language :: Python :: 3.7
Classifier: Programming Language :: Python :: 3.8
Classifier: Programming Language :: Python :: 3.9
Classifier: Programming Language :: Python :: 3.10
Classifier: Programming Language :: Python :: 3.11
Classifier: Programming Language :: Python :: 3.12
Classifier: Programming Language :: Python :: Implementation :: CPython
Classifier: Topic :: Scientific/Engineering
Classifier: Operating System :: Microsoft :: Windows
Classifier: Operating System :: POSIX
Classifier: Operating System :: Unix
Classifier: Operating System :: MacOS
Requires-Dist: numpy

CasADi is a symbolic framework for numeric optimization implementing automatic differentiation in forward and reverse modes on sparse matrix-valued computational graphs. It supports self-contained C-code generation and interfaces state-of-the-art codes such as SUNDIALS, IPOPT etc. It can be used from C++, Python or Matlab/Octave.

Example pack can be downloaded from http://install.casadi.org

""" % (version))


#wheel_dist_name = wheel_dist_name+".post1"
if isinstance(tag,list):
  fullname = wheel_dist_name+"-"+tag[0]
  for t in tag[1:]:
    fullname+="."+t.split("-")[-1]
else:
  fullname = wheel_dist_name+"-"+tag

write_record(bdist_dir, distinfo_dir)
archive_wheelfile(fullname,dir_name)

import sys
sys.stdout.write(fullname+".whl")
