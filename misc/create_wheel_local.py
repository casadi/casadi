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
import json
import fnmatch

def open_for_csv(name, mode):
    if sys.version_info[0] < 3:
        kwargs = {}
        mode += 'b'
    else:
        kwargs = {'newline': '', 'encoding': 'utf-8'}

    return open(name, mode, **kwargs)

# Optional --manifest <path> trailing argument
manifest_path = None
argv = list(sys.argv)
if "--manifest" in argv:
  i = argv.index("--manifest")
  manifest_path = argv[i+1]
  del argv[i:i+2]

version   = argv[1]
pyversion = argv[2]
os_name   = argv[3]
bitness   = argv[4]
arch      = argv[5]
dir_name  = argv[6]


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
  elif arch=="manylinux2014-x64-modern":
    arch = "manylinux2014_x86_64"
  elif arch=="manylinux2014-x86-modern":
    arch = "manylinux2014_i686"
  tag = "cp%s-none-%s" % (pyversion,arch.replace("-","_"))
elif os_name=="osx":
  if arch=="x86_64":
    tag = ["cp%s-none-macosx_11_0_x86_64" % (pyversion),
         "cp%s-none-macosx_11_0_intel" % (pyversion)]
  elif arch=="arm64":
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

def normalize_dist(name):
  # PEP 427: replace runs of non-alphanumerics with a single underscore
  return re.sub(r"[-_.]+", "_", name)

def fullname_for(dist_name):
  base = dist_name + "-" + version
  if isinstance(tag, list):
    fn = base + "-" + tag[0]
    for t in tag[1:]:
      fn += "." + t.split("-")[-1]
    return fn
  return base + "-" + tag

def write_wheel_metafiles(bdist_dir, dist_name, summary, description,
                          requires_dist=None, provides_extra=None):
  wheel_dist = dist_name + "-" + version
  distinfo_dir = os.path.join(bdist_dir, '%s.dist-info' % wheel_dist)
  if not os.path.exists(distinfo_dir):
    os.mkdir(distinfo_dir)

  msg = Message()
  msg['Wheel-Version'] = '1.0'
  msg['Root-Is-Purelib'] = "false"
  if isinstance(tag, list):
    for t in tag: msg["Tag"] = t
  else:
    msg["Tag"] = tag
  with open(os.path.join(distinfo_dir, 'WHEEL'), 'w') as f:
    Generator(f, maxheaderlen=0).flatten(msg)

  lines = [
    "Metadata-Version: 2.1",
    "Name: " + dist_name,
    "Version: " + version,
    "Summary: " + summary,
    "Home-page: http://casadi.org",
    "Author: Joel Andersson, Joris Gillis, Greg Horn",
    "Author-email:  developer_first_name@casadi.org",
    "License: GNU Lesser General Public License v3 or later (LGPLv3+)",
    "Download-URL: http://install.casadi.org",
    "Project-URL: Documentation, http://docs.casadi.org",
    "Project-URL: Bug Tracker, https://github.com/casadi/casadi/issues",
    "Project-URL: Source Code, https://github.com/casadi/casadi",
    "Platform: Windows",
    "Platform: Linux",
    "Platform: Mac OS-X",
    "Classifier: Development Status :: 5 - Production/Stable",
    "Classifier: Intended Audience :: Science/Research",
    "Classifier: License :: OSI Approved",
    "Classifier: Programming Language :: C++",
    "Classifier: Programming Language :: Python",
    "Classifier: Programming Language :: Python :: 3",
    "Classifier: Programming Language :: Python :: Implementation :: CPython",
    "Classifier: Topic :: Scientific/Engineering",
    "Classifier: Operating System :: Microsoft :: Windows",
    "Classifier: Operating System :: POSIX",
    "Classifier: Operating System :: Unix",
    "Classifier: Operating System :: MacOS",
  ]
  for r in (requires_dist or []):
    lines.append("Requires-Dist: " + r)
  for e in (provides_extra or []):
    lines.append("Provides-Extra: " + e)
  lines.append("")
  lines.append(description)
  lines.append("")
  with open(os.path.join(distinfo_dir, 'METADATA'), 'w') as f:
    f.write("\n".join(lines))

  return distinfo_dir

# Consolidate non-casadi-prefixed top-level dirs into casadi/ (existing behavior)
bdist_dir = dir_name
for d in os.listdir(dir_name):
  if not d.startswith("casadi") and os.path.isdir(os.path.join(dir_name, d)):
    shutil.copytree(os.path.join(dir_name, d), os.path.join(bdist_dir, "casadi", d), dirs_exist_ok=True)
    shutil.rmtree(os.path.join(dir_name, d))

# Parse manifest
plugins = []  # list of dicts: {name, dist_name, normalized, patterns, summary, extras}
if manifest_path:
  with open(manifest_path) as f:
    raw = json.load(f)
  for name, spec in raw.items():
    patterns = spec.get("patterns") or spec.get("match") or []
    if isinstance(patterns, str):
      patterns = [patterns]
    plugins.append({
      "dist_name": name,
      "normalized": normalize_dist(name),
      "patterns": patterns,
      "summary": spec.get("summary", "CasADi plugin: " + name),
      "extras": spec.get("extras"),  # extras key for main wheel; defaults below
    })

casadi_pkg = os.path.join(bdist_dir, "casadi")

# For each plugin, carve out matching files into a separate staging dir and build a wheel.
# Matching is by file basename via fnmatch (case-insensitive on Windows-style; we use plain fnmatch).
plugin_assignments = {}  # absolute file path -> plugin index (first match wins)
for idx, plug in enumerate(plugins):
  for root, _dirs, files in os.walk(casadi_pkg):
    for fname in files:
      if any(fnmatch.fnmatchcase(fname, pat) for pat in plug["patterns"]):
        absp = os.path.join(root, fname)
        plugin_assignments.setdefault(absp, idx)

built_wheels = []
plugin_stage_dirs = []

for idx, plug in enumerate(plugins):
  stage = os.path.join(os.path.dirname(os.path.abspath(bdist_dir)),
                       "_plugin_stage_" + plug["normalized"])
  if os.path.exists(stage):
    shutil.rmtree(stage)
  os.makedirs(os.path.join(stage, "casadi"))
  plugin_stage_dirs.append(stage)

  moved_any = False
  for absp, owner in list(plugin_assignments.items()):
    if owner != idx:
      continue
    rel = os.path.relpath(absp, casadi_pkg)
    dst = os.path.join(stage, "casadi", rel)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.move(absp, dst)
    moved_any = True

  if not moved_any:
    sys.stderr.write("[create_wheel] WARNING: plugin %r matched no files\n" % plug["dist_name"])

  desc = "Plugin package for CasADi: %s. Install via casadi[%s]." % (
    plug["dist_name"], plug["extras"] or plug["dist_name"].replace("casadi-plugin-", ""))
  write_wheel_metafiles(
    stage,
    plug["normalized"],
    plug["summary"],
    desc,
    requires_dist=["casadi==" + version.split("+")[0]],
  )
  fullname = fullname_for(plug["normalized"])
  distinfo_dir = os.path.join(stage, plug["normalized"] + "-" + version + ".dist-info")
  write_record(stage, distinfo_dir)
  archive_wheelfile(os.path.join(os.path.dirname(stage), fullname), stage)
  built_wheels.append(os.path.join(os.path.dirname(stage), fullname) + ".whl")

# Now build the main casadi wheel from whatever remains
extras_lines = []
for plug in plugins:
  extra = plug["extras"] or plug["dist_name"].replace("casadi-plugin-", "")
  extras_lines.append((extra, plug["dist_name"]))

requires_dist = ["numpy"] + [
  '%s ; extra == "%s"' % (dist, extra) for extra, dist in extras_lines
]
provides_extra = [extra for extra, _ in extras_lines]

main_description = (
  "CasADi is a symbolic framework for numeric optimization implementing automatic "
  "differentiation in forward and reverse modes on sparse matrix-valued computational "
  "graphs. It supports self-contained C-code generation and interfaces state-of-the-art "
  "codes such as SUNDIALS, IPOPT etc. It can be used from C++, Python or Matlab/Octave.\n\n"
  "Example pack can be downloaded from http://install.casadi.org"
)

write_wheel_metafiles(
  bdist_dir,
  "casadi",
  "CasADi -- framework for algorithmic differentiation and numeric optimization",
  main_description,
  requires_dist=requires_dist,
  provides_extra=provides_extra,
)

main_distinfo = os.path.join(bdist_dir, "casadi-" + version + ".dist-info")
main_fullname = fullname_for("casadi")
write_record(bdist_dir, main_distinfo)
archive_wheelfile(main_fullname, dir_name)
built_wheels.append(main_fullname + ".whl")

if manifest_path:
  # multi-wheel mode: print one wheel per line
  sys.stdout.write("\n".join(built_wheels))
else:
  # back-compat: single wheel filename, no trailing newline
  sys.stdout.write(main_fullname + ".whl")
