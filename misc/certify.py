import os
import sys
import subprocess
import shutil

identity = sys.argv[1]
archive = sys.argv[2]
apple_id = sys.argv[3]
apple_team_id = sys.argv[4]
notarize_password = sys.argv[5]

# Unzip into dir
dir = "certify_dir"
notarize_dir = "notarize_dir"
shutil.rmtree(notarize_dir,ignore_errors=True)
os.mkdir(notarize_dir)
shutil.rmtree(dir,ignore_errors=True)
assert subprocess.call(["unzip",archive,"-d",dir])==0


def each_signable(dir):
    for root, dirs, files in os.walk(dir):
        for file in files:
            if file.endswith(".dylib") or file.endswith(".mexmaci64") or file.endswith(".mexmaca64") or file.endswith(".so") or file.endswith(".mex") or os.access(file, os.X_OK):
                yield os.path.join(root, file)
            

# Recursively look for all shared libraries in `dir`
for path in each_signable(dir):
    print("Signing",path)
    assert subprocess.call(["codesign", "--remove-signature", path])==0
    assert subprocess.call(["codesign", "--force", "--sign", identity, path])==0
    shutil.copy(path,notarize_dir)


assert subprocess.call(["zip","-rq","../" + "notarize_" + archive, "."],cwd=notarize_dir)==0

assert subprocess.call(["xcrun","notarytool","submit","notarize_" + archive, "--apple-id", apple_id, "--team-id", apple_team_id, "--password", notarize_password, "--wait"])==0
for path in each_signable(dir):
    assert subprocess.call(["xcrun", "stapler", "staple", path])==0

# Zip again
os.remove(archive)
subprocess.call(["zip","-rq","../" + archive, "."],cwd=dir)

shutil.rmtree(dir)


