#!/bin/bash
echo "This script will synchronize the html doc/folder with http://casadi.sourceforge.net/api/html/"
if [ -z "$sfaccount" ]; then
echo -n "Please type your sourceforge username:"
read username
else
username="$sfaccount"
fi
rsync -avP -e ssh html "$username,casadi@web.sourceforge.net:/home/groups/c/ca/casadi/htdocs/api"
rsync -avP -e ssh ../tutorials/python/pdf/ "$username,casadi@web.sourceforge.net:/home/groups/c/ca/casadi/htdocs/tutorials/"
rsync -avP -e ssh ../documents/*.pdf "$username,casadi@web.sourceforge.net:/home/groups/c/ca/casadi/htdocs/documents/"
#rsync -avP -e ssh ../tutorials/cpp/ "$username,casadi@web.sourceforge.net:/home/groups/c/ca/casadi/htdocs/tutorials/cpp"
rsync -avP -e ssh ../cheatsheet/*.pdf "$username,casadi@web.sourceforge.net:/home/groups/c/ca/casadi/htdocs/cheatsheets/"
#rsync -avP -e ssh sphinx/_build/html "$username,casadi@web.sourceforge.net:/home/groups/c/ca/casadi/htdocs/api/python/"
rsync -avP -e ssh ../users_guide/*.pdf "$username,casadi@web.sourceforge.net:/home/groups/c/ca/casadi/htdocs/users_guide"
rsync -avP -e ssh ../users_guide/casadi-users_guide/* "$username,casadi@web.sourceforge.net:/home/groups/c/ca/casadi/htdocs/users_guide/html"


echo -e "cd tested\nput ../example_pack/example_pack.zip" | sftp -- casaditestbot,casadi@web.sourceforge.net:/home/pfs/project/c/ca/casadi/CasADi
