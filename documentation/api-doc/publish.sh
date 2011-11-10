#!/bin/bash
echo "This script will synchronize the html doc/folder with http://casadi.sourceforge.net/api/html/"
echo -n "Please type your sourceforge username:"
read username
rsync -avP -e ssh html "$username,casadi@web.sourceforge.net:/home/groups/c/ca/casadi/htdocs/api/"
rsync -avP -e ssh ../tutorials/python/pdf/ "$username,casadi@web.sourceforge.net:/home/groups/c/ca/casadi/htdocs/tutorials/"
rsync -avP -e ssh ../cheatsheet/*.pdf "$username,casadi@web.sourceforge.net:/home/groups/c/ca/casadi/htdocs/cheatsheets/"
#rsync -avP -e ssh sphinx/_build/html "$username,casadi@web.sourceforge.net:/home/groups/c/ca/casadi/htdocs/api/python/"
