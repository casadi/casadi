echo "Please enter sf username"
read username
rsync -avP -e ssh pdf/ "$username,casadi@web.sourceforge.net:/home/groups/c/ca/casadi/htdocs/tutorials/"
