#!/bin/bash

function fileupdate {
myfile=$1

TMP="/tmp/FC.tmp"

cat > $TMP

DIFF=`diff -q $TMP $myfile`

if [ ! "$DIFF" = "" ]; then
  if [ -s $TMP ]; then
	  if [ "$2" = "" ]; then
		  mv $myfile  $myfile~
	  elif [ ! "$2" = "none" ]; then
		  mv $myfile  $myfile$2
	  fi
	  mv $TMP $myfile
	  echo "Changed $1";
	fi
fi
}

function findreplace {
  echo "This is findrep"
  echo "perl regex: ${FGRED}$3${NORMAL}"
  echo "(like it would appear in perl file)"
  for file in `find -name "$1" -type f | grep -v .svn | grep -v .git | grep -v .findrep~`;
    do perl -p -e "chomp;my \$in = \$_;my \$out=\$_; if (\$out=~$2) {\$_ = qq<${FGMAGENTA}$file${NORMAL}:${FGGREEN}\$.${FGCYAN}:\n${NORMAL}\$in\n${FGRED}\$out${NORMAL}\n> } else {\$_=''}" $file
    perl -p -w -e "$2" $file | fileupdate $file .findrep~ ;
  done;
}
           
findreplace "*.py" "s/(?<method>\.(s|g)et(Input|Output|FwdSeed|AdjSeed|FwdSens|AdjSens)\s*)\((?<bulk>.*)?,(?<white>\s*)(DAE|INTEGRATOR|NLP_SOLVER|QP|RDAE|ACADO|MAYER|OCP|MUSCOD|CONTROL_DAE|CONTROLSIMULATOR|SDP)_(?<scheme>[A-Z_0]*\s*)\)/qq|\$+{method}(\$+{bulk},\$+{white}\"|.lc(\$+{scheme}).qq|\")|/ge"

findreplace "*.py" "s/(?<method>\.(input|output|fwdSeed|adjSeed|fwdSens|adjSens)\s*)\((?<white>\s*)(DAE|INTEGRATOR|NLP_SOLVER|QP|RDAE|ACADO|MAYER|OCP|MUSCOD|CONTROL_DAE|CONTROLSIMULATOR|SDP)_(?<scheme>[A-Z_0]*\s*)\)/qq|\$+{method}(\$+{white}\"|.lc(\$+{scheme}).qq|\")|/ge"

findreplace "*.py" "s/(?<method>\.(jacobian|hessian|gradient|jac|hess|grad)\s*)\((?<white>\s*)(DAE|INTEGRATOR|NLP_SOLVER|QP|RDAE|ACADO|MAYER|OCP|MUSCOD|CONTROL_DAE|CONTROLSIMULATOR|SDP)_(?<scheme>[A-Z_0]*\s*),(?<white2>\s*)(DAE|INTEGRATOR|NLP_SOLVER|QP|RDAE|ACADO|MAYER|OCP|MUSCOD|CONTROL_DAE|CONTROLSIMULATOR|SDP)_(?<scheme2>[A-Z_0]*\s*)\)/qq|\$+{method}(\$+{white}\"|.lc(\$+{scheme}).qq|\",\$+{white2}\"| .lc(\$+{scheme2}) . qq|\")|/ge"
