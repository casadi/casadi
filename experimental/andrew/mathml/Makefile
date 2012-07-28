.DEFAULT:examples.m
example1.xml:examples.m
	matlab -nodesktop -nojvm -nosplash -r 'run examples; exit'
	xmllint --format $@ >> temp_pretty_xml_file
	mv temp_pretty_xml_file $@
