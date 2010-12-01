dot -Tpng < graph.dot > overview.png
echo "<html>" > overview.html
dot -Tcmapx < graph.dot >> overview.html
echo "<IMG SRC='overview.png' USEMAP='#overview' /></html>" >> overview.html
