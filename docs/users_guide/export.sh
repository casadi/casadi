
rm -rf build && rm -f *.in && make bare_body
sed -i 's/{{</{\&lbrace;</g' build/singlehtml/index.html
sed -i 's/>}}/>}\&rbrace;/g' build/singlehtml/index.html
cp build/singlehtml/index.html ~/work/casadi.github.io/content/docs/_index.html
make bare_toc
sed -i 's/index\.html//g' build/singlehtml/index.html
cat addendum.html build/singlehtml/index.html > /home/jgillis/work/casadi.github.io/themes/casadi-theme/layouts/partials/docs-sidebar.html

