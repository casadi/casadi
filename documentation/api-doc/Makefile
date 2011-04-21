all:
	doxygen Doxyfile.in

full:
	doxygen Doxyfile.in
	python generate_options_doc.py
	python fetch_examples_pdf.py
	doxygen Doxyfile.in

publish: 
	./publish.sh
	
clean: 
	rm -rf XML/ html/ 
