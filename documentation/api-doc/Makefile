all:
	doxygen Doxyfile.in

full:
	doxygen Doxyfile.in
	python extra/generate_options_doc.py
	python extra/fetch_examples_pdf.py
	doxygen Doxyfile.in

publish: 
	./publish.sh
	
clean: 
	rm -rf XML/ html/ 
