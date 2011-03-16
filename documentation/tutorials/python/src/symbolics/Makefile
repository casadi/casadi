SRCS = \
	sx.py \
	sxfunction.py\
	mx.py\
	integration.py



PDFS = $(SRCS:.py=.pdf)

all: ${PDFS}

%.pdf : %.py
	@echo "Creating pyreport" $@
	pyreport -d -l $< 
