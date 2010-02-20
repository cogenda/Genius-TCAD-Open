All:contrib

SUBDIRS := ${shell ls -F | grep /}

#SUBDIRS := ann brkpnts ctri csabathy expreval dfise metis \
#           qshep silvaco spice spline tif triangle tetgen

contrib:
	@for dir in $(SUBDIRS) ; do \
		(cd $$dir ; make) || exit 1 ; \
	done


.PHONY : clean distclean

clean:
	@for dir in $(SUBDIRS) ; do \
		(cd $$dir ; make clean) || exit 1 ; \
	done



distclean:
	@for dir in $(SUBDIRS) ; do \
		(cd $$dir ; make distclean) || exit 1 ; \
	done


