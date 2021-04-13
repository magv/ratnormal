XCXXFLAGS=${CXXFLAGS} \
	  -std=c++14 -Os \
	  -pedantic -Wall -Wno-deprecated-declarations -Wno-format

XLDFLAGS=${LDFLAGS} -lflint -lgmp -lginac -lcln

ratnormal: ratnormal.cpp Makefile
	${CXX} ${XCXXFLAGS} -o $@ ratnormal.cpp ${XLDFLAGS}

README.md: ratnormal.cpp mkmanual.sh
	sed '/MANUAL/{n;q}' $@ >$@.tmp
	./mkmanual.sh >>$@.tmp <$<
	mv $@.tmp $@
