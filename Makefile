XCXXFLAGS=${CXXFLAGS} \
	  -std=c++14 -Os \
	  -pedantic -Wall -Wno-deprecated-declarations -Wno-format

XCXXFLAGS_STATIC=${XCXXFLAGS} -Os -s -static -pthread

XLDFLAGS=${LDFLAGS} -lflint -lmpfr -lgmp -lginac -lcln

XLDFLAGS_STATIC=${XLDFLAGS}

ratnormal: ratnormal.cpp Makefile
	@date "+static const char VERSION[] = \"Ratnormal $$(git --git-dir=.git rev-parse --short=12 HEAD), built on %Y-%m-%d\n\";" >version.h
	${CXX} ${XCXXFLAGS} -include version.h -o $@ ratnormal.cpp ${XLDFLAGS}

ratnormal.static: ratnormal.cpp Makefile mkversion.sh
	env CXX="${CXX}" ./mkversion.sh >version.h
	${CXX} ${XCXXFLAGS_STATIC} -include version.h -o $@ ratnormal.cpp ${XLDFLAGS_STATIC}
	@upx --best "$@"

README.md: ratnormal.cpp mkmanual.sh
	sed '/MANUAL/{n;q}' $@ >$@.tmp
	./mkmanual.sh >>$@.tmp <$<
	mv $@.tmp $@
