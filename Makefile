.PHONY: all clean install uninstall

all:
	omake

clean:
	omake clean

install:
	cd lib && omake install

uninstall:
	cd lib && omake uninstall
