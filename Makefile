CC      ?= gcc
CFLAGS  ?= -O3 -Wall
PREFIX  ?= /usr/local

tandex: src/tandex.c
	$(CC) $(CFLAGS) -o $@ $< -lm -lz

test: tandex
	python3 tests/run_tests.py ./tandex

install: tandex
	install -d $(DESTDIR)$(PREFIX)/bin
	install -m 755 tandex $(DESTDIR)$(PREFIX)/bin/

clean:
	rm -f tandex

.PHONY: test install clean
