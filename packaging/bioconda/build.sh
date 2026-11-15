#!/bin/bash
set -euo pipefail
mkdir -p "${PREFIX}/bin"
${CC} ${CFLAGS} ${CPPFLAGS} -O3 -o tandex src/tandex.c ${LDFLAGS} -lm -lz
install -m 755 tandex "${PREFIX}/bin/"
