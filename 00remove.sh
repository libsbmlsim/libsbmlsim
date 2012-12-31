#!/bin/zsh
# Run this script before making a release tar ball.
# ./00remove.sh
# rm 00remove.sh
rm -rf **/CVS
rm -rf bindings diff simulation_results
rm -f l?v*
rm -f header.txt
rm -f src/ev_func_piece.xml
rm -f src/Makefile.dist
rm -f src/bindings/java/build_java.sh
rm -rf src/bindings/java/bin
rm -f src/bindings/python/build_python.sh
rm -f **/.DS_Store(N)
