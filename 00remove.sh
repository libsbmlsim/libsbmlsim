# run this script before making a release tar ball.
# ./00remove.sh
# rm 00remove.sh
rm -rf **/CVS
rm -f **/.DS_*
rm -rf bindings diff simulation_results
rm -f l?v*
rm -f src/ev_func_piece.xml
rm -f src/Makefile.dist
rm -f header.txt
