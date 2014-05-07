# Make NBODY6
cd src
make
rm *.o
mv nbody6 ../

# Make MCluster
gcc -O2 -fopenmp -ffast-math -I /usr/local/Cellar/gcc47/4.7.3/lib/gcc/x86_64-apple-darwin13.1.0/4.7.3/include/ -o mcluster main.c -lm
mv mcluster ../
