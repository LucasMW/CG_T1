
#t1: Cluster.o Image.o ImageLoader.o Pixel.o main.o
#	g++ -o sp Cluster.o Image.o ImageLoader.o main.o Pixel.o
sp1: 
	g++ -o sp Cluster.cpp Image.cpp ImageLoader.cpp main.cpp Pixel.cpp
ssp: 
	g++ -o ssp Cluster.cpp Image.cpp ImageLoader.cpp main.cpp Pixel.cpp -O2
Cluster.o:
	g++ -c Cluster.cpp
Image.o:
	g++ -c Image.cpp
ImageLoader.o:
	g++ -c ImageLoader.cpp
Pixel.o:
	g++ -c Pixel.cpp
main.o:
	g++ -c main.cpp
clean:
	rm -f *.o
	rm -f AB_ufv_06752.bmp 
	rm -f sp
	rm -f ssp