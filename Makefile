hybridizer: hybridizer.o read_data.o
	gcc -g -Wall hybridizer.o read_data.o -o hybridizer `pkg-config --cflags --libs glib-2.0` `pkg-config --cflags --libs gsl`

hybridizer.o: hybridizer.c
	gcc -c -g -Wall hybridizer.c `pkg-config --cflags --libs glib-2.0` `pkg-config --cflags --libs gsl`

read_data.o: read_data.c read_data.h
	gcc -c -g -Wall read_data.c `pkg-config --cflags --libs glib-2.0` `pkg-config --cflags --libs gsl`
