SABRModel:LcgGenerator.o BoxMullerGenerator.o BSModel.o BSpline.o SABR.o Swaption.o execution.cpp
	g++ -Wall -std=c++0x LcgGenerator.o BoxMullerGenerator.o BSModel.o BSpline.o SABR.o Swaption.o execution.cpp -o SABRModel

LcgGenerator.o:LcgGenerator.cpp LcgGenerator.hpp
        g++ -Wall -std=c++0x LcgGenerator.cpp -o LcgGenerator.o

BoxMullerGenerator.o:BoxMullerGenerator.hpp BoxMullerGenerator.cpp
        g++ -Wall -c -std=c++0x BoxMullerGenerator.cpp -o BoxMullerGenerator.o

SABR.o:SABR.h SABR.cpp BSModel.h
	g++ -Wall -c -std=c++0x SABR.cpp -o SABR.o	

BSModel.o:BSModel.cpp BSModel.h
	g++ -Wall -c -std=c++0x BSModel.cpp -o BSModel.o

BSpline.o:BSpline.cpp BSpline.h
	g++ -Wall -c -std=c++0x BSpline.cpp -o BSpline.o

Swaption.o:Swaption.cpp Swaption.h BSModel.h
	g++ -Wall -c -std=c++0x Swaption.cpp -o Swaption.o

clean:
	rm -rf *o SABRModel


       
