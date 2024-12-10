CXX = g++ -O3 -std=c++17

all: mdtraj diffusion

mdtraj: src/mdtraj.cpp
	cd src && bash collect_all_headers.sh && cd .. && \
	$(CXX) src/mdtraj.cpp -o bin/mdtraj

diffusion: src/lib/diffusion.cpp
	$(CXX) src/lib/diffusion.cpp -o bin/difflux1d

clean:
	rm -f core bin/mdtraj bin/difflux1d
