CXX = g++ -O3 -std=c++17

mdtraj: src/mdtraj.cpp
	cd src && bash collect_all_headers.sh && cd .. && \
	$(CXX) src/mdtraj.cpp -o bin/mdtraj

clean:
	rm -f core bin/mdtraj
