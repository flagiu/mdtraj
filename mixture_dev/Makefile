CXX = g++ -O3

mdtraj: src/mdtraj.cpp
	cd src && bash collect_all_headers.sh && cd .. && \
	$(CXX) src/mdtraj.cpp -o bin/mdtraj

clean:
	rm -f core bin/mdtraj
