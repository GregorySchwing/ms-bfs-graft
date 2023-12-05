# README #
Multithreaded code for computing maximum cardinality matching.
Author: Ariful Azad, Aydin Buluç, Lawrence Berkeley National Laboratory

References: 


1. "Computing Maximum Cardinality Matchings in Parallel on Bipartite Graphs via Tree-Grafting", A. Azad, A. Buluç, A. Pothen, IEEE Transactions on Parallel and Distributed Systems (TPDS), 2016, DOI: 10.1109/TPDS.2016.2546258.

2. "A Parallel Tree Grafting Algorithm for Maximum Cardinality Matching in Bipartite Graphs", A. Azad, A. Buluç, A. Pothen, IEEE International Parallel & Distributed Processing Symposium (IPDPS), 2015.




*** Copyright Notice ***

"Parallel Maximum Cardinality Matchings via Tree-Grafting" Copyright (c) 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy).  All rights reserved.

If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Innovation & Partnerships Office at  IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department of Energy and the U.S. Government consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, prepare derivative works, and perform publicly and display publicly, and to permit other to do so.


To compile: 

```
#!c++

make all
```

**To run:**

```
#!c++

./msBFSGraft file_name (matrix market format) nthreads
```

**Example serial run (assume the input file amazon0312.mtx is stored in the current directory):**

```
#!c++

./msBFSGraft amazon0312.mtx 1
```

**Example parallel run with 4 threads (assume the input file amazon0312.mtx is stored in the current directory):**

```
#!c++

./msBFSGraft amazon0312.mtx 4
```