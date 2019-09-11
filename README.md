# SimpleMVPTree
An efficient data structure for distance queries for high-dimensional data in metric spaces.

My implementation is based on the "two-level mvp-tree that keeps all vantage points in a single directory", which performed best in the experiments in the paper referenced below.

# Basic Usage
```c++
#include "SimpleMVPTree.hpp"

class MyClass {
// ... implement MyClass.
};

// Any distance function that is a metric is fine here.
double myDistance(const MyClass& a, const MyClass& b) {
// ... implement the distance function.
}

void demonstrateUsage() {
  std::vector<MyClass> myData;
  // ... fill myData with elements.
  SimpleMVPTree<MyClass> mvp(myData, myDistance);
  
  MyClass query;
  // ... define the query object.
  
  // Use Case 1: Find the closest distance to the query object
  double closestDistance = mvp.search_mindist(query);

  // Use Case 2: Find the closest object to the query object
  MyClass closestObject = mvp.search_mindist_elem(query);
}
```

# References

Bozkaya, Tolga; Ozsoyoglu, Meral 1999."Indexing Large Metric Spaces for Similarity Search Queries". ACM Transactions in Database Systems, Vol. 24, No. 3, September 1999, pg. 361-404.
