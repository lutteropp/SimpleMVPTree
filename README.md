# SimpleMVPTree
An efficient data structure for nearest neighbor queries for high-dimensional data in metric spaces.

My implementation is based on the "two-level mvp-tree that keeps all vantage points in a single directory", which performed best in the experiments in the paper referenced below.

Using the SimpleMVPTree in your C++ project is very simple. Just include the `SimpleMVPTree.hpp` header.

# Please Note
The distance function `myDistance(a, b)` provided by you must be a metric! This means, the following conditions must hold:
1. `myDistance(a, b) >= 0` and `myDistance(a,b) == 0, if and only if a == b` for all a, b
2. `myDistance(a, b) == myDistance(b, a)` for all a, b
3. `myDistance(a, c) <= myDistance(a, b) + myDistance(b, c)` for all a, b, c

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
  
  // build the data structure
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
