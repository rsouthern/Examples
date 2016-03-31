This is a general dimensional Moving Least Squares projection library which I developed for use in another package. 
It is essentially a header only C++ library consisting of the single file pointcloud.hpp. 
Check out main.cpp for an example of usage.

This is still very much in development, and I will write a blog entry to explain what is going on shortly.

It makes use of nanoflann, a nearest neighbour library using KD trees to improve performance (included). 
It also makes heavy use of Eigen3 and of boost (although this dependency can be trivially removed).
