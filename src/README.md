### Installation

To install the C++ version of cdc in Linux or Mac OS X you will need a compiler supporting C++11 (i.e. gcc >= 4.7 or Clang >= 3.0) and Cmake. To build start a terminal from the ranger main directory and run the following commands

```bash
mkdir build
cd build
cmake ..
make
```
After compilation there should be an executable called "cdc" in the build directory.

To run the C++ version in Microsoft Windows please cross compile or ask for a binary.

### Usage
First you need a dataset in a file. Values can be seperated by whitespace. A typical call of cdc would be for example

```bash
./cdc --file filepath/filename.txt --xindex 1 --yindex 2 --zindex 3
```

If you find any bugs, or if you experience any crashes, please report to us. If you have any questions just ask, we won't bite.
