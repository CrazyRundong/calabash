# Calabash
This is Rundong Li's C++ implementation of `CS240`: *Angorithm Design and Analysis* 
[course project](https://github.com/meijun/calabash).

## How to use
Requirements:
- `Eigen3`
- `Intel MKL`
- `Visual Studio 2017` or `cmake`

### Windows
Please install `Visual Studio 2017` with `C++` support, compile the `Release` or `ReleaseWithInfo`
binary, then execute
```bash
calabash.exe <num_of_steps> <graph_file_path> <output_file_path>
```

### Linux
Please install `cmake` and `libeigen3-dev` via
```bash
apt install cmake libeigen3-dev
```
then build the project under ignored dictionary
```bash
mkdir build && cd build
cmake -DDEBUG_DIAMOND .. # use DEBUG_DIAMOND symbol to activate file stream
./calabash <num_of_steps> <graph_file_path> <output_file_path>
```

> This repostory will remaining private before this course ends.
