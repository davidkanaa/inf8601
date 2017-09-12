#Lab 1
Sources for the first lab.

## Dependencies

#### Ubuntu
`apt-get install build-essential libtbb-dev pkg-config`

#### Fedora
`yum gcc gcc-c++ automake glibc-devel tbb-devel`

## Notes for compilation
By default the opimisation `-O2` is used. It is necessary to get performance results. If `gdb` is used for debugging then the source code won't match the instructions.

To compile with debugging options:
`./configure --enable-debug`


