\page metatensor Metatensor interface 

<!--
description: Module that implements interface with metatensor 
authors: Guillaume Fraux
reference: 
-->

The metatensor module allows one to call metatensor, which is a specialized dat storage format for atomistic machine learning, from PLUMED.
You can learn more about metatensor by reading [the metatensor manual](https://lab-cosmo.github.io/metatensor/latest/).

To use metatensor from plumed you use the \ref METATENSOR action.

\section metatensor-build Building PLUMED with metatensor

You'll need to fist install libtorch, either by installing PyTorch itself
with Python, or by downloading the prebuilt C++ library from
[here](https://pytorch.org/get-started/locally/).  Once you have installed libtorch you 
need to set the following options in bash.

\verbatim
# point this to the path where you extracted the C++ libtorch
TORCH_PREFIX=../../..
# if you used Python to install torch, you can do this:
TORCH_CMAKE_PREFIX=$(python -c "import torch; print(torch.utils.cmake_prefix_path)")
TORCH_PREFIX=$(cd "$TORCH_CMAKE_PREFIX/../.." && pwd)

TORCH_INCLUDES="-I$TORCH_PREFIX/include -I$TORCH_PREFIX/include/torch/csrc/api/include"
\endverbatim

You can then either build and install metatensor-torch from source. You'll need a rust
compiler on your system.  The easiest way to install a rust compiler is to use [the following site](https://rustup.rs/).
Once the rust compiler is installed you should then follow the instructions that follow: 

\verbatim
# patch a bug from torch's MKL detection
cd <PLUMED/DIR>
./src/metatensor/patch-torch.sh "$TORCH_PREFIX"

cd <SOME/PLACE/WHERE/TO/PUT/METATENSOR/SOURCES>

# define a location where metatensor should be installed
METATENSOR_PREFIX=<...>

METATENSOR_TORCH_PREFIX="$METATENSOR_PREFIX"

git clone https://github.com/lab-cosmo/metatensor --branch=metatensor-torch-v0.4.0
cd metatensor

mkdir build && cd build
cmake -DBUILD_SHARED_LIBS=ON \
      -DCMAKE_INSTALL_PREFIX="$METATENSOR_PREFIX" \
      -DCMAKE_PREFIX_PATH="$TORCH_PREFIX" \
      -DBUILD_METATENSOR_TORCH=ON \
      -DMETATENSOR_INSTALL_BOTH_STATIC_SHARED=OFF \
      ..

cmake --build . --target install --parallel
\endverbatim

Instead of doing the above, you can alternatively, use metatensor-torch from Python (`pip install metatensor[torch]`).
Once the installation has completed you then enter the following commands:

\verbatim
METATENSOR_CMAKE_PREFIX=$(python -c "import metatensor; print(metatensor.utils.cmake_prefix_path)")
METATENSOR_PREFIX=$(cd "$METATENSOR_CMAKE_PREFIX/../.." && pwd)

METATENSOR_TORCH_CMAKE_PREFIX=$(python -c "import metatensor.torch; print(metatensor.torch.utils.cmake_prefix_path)")
METATENSOR_TORCH_PREFIX=$(cd "$METATENSOR_TORCH_CMAKE_PREFIX/../.." && pwd)
\endverbatim

Once these two steps are completed you can then install plumed itself using the following commands:

\verbatim
cd <PLUMED/DIR>

# set the rpath to make sure plumed executable will be able to find the right libraries
RPATH="-Wl,-rpath,$TORCH_PREFIX/lib -Wl,-rpath,$METATENSOR_PREFIX/lib -Wl,-rpath,$METATENSOR_TORCH_PREFIX/lib"

# configure PLUMED with metatensor
./configure --enable-libtorch --enable-metatensor --enable-modules=+metatensor \
    LDFLAGS="-L$TORCH_PREFIX/lib -L$METATENSOR_PREFIX/lib -L$METATENSOR_TORCH_PREFIX/lib $RPATH" \
    CPPFLAGS="$TORCH_INCLUDES -I$METATENSOR_PREFIX/include -I$METATENSOR_TORCH_PREFIX/include"

# If you are on Linux and use a pip-installed version of libtorch, or the
# pre-cxx11-ABI build of libtorch, you'll need to add "-D_GLIBCXX_USE_CXX11_ABI=0"
# to the compilation flags:
./configure --enable-libtorch --enable-metatensor --enable-modules=+metatensor \
    LDFLAGS="-L$TORCH_PREFIX/lib -L$METATENSOR_PREFIX/lib -L$METATENSOR_TORCH_PREFIX/lib $RPATH" \
    CPPFLAGS="$TORCH_INCLUDES -I$METATENSOR_PREFIX/include -I$METATENSOR_TORCH_PREFIX/include" \
    CXXFLAGS="-D_GLIBCXX_USE_CXX11_ABI=0"

make -j && make install
\endverbatim


<!-- TODO: explain vesin update process -->
