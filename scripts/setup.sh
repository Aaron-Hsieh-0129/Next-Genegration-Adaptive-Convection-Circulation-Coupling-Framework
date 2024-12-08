#!/bin/bash

# Define installation directories
INSTALL_DIR="$HOME/local"
BIN_DIR="$INSTALL_DIR/bin"
LIB_DIR="$INSTALL_DIR/lib"
LIB64_DIR="$INSTALL_DIR/lib64"
INCLUDE_DIR="$INSTALL_DIR/include"

# Ensure directories exist
mkdir -p "$BIN_DIR" "$LIB_DIR" "$INCLUDE_DIR"

# Function to check if a command exists
check_command() {
    command -v "$1" >/dev/null 2>&1
}

# Function to compare versions
# Returns 0 if the installed version is greater or equal to the required version, 1 otherwise.
version_ge() {
    local installed_version="$1"
    local required_version="$2"
    [ "$(printf '%s\n' "$required_version" "$installed_version" | sort -V | head -n1)" = "$required_version" ]
}

# Prompt user for yes/no
prompt_yes_no() {
    while true; do
        read -p "$1 (y/n): " choice
        case "$choice" in
            y|Y ) return 0 ;;
            n|N ) return 1 ;;
            * ) echo "Please enter y or n." ;;
        esac
    done
}

# Install CMake
install_cmake() {
    local CMAKE_VERSION=3.27.2
    if check_command cmake; then
        local INSTALLED_VERSION
        INSTALLED_VERSION=$(cmake --version | head -n1 | awk '{print $3}')
        if version_ge "$INSTALLED_VERSION" "$CMAKE_VERSION"; then
            echo "CMake version $INSTALLED_VERSION is already installed and satisfies the requirement."
            return
        fi
    fi

    echo "Installing CMake (required version: $CMAKE_VERSION)..."
    CMAKE_URL="https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-linux-x86_64.tar.gz"
    wget "$CMAKE_URL" -O cmake.tar.gz
    tar -xzf cmake.tar.gz --strip-components=1 -C "$INSTALL_DIR"
    rm cmake.tar.gz
    echo "CMake $CMAKE_VERSION installed successfully."
}

# Install GCC
install_gcc() {
    local GCC_VERSION=13.2.0
    if check_command gcc; then
        local INSTALLED_VERSION
        INSTALLED_VERSION=$(gcc --version | head -n1 | awk '{print $NF}')
        if version_ge "$INSTALLED_VERSION" "$GCC_VERSION"; then
            echo "GCC version $INSTALLED_VERSION is already installed and satisfies the requirement."
            return
        fi
    fi

    echo "Installing GCC (required version: $GCC_VERSION)..."
    GCC_URL="https://ftp.mirrorservice.org/sites/sourceware.org/pub/gcc/releases/gcc-${GCC_VERSION}/gcc-${GCC_VERSION}.tar.gz"
    wget "$GCC_URL" -O gcc.tar.xz
    tar -xf gcc.tar.xz -C "$INSTALL_DIR"
    rm gcc.tar.xz
    cd "$INSTALL_DIR/gcc-$GCC_VERSION"
    ./contrib/download_prerequisites
    ./configure --prefix=$INSTALL_DIR --enable-languages=c,c++,fortran --disable-multilib
    make -j 32
    make install
    echo "GCC $GCC_VERSION installed successfully."
}

install_hdf5() {
    local HDF5_VERSION=1.14.4
    HDF5_VERSION_UNDERSCORE=${HDF5_VERSION//./_}
    if [ ! -d "$INCLUDE_DIR/hdf5" ]; then
        echo "Installing hdf5 (version: $HDF5_VERSION)..."
        HDF5_URL="https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5_${HDF5_VERSION}.tar.gz"
        wget "$HDF5_URL" -O hdf5.tar.gz
        tar xf hdf5.tar.gz
        cd hdf5-hdf5_${HDF5_VERSION}
        ./configure --prefix=${INSTALL_DIR} --enable-fortran --enable-cxx --enable-shared --enable-hl 
        make -j 32
        make install
        cd ..
	rm -rf hdf5.tar.gz hdf5-hdf5_${HDF5_VERSION}
        echo "hdf5 $HDF5_VERSION installed successfully."
    else
        echo "hdf5 is already installed in $INCLUDE_DIR/hdf5."
    fi
}


install_netcdfc() {
    local NETCDFC_VERSION=4.3.3.1
    NETCDFC_VERSION_UNDERSCORE=${NETCDFC_VERSION//./_}
    if [ ! -d "$INCLUDE_DIR/netcdf-c" ]; then
        echo "Installing netcdf-c (version: $NETCDFC_VERSION)..."
        NETCDFC_URL="https://github.com/Unidata/netcdf-c/archive/v${NETCDFC_VERSION}.tar.gz"
        wget "$NETCDFC_URL" -O netcdfc.tar.gz
        tar xf netcdfc.tar.gz
        cd netcdf-c-${NETCDFC_VERSION}
        ./configure --prefix=${INSTALL_DIR} --enable-netcdf-4 LDFLAGS=-L$INSTALL_DIR/lib CFLAGS=-I$INSTALL_DIR/include 
        make -j 32
        make install
        cd ..
	rm -rf netcdfc.tar.gz netcdf-c-${NETCDFC_VERSION}
        echo "netcdf-c $NETCDFC_VERSION installed successfully."
    else
        echo "netcdf-c is already installed in $INCLUDE_DIR/netcdf-c."
    fi
}

install_netcdfcxx() {
    local NETCDFCXX_VERSION=4.2.1
    NETCDFCXX_VERSION_UNDERSCORE=${NETCDFCXX_VERSION//./_}
    if [ ! -d "$INCLUDE_DIR/netcdf-cxx" ]; then
        echo "Installing netcdf-cxx (version: $NETCDFCXX_VERSION)..."
        NETCDFCXX_URL="https://github.com/Unidata/netcdf-cxx4/archive/refs/tags/v${NETCDFCXX_VERSION}.tar.gz"
        wget "$NETCDFCXX_URL" -O netcdfcxx.tar.gz
        tar xf netcdfcxx.tar.gz
        cd netcdf-cxx4-${NETCDFCXX_VERSION}
        ./configure --prefix=${INSTALL_DIR} 
        make -j 32
        make install
        cd ..
	rm -rf netcdfcxx.tar.gz netcdf-cxx4-${NETCDFCXX_VERSION}
        echo "netcdf-cxx $NETCDFCXX_VERSION installed successfully."
    else
        echo "netcdf-cxx is already installed in $INCLUDE_DIR/netcdf-cxx."
    fi
}

# Update shell configuration file
update_shell_config() {
    SHELL_CONFIG_FILE="$HOME/.bashrc"
    if [ "$SHELL" = "/bin/zsh" ] || [ "$SHELL" = "/usr/bin/zsh" ]; then
        SHELL_CONFIG_FILE="$HOME/.zshrc"
    fi

    echo "Updating $SHELL_CONFIG_FILE..."
    {
        echo "# Added by setup script"
        echo "export PATH=\"$BIN_DIR:\$PATH\""
        echo "export LD_LIBRARY_PATH=\"$LIB_DIR:\$LD_LIBRARY_PATH\""
        echo "export LD_LIBRARY_PATH=\"$LIB64_DIR:\$LD_LIBRARY_PATH\""
        echo "export LIBRARY_PATH=\"$LIB_DIR:\$LIBRARY_PATH\""
        echo "export LIBRARY_PATH=\"$LIB64_DIR:\$LIBRARY_PATH\""
        echo "export CMAKE_PREFIX_PATH=\"$INSTALL_DIR:\$CMAKE_PREFIX_PATH\""
        echo "export C_INCLUDE_PATH=\"$INCLUDE_DIR:\$C_INCLUDE_PATH\""
        echo "export CPLUS_INCLUDE_PATH=\"$INCLUDE_DIR:\$CPLUS_INCLUDE_PATH\""
    } >> "$SHELL_CONFIG_FILE"
    echo "Environment variables added to $SHELL_CONFIG_FILE. Please restart your shell or run 'source $SHELL_CONFIG_FILE' to apply changes."
}

# Main script execution
echo "Welcome to the setup script! Please select the tools to install."
if prompt_yes_no "Update your shell configuration file with environment variables?"; then
    update_shell_config
    source $SHELL_CONFIG_FILE 
fi

if prompt_yes_no "Install or update CMake?"; then
    install_cmake
fi

if prompt_yes_no "Install or update GCC? (This would take about 10~30 mins depending on the core numbers)"; then
    install_gcc
fi

if prompt_yes_no "Install hdf5?"; then
    install_hdf5
fi

if prompt_yes_no "Install netcdf-c?"; then
    install_netcdfc
fi

if prompt_yes_no "Install netcdf-cxx?"; then
    install_netcdfcxx
fi

echo "Setup completed!"

