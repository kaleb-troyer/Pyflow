
cmake_minimum_required(VERSION 3.15)
project(core VERSION 1.0 LANGUAGES CXX)

# ------------------------------
# Compiler settings
# ------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optional: treat warnings as errors
# add_compile_options(-Wall -Wextra -Werror)

# ------------------------------
# Include directories
# ------------------------------
include_directories(include)

# ------------------------------
# Source files
# ------------------------------
file(GLOB SRC_FILES
    src/*.cpp
)

# ------------------------------
# Precompiled header (optional)
# ------------------------------
# target_precompile_headers(core PRIVATE include/pch.h)  # CMake 3.16+

# ------------------------------
# Build shared library
# ------------------------------
add_library(core SHARED ${SRC_FILES})

# ------------------------------
# Define export macro for DLL
# ------------------------------
target_compile_definitions(core PRIVATE "CORE_EXPORTS")

# ------------------------------
# Platform-specific settings
# ------------------------------
if(WIN32)
    message(STATUS "Building on Windows")
else()
    message(STATUS "Building on Linux/macOS")
endif()
