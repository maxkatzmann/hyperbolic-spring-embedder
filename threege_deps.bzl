"""Load dependencies needed to compile the spring embedder."""

load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

def threege_deps():
    """Loads common dependencies needed to compile the embedder."""

    """gflags"""
    if not native.existing_rule("com_github_gflags_gflags"):
        http_archive(
            name = "com_github_gflags_gflags",
            sha256 = "34af2f15cf7367513b352bdcd2493ab14ce43692d2dcd9dfc499492966c64dcf",
            strip_prefix = "gflags-2.2.2",
            urls = [
                "https://github.com/gflags/gflags/archive/v2.2.2.tar.gz",
            ],
        )

    """glog"""
    if not native.existing_rule("com_google_glog"):
        http_archive(
            name = "com_google_glog",
            sha256 = "9e1b54eb2782f53cd8af107ecf08d2ab64b8d0dc2b7f5594472f3bd63ca85cdc",
            strip_prefix = "glog-0.4.0",
            urls = ["https://github.com/google/glog/archive/v0.4.0.zip"],
        )

    if not native.existing_rule("eigen"):
        http_archive(
            name = "eigen",
            build_file_content = """
# Description:
#   Eigen is a C++ template library for linear algebra: vectors,
#   matrices, and related algorithms.

licenses([
    # Note: Eigen is an MPL2 library that includes GPL v3 and LGPL v2.1+ code.
    #       We've taken special care to not reference any restricted code.
    "reciprocal",  # MPL2
    "notice",  # Portions BSD
])

exports_files(["COPYING.MPL2"])

EIGEN_FILES = [
    "Eigen/**",
    "unsupported/Eigen/CXX11/**",
    "unsupported/Eigen/FFT",
    "unsupported/Eigen/KroneckerProduct",
    "unsupported/Eigen/src/FFT/**",
    "unsupported/Eigen/src/KroneckerProduct/**",
    "unsupported/Eigen/MatrixFunctions",
    "unsupported/Eigen/SpecialFunctions",
    "unsupported/Eigen/src/MatrixFunctions/**",
    "unsupported/Eigen/src/SpecialFunctions/**",
]

# Files known to be under MPL2 license.
EIGEN_MPL2_HEADER_FILES = glob(
    EIGEN_FILES,
    exclude = [
        # Guarantees that any non-MPL2 file added to the list above will fail to
        # compile.
        "Eigen/src/Core/util/NonMPL2.h",
        "Eigen/**/CMakeLists.txt",
    ],
)

cc_library(
    name = "eigen",
    hdrs = EIGEN_MPL2_HEADER_FILES,
    defines = [
        # This define (mostly) guarantees we don't link any problematic
        # code. We use it, but we do not rely on it, as evidenced above.
        "EIGEN_MPL2_ONLY",
        "EIGEN_MAX_ALIGN_BYTES=64",
        "EIGEN_HAS_TYPE_TRAITS=0",
    ],
    includes = ["."],
    visibility = ["//visibility:public"],
)

filegroup(
    name = "eigen_header_files",
    srcs = EIGEN_MPL2_HEADER_FILES,
    visibility = ["//visibility:public"],
)
""",
            sha256 = "e09b89aae054e9778ee3f606192ee76d645eec82c402c01c648b1fe46b6b9857",
            strip_prefix = "eigen-3.3.7",
            url = "https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.zip",
        )

    """GSL"""
    if not native.existing_rule("gsllib"):
        native.new_local_repository(
            name = "gsllib",
            build_file_content = """
cc_library(
name = "lib",
srcs = ["libgsl.a"],
visibility = ["//visibility:public"],
)
""",
            path = "/usr/local/Cellar/gsl/2.6/lib",
        )

    if not native.existing_rule("gslinc"):
        native.new_local_repository(
            name = "gslinc",
            build_file_content = """
package(default_visibility = ["//visibility:public"])
cc_library(
name = "headers",
hdrs = glob(["**/*.h"]),
includes = ["."],
)
""",
            path = "/usr/local/Cellar/gsl/2.6/include",
        )
