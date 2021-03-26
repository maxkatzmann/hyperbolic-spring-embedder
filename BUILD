# Legacy hyperbolic spring embedder.
load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_library")

cc_library(
    name = "random",
    srcs = ["random.cc"],
    hdrs = ["random.h"],
    deps = [
        "@gslinc//:headers",
        "@gsllib//:lib",
    ],
)

cc_library(
    name = "powerlawCommon",
    srcs = ["powerlawCommon.cc"],
    hdrs = ["powerlawCommon.h"],
    deps = [
        ":random",
        "@gslinc//:headers",
        "@gsllib//:lib",
    ],
)

cc_library(
    name = "coordinate",
    srcs = ["coordinate.cc"],
    hdrs = ["coordinate.h"],
)

cc_library(
    name = "gnuplot",
    hdrs = ["gnuplot_i.h"],
)

cc_library(
    name = "force",
    srcs = ["force.cc"],
    hdrs = ["force.h"],
    deps = [
        ":coordinate",
    ],
)

cc_library(
    name = "constantforce",
    srcs = ["constantforce.cc"],
    hdrs = ["constantforce.h"],
    deps = [
        ":coordinate",
        ":force",
    ],
)

cc_library(
    name = "surfaceforce",
    srcs = ["surfaceforce.cc"],
    hdrs = ["surfaceforce.h"],
    deps = [
        ":coordinate",
        ":force",
    ],
)

cc_library(
    name = "hyperbolicspace",
    hdrs = ["hyperbolicspace.h"],
    deps = [
    ],
)

cc_library(
    name = "graph",
    srcs = ["graph.cc"],
    hdrs = ["graph.h"],
    deps = [
        ":powerlawCommon",
        ":random",
        "@com_google_glog//:glog",
    ],
)

cc_library(
    name = "NLEHelper",
    srcs = ["NLEHelper.cc"],
    hdrs = ["NLEHelper.h"],
    deps = [
        ":graph",
        ":powerlawCommon",
        ":random",
    ],
)

cc_library(
    name = "embedding",
    srcs = ["embedding.cc"],
    hdrs = ["embedding.h"],
    deps = [
        ":coordinate",
        ":force",
        ":graph",
        ":hyperbolicspace",
        "@com_google_glog//:glog",
        "@eigen",
    ],
)

cc_binary(
    name = "main",
    srcs = ["main.cc"],
    visibility = ["//visibility:public"],
    deps = [
        ":NLEHelper",
        ":constantforce",
        ":coordinate",
        ":embedding",
        ":gnuplot",
        ":hyperbolicspace",
        ":surfaceforce",
        "@com_github_gflags_gflags//:gflags",
        "@com_google_glog//:glog",
    ],
)
