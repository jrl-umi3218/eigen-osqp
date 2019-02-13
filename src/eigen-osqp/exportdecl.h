#pragma once

#if defined _WIN32 || defined __CYGWIN__
#define EIGEN_OSQP_DLLIMPORT __declspec(dllimport)
#define EIGEN_OSQP_DLLEXPORT __declspec(dllexport)
#define EIGEN_OSQP_DLLLOCAL
#else
#if __GNUC__ >= 4
#define EIGEN_OSQP_DLLIMPORT __attribute__((visibility("default")))
#define EIGEN_OSQP_DLLEXPORT __attribute__((visibility("default")))
#define EIGEN_OSQP_DLLLOCAL __attribute__((visibility("hidden")))
#else
#define EIGEN_OSQP_DLLIMPORT
#define EIGEN_OSQP_DLLEXPORT
#define EIGEN_OSQP_DLLLOCAL
#endif
#endif

#ifdef EIGEN_OSQP_STATIC
#define EIGEN_OSQP_DLLAPI
#define EIGEN_OSQP_LOCAL
#else
#ifdef EIGEN_OSQP_EXPORTS
#define EIGEN_OSQP_DLLAPI EIGEN_OSQP_DLLEXPORT
#else
#define EIGEN_OSQP_DLLAPI EIGEN_OSQP_DLLIMPORT
#endif
#define EIGEN_OSQP_LOCAL EIGEN_OSQP_DLLLOCAL
#endif