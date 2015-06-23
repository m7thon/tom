#ifndef MACROS_H
#define MACROS_H

// NOTE: The lowest and recommended checking level, especially when working from python, is
//       TOM_CHECK. Enable TOM_DEBUG for full checking (also in libraries) and TOM_CHECK to
//       turn on asserts only in the tom library.

#ifdef TOM_DEBUG
    #undef NDEBUG
    #define TOM_CHECK
#endif

#ifdef TOM_CHECK
    #define CHECK(...) __VA_ARGS__
#else
    #define CHECK(...)
#endif

#ifdef SWIG
    #define SWIGCODE(...) __VA_ARGS__
    #define NEWOBJECT(func) %newobject func;
    #define TEMPLATE(pyname,name,...) %template(pyname) name<__VA_ARGS__ >;
    #define OUTMATRIXXD(mat) %apply Eigen::MatrixXd& OUTPUT { Eigen::MatrixXd& mat };
    #define CLEAROUTMATRIXXD(mat) %clear Eigen::MatrixXd& mat;
    #define OUTMAP(mat) %apply const Eigen::MatrixBase<Eigen::MatrixXd>& OUTPUT { const Eigen::MatrixBase<Eigen::MatrixXd>& mat };
    #define CLEAROUTMAP(mat) %clear const Eigen::MatrixBase<Eigen::MatrixXd>& mat;
    #define NOSWIG(...)
#else
    #define SWIGCODE(...)
    #define NEWOBJECT(func)
    #define TEMPLATE(pyname,name,...)
    #define OUTMATRIXXD(mat)
    #define CLEAROUTMATRIXXD(mat)
    #define OUTMAP(mat)
    #define CLEAROUTMAP(mat)
    #define NOSWIG(...) __VA_ARGS__
#endif

#ifdef VERBOSE
    #define PG_INFO(s) {std::cerr << s; std::cerr.flush();}
    #define PG_MARK(i,max) {if ((i + 1) % ((max) / 50 + 1) == 0) {std::cerr << "."; std::cerr.flush();}}
    #define PG_DONE {std::cerr << "done!" << std::endl; std::cerr.flush();}
#else
    #define PG_INFO(s)
    #define PG_MARK(i,max)
    #define PG_DONE
#endif

#define ERR(s) {std::cerr << s << std::endl; std::cerr.flush();}

#endif // MACROS_H
