// Generic*Stream code from https://code.google.com/p/rapidjson/issues/detail?id=20
// Modified to match the example code given in the rapidjson documentation:
// https://github.com/miloyip/rapidjson/blob/master/doc/stream.md

#ifndef RAPIDJSON_GENERICSTREAM_H_
#define RAPIDJSON_GENERICSTREAM_H_

#include "rapidjson.h"
#include <iostream>

#ifdef _MSC_VER
  // Is disabling these warnings still required?
  #pragma warning(push)
  #pragma warning(disable: 4127) // conditional expression is constant
  #pragma warning(disable: 4512) // assignment operator could not be generated
  #pragma warning(disable: 4100) // unreferenced formal parameter
#endif

RAPIDJSON_NAMESPACE_BEGIN

//! Wrapper of std::istream for input.
class GenericReadStream {
public:
  typedef char Ch;    //!< Character type (byte).

  //! Constructor.
  /*!
    \param is Input stream.
    */
  GenericReadStream(std::istream& is) : is_(is) {
  }

  Ch Peek() const {
    int c = is_.peek();
    return c == std::char_traits<char>::eof() ? '\0' : (Ch)c;
  }

  Ch Take() {
    int c = is_.get();
    return c == std::char_traits<char>::eof() ? '\0' : (Ch)c;
  }

  size_t Tell() const {
    return (size_t)is_.tellg();
  }

  // Not implemented
  void Put(Ch)       { RAPIDJSON_ASSERT(false); }
  void Flush()       { RAPIDJSON_ASSERT(false); }
  Ch* PutBegin()     { RAPIDJSON_ASSERT(false); return 0; }
  size_t PutEnd(Ch*) { RAPIDJSON_ASSERT(false); return 0; }

private:
  GenericReadStream(const GenericReadStream&);
  GenericReadStream& operator=(const GenericReadStream&);
  
  std::istream& is_;
};


//! Wrapper of std::ostream for output.
class GenericWriteStream {
public:
  typedef char Ch;    //!< Character type. Only support char.

  //! Constructor
  /*!
    \param os Output stream.
    */
  GenericWriteStream(std::ostream& os) : os_(os) {
  }

  void Put(char c) {
    os_.put(c);
  }

  void Flush() {
    os_.flush();
  }

  size_t Tell() const {
    return (size_t)os_.tellp();
  }

  // Not implemented
  char Peek() const    { RAPIDJSON_ASSERT(false); return '\0'; }
  char Take()          { RAPIDJSON_ASSERT(false); return '\0'; }
  char* PutBegin()     { RAPIDJSON_ASSERT(false); return 0; }
  size_t PutEnd(char*) { RAPIDJSON_ASSERT(false); return 0; }

private:
  GenericWriteStream(const GenericWriteStream&);
  GenericWriteStream& operator=(const GenericWriteStream&);
  
  std::ostream& os_;
};

RAPIDJSON_NAMESPACE_END

// On MSVC, restore warnings state
#ifdef _MSC_VER
    #pragma warning(pop)
#endif
#endif // RAPIDJSON_GENERICSTREAM_H_
