#ifndef BOOST_CORE_DEMANGLE_HPP
#define BOOST_CORE_DEMANGLE_HPP

#include <string>

namespace boost {
namespace core {
    inline std::string demangle(const char* name) {
        return std::string(name);
    }
}
}

#endif
