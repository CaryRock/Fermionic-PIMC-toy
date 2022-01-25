#include <cstdlib>
#include <cstdio>
#include <string>

#include "jlcxx/jlcxx.hpp"

std::string greet()
{
    return "Hello, world!\n";
}

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
    mod.method("greet", &greet);
}
