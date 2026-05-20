#define CATCH_CONFIG_MAIN
#include "../catch.hpp"
#include <iostream>

void Two() 
{
    std::cout << "\nOk\n"; 
}

TEST_CASE("Test Two function", "[two]") {
    Two();
    SUCCEED();
}