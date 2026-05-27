#include "catch.hpp"

extern "C"
{
    int uadd(unsigned, unsigned);
}

TEST_CASE( "TEST BX: " )
{    
    SECTION( "test BX" )
    {
        REQUIRE( uadd(5, 4) == 9 );
        REQUIRE( uadd(6, 4) == 9 );
        
    }
}