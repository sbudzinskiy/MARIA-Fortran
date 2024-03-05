#ifndef _assert_mod_h_
#define _assert_mod_h_

#define ASSERT(x) call assert_linefile(x, __FILE__, __LINE__)

use maria_assert_mod, only: &
    assert_linefile

#endif
