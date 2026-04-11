#include "SpcStats_fns.hpp"

namespace saint_spc {

const double log_factorial_table[log_factorial_max + 1] = 
#include "SpcLogFactorial.hpp"
;

} // namespace saint_spc
