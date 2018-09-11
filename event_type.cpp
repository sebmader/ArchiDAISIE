#include "event_type.h"

bool is_global(const event_type t) noexcept
{
    return t == event_type::global_cladogenesis
      || t == event_type::global_anagenesis
      || t == event_type::global_extinction
    ;
}

bool is_local(const event_type t) noexcept
{
    return !is_global(t);
}
