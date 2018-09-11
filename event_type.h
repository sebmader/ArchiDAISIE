#ifndef EVENT_TYPE_H
#define EVENT_TYPE_H


enum class event_type
{
    local_migration,
    local_cladogenesis,
    local_anagenesis,
    local_extinction,
    immigration,
    global_anagenesis,
    global_cladogenesis,
    global_extinction,
};

bool is_global(const event_type t) noexcept;
bool is_local(const event_type t) noexcept;


#endif // EVENT_TYPE_H
