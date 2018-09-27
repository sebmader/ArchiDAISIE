#ifndef EVENT_TYPE_H
#define EVENT_TYPE_H


enum class event_type
{
    immigration = 0,
    local_migration = 1,
    local_cladogenesis = 2,
    local_anagenesis = 3,
    local_extinction = 4,
    global_cladogenesis = 5,
    global_anagenesis = 6,
    global_extinction =7,
};

bool is_global(const event_type t) noexcept;
bool is_local(const event_type t) noexcept;


#endif // EVENT_TYPE_H
