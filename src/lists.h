#ifndef LISTS_H
#define LISTS_H

#include <string.h>
#include <stdlib.h>

#define LIST_INIT_DEFAULT_NMAX 10

typedef struct list {
	void  *items;        /* pointer to raw bytes */
	size_t elem_size;    /* sizeof(element type) */
	unsigned int N;
	unsigned int Nmax;
} s_list;


static inline s_list list_initialize(size_t elem_size, size_t Nmax);
static inline int list_ensure_capacity(s_list *list, size_t need);
static inline int list_push(s_list *list, const void *elem);
static inline int list_change_entry(s_list *list, unsigned id, const void *elem);
static inline void *list_get_ptr(s_list *list, size_t idx);
static inline int list_get_value(const s_list *list, size_t id, void *out);
static inline int list_pop(s_list *list, void *out_elem);
static inline void free_list(s_list *list);


static inline s_list list_initialize(size_t elem_size, size_t Nmax)
{
	if (elem_size == 0) elem_size = 1;
	if (Nmax == 0) Nmax = LIST_INIT_DEFAULT_NMAX;
	s_list out = { 
        .items     = malloc(Nmax * elem_size),
		.elem_size = elem_size,
		.N         = 0,
		.Nmax      = Nmax,
	};
	return out;
}

static inline int list_ensure_capacity(s_list *list, size_t need) 
{   
	if (!list) return 0;
	if (need < list->Nmax) return 1;

	size_t Nmax_new = list->Nmax ? list->Nmax : 1;
	while (need >= Nmax_new) Nmax_new *= 2;

	void *tmp = realloc(list->items, Nmax_new * list->elem_size);
	if (!tmp) return 0;

	list->items = tmp;
	list->Nmax = Nmax_new;
	return 1;
}

static inline int list_push(s_list *list, const void *elem) 
{
	if (!list) return 0;
	if (!list_ensure_capacity(list, list->N + 1)) return 0;
	void *dst = (char *)list->items + list->N * list->elem_size;
	memcpy(dst, elem, list->elem_size);
	list->N += 1;
	return 1;
}

static inline int list_change_entry(s_list *list, unsigned id, const void *elem)
{
    if (!list || id >= list->N) return 0;

    void *dst = (uint8_t*)list->items + id * list->elem_size;
    memmove(dst, elem, list->elem_size);
    return 1;
}

static inline void *list_get_ptr(s_list *list, size_t id)
{
	if (!list || id >= list->N) return NULL;
	return (uint8_t*)list->items + id * list->elem_size;
}

static inline int list_get_value(const s_list *list, size_t id, void *out)
{
	if (!list || !out || id >= list->N) return 0;

	const void *src = (uint8_t*)list->items + id * list->elem_size;
	memmove(out, src, list->elem_size);
	return 1;
}

static inline int list_pop(s_list *list, void *out_elem)
{
	if (!list || list->N == 0) return 0;
	list->N -= 1;
	void *src = (uint8_t*)list->items + list->N * list->elem_size;
	if (out_elem) memcpy(out_elem, src, list->elem_size);
	return 1;
}

static inline void free_list(s_list *list) 
{
	if (!list) return;
	if (list->items) free(list->items);
	memset(list, 0, sizeof(*list));
}

#endif
