#ifndef FERXINII_HASH_H
#define FERXINII_HASH_H

/* TODO: Right now, nbuckets is fixed. In the future, might allow for the table to grow (and rehash) if it becomes too full. */

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stddef.h>
#include <stdalign.h>

/* To be defined by the user: */
typedef size_t (*f_hash_func)(const void *key);  /* Hash function */
typedef bool   (*f_hash_key_cmp)(const void *key1, const void *key2);  /* Compares two keys. TRUE if equal, FALSE if not */
typedef void (*f_hash_value_free)(void *value);  /* OPTIONAL function to free value */


/* The hash table is an array of buckets. Each bucket is a linked list entries. Each entry is a key-value pair. */
typedef struct hash_entry { 
    struct hash_entry *next;
    /* Memory buffer follows immediately in memory */
    /* [ key bytes ][ value bytes ] */
} s_hash_entry;

typedef struct hash_table {
    s_hash_entry **buckets;  /* Array of linked lists */
    size_t nbuckets;

    size_t key_size;
    size_t value_size;

    f_hash_func hash;
    f_hash_key_cmp equals;
    f_hash_value_free value_free;

    size_t size;   /* number of stored entries */

    /* Arena support (optional) */
	void *arena;            /* base pointer of arena block, NULL if not used */
	size_t arena_capacity;  /* number of slots in arena */
	size_t arena_used;      /* how many used so far */
	size_t arena_entry_stride;    /* stride (bytes) for each entry in arena */
} s_hash_table;


static inline int hash_init(s_hash_table *ht, size_t key_size, size_t value_size, size_t nbuckets, size_t expected_entries, f_hash_func hash, f_hash_key_cmp equals, f_hash_value_free value_free);  /* If expected_entries > 0, arena is enabled */ /* 0 ERROR, 1 0K */
static inline void hash_free(s_hash_table *ht);
static inline void *hash_get(s_hash_table *ht, const void *key);  /* value* (void*) FOUND, null NOT FOUND */
static inline int hash_insert(s_hash_table *ht, const void *key, const void *value);  /* -1 DUPLICATE, 0 ERROR, 1 OK */
static inline void *hash_get_or_create(s_hash_table *ht, const void *key);  /* value* (void*) OK, null ERROR. If created, initializes entry to 0. */



/* IMPLEMENTATION */

#define HASH_ARENA_ALIGN alignof(max_align_t)

static inline size_t compute_arena_entry_stride(size_t key_size, size_t value_size)
{
	size_t header = sizeof(s_hash_entry) + key_size + value_size;
	size_t stride = (header + HASH_ARENA_ALIGN - 1) & ~(HASH_ARENA_ALIGN - 1);
	return stride;
}

static inline s_hash_entry *entry_alloc(s_hash_table *ht)
{   /* Allocate one memory blob for entry + key + value, initialized to 0 */
    if (ht->arena) {
		if (ht->arena_used < ht->arena_capacity) {
			char *base = (char*)ht->arena;
			s_hash_entry *e = (s_hash_entry*)(base + ht->arena_used * ht->arena_entry_stride);
			ht->arena_used++;
			return e;
		}
		/* Arena exhausted: fall back to heap allocation */
	}

	/* Heap allocation: exact-sized blob */
	size_t total = sizeof(s_hash_entry) + ht->key_size + ht->value_size;
	s_hash_entry *e = malloc(total);
	if (!e) return NULL;
	memset(e, 0, total);
	return e;
}

static inline bool entry_is_in_arena(const s_hash_table *ht, const s_hash_entry *e)
{
	if (!ht->arena) return false;
	char *base = (char*)ht->arena;
	char *ptr = (char*)e;
	size_t block = ht->arena_capacity * ht->arena_entry_stride;
	return (ptr >= base) && (ptr < base + block);
}

static inline void *entry_key(const s_hash_entry *e)
{
	return (void*)( (char*)e + sizeof(s_hash_entry) );
}

static inline void *entry_value(const s_hash_table *ht, const s_hash_entry *e)
{
	return (void*)( (char*)e + sizeof(s_hash_entry) + ht->key_size );
}


static inline int hash_init(s_hash_table *ht, size_t key_size, size_t value_size, size_t nbuckets, size_t expected_entries, f_hash_func hash, f_hash_key_cmp equals, f_hash_value_free value_free)
{
    if (nbuckets < 1) { fprintf(stderr, "hash_init: nbuckets needs to be >= 1.\n"); return 0; }
    ht->buckets = calloc(nbuckets, sizeof(s_hash_entry*));
    if (!ht->buckets) return 0;

    ht->nbuckets = nbuckets;
    ht->key_size = key_size;
    ht->value_size = value_size;
    ht->hash = hash;
    ht->equals = equals;
    ht->value_free = value_free;
    ht->size = 0;

    /* Arena */
	ht->arena = NULL;
	ht->arena_capacity = 0;
	ht->arena_used = 0;
	ht->arena_entry_stride = 0;
    if (expected_entries > 0) {  /* compute stride and allocate arena block */
        ht->arena_entry_stride = compute_arena_entry_stride(key_size, value_size);

        size_t total = ht->arena_entry_stride * expected_entries;
        void *arena = malloc(total);
        if (!arena) { fprintf(stderr, "hash_init: Could not initialize arena.\n"); free(ht->buckets); return 0; }
        ht->arena = arena;
        ht->arena_capacity = expected_entries;
        ht->arena_used = 0;

        memset(arena, 0, total);  /* Zero the whole arena memory for safety */
    }

    return 1;
}

static inline void hash_free(s_hash_table *ht)
{
    for (size_t i = 0; i < ht->nbuckets; i++) {
        s_hash_entry *e = ht->buckets[i];
        while (e) {
            s_hash_entry *next = e->next;
            
            /* Free value if user provided special free function */
            if (ht->value_free) { 
                void *vptr = entry_value(ht, e);
                ht->value_free(vptr);   
            }

            /* Free entry memory if it does not belong to the arena */
			if (!entry_is_in_arena(ht, e)) {
				free(e);
			}

            e = next;
        }
    }

    if (ht->arena) free(ht->arena);
    free(ht->buckets);

    memset(ht, 0, sizeof(s_hash_table));
}

static inline size_t bucket_index(s_hash_table *ht, const void *key)
{
    return ht->hash(key) % ht->nbuckets;
}

static inline void *hash_get(s_hash_table *ht, const void *key)
{
    size_t idx = bucket_index(ht, key);
    s_hash_entry *e = ht->buckets[idx];

    while (e) {
        if (ht->equals(entry_key(e), key)) { return entry_value(ht, e); }
        e = e->next;
    }
    return NULL;
}

static inline int hash_insert(s_hash_table *ht, const void *key, const void *value)
{
    size_t idx = bucket_index(ht, key);

    /* Check if entry already exists in linked list */
    for (s_hash_entry *e = ht->buckets[idx]; e; e = e->next)
        if (ht->equals(entry_key(e), key)) return -1;

    /* Malloc entry */
    s_hash_entry *e = entry_alloc(ht);
    if (!e) return 0;

    void *kptr = entry_key(e);
    void *vptr = entry_value(ht, e);
    memcpy(kptr, key, ht->key_size);
    memcpy(vptr, value, ht->value_size);

    e->next = ht->buckets[idx];
    ht->buckets[idx] = e;
    ht->size++;

    return 1;
}

static inline void *hash_get_or_create(s_hash_table *ht, const void *key)
{
    size_t idx = bucket_index(ht, key);

    for (s_hash_entry *e = ht->buckets[idx]; e; e = e->next)
        if (ht->equals(entry_key(e), key)) return entry_value(ht, e);

    /* Not found: create new entry */
    s_hash_entry *e = entry_alloc(ht);
    if (!e) return NULL;

    void *kptr = entry_key(e);
    memcpy(kptr, key, ht->key_size);
    void *vptr = entry_value(ht, e);  /* Set value to 0 */
    memset(vptr, 0, ht->value_size);
    
    e->next = ht->buckets[idx];
    ht->buckets[idx] = e;
    ht->size++;

    return entry_value(ht, e);
}




#endif
