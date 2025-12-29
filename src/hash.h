#ifndef HASH_H
#define HASH_H

#include <stdlib.h>
#include <string.h>

typedef struct hash_entry {
    void *key;
    void *value;
    struct hash_entry *next;
} s_hash_entry;

typedef size_t (*hash_func_t)(const void *key);
typedef int    (*hash_eq_t)(const void *key1, const void *key2);  /* Compares two keys */
typedef void (*hash_value_free_t)(void *value);

typedef struct hash_table {
    s_hash_entry **buckets;  /* Array of linked lists */
    size_t nbuckets;

    size_t key_size;
    size_t value_size;

    hash_func_t hash;
    hash_eq_t   equals;
    hash_value_free_t value_free;

    size_t size;   /* number of stored entries */
} s_hash_table;



static inline int hash_init(s_hash_table *ht, size_t key_size, size_t value_size, size_t nbuckets, hash_func_t hash, hash_eq_t equals, hash_value_free_t value_free);
static inline void hash_free(s_hash_table *ht);
static inline void *hash_get(s_hash_table *ht, const void *key);
static inline int hash_insert(s_hash_table *ht, const void *key, const void *value);
static inline void *hash_get_or_insert(s_hash_table *ht, const void *key);





static inline int hash_init(s_hash_table *ht, size_t key_size, size_t value_size, size_t nbuckets, hash_func_t hash, hash_eq_t equals, hash_value_free_t value_free)
{
    ht->buckets = calloc(nbuckets, sizeof(s_hash_entry *));
    if (!ht->buckets) return 0;

    ht->nbuckets = nbuckets;
    ht->key_size = key_size;
    ht->value_size = value_size;
    ht->hash = hash;
    ht->equals = equals;
    ht->value_free = value_free;
    ht->size = 0;

    return 1;
}

static inline void hash_free(s_hash_table *ht)
{
    for (size_t i = 0; i < ht->nbuckets; ++i) {
        s_hash_entry *e = ht->buckets[i];
        while (e) {
            s_hash_entry *next = e->next;
            free(e->key);

            if (ht->value_free) ht->value_free(e->value);   
            else free(e->value);

            free(e);
            e = next;
        }
    }
    free(ht->buckets);
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
        if (ht->equals(e->key, key)) return e->value;
        e = e->next;
    }
    return NULL;
}

static inline int hash_insert(s_hash_table *ht, const void *key, const void *value)
{
    size_t idx = bucket_index(ht, key);

    /* reject duplicates */
    for (s_hash_entry *e = ht->buckets[idx]; e; e = e->next)
        if (ht->equals(e->key, key)) return 0;

    s_hash_entry *e = malloc(sizeof(*e));
    if (!e) return 0;

    e->key = malloc(ht->key_size);
    e->value = malloc(ht->value_size);
    if (!e->key || !e->value) 
        { free(e->key); free(e->value); free(e); return 0; }

    memcpy(e->key, key, ht->key_size);
    memcpy(e->value, value, ht->value_size);

    e->next = ht->buckets[idx];
    ht->buckets[idx] = e;

    ht->size++;
    return 1;
}

static inline void *hash_get_or_insert(s_hash_table *ht, const void *key)
{
    size_t idx = bucket_index(ht, key);

    for (s_hash_entry *e = ht->buckets[idx]; e; e = e->next)
        if (ht->equals(e->key, key))
            return e->value;

    /* not found → create new entry */
    s_hash_entry *e = malloc(sizeof(*e));
    if (!e) return NULL;

    e->key = malloc(ht->key_size);
    e->value = calloc(1, ht->value_size);   /* ← important */
    if (!e->key || !e->value) {
        free(e->key);
        free(e->value);
        free(e);
        return NULL;
    }

    memcpy(e->key, key, ht->key_size);

    e->next = ht->buckets[idx];
    ht->buckets[idx] = e;
    ht->size++;

    return e->value;
}



#endif
