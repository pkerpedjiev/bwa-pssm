#include "bwtge.h"

/*
void *__real_malloc (size_t);
void *__real_calloc (size_t, size_t);
void *__real_realloc(void *, size_t);

void *
__wrap_realloc (void *b, size_t c)
{
    fprintf (stderr,"realloc called with %zu\n", c);
    return (void *)__real_realloc (b, c);
}


void *
__wrap_calloc (size_t c, size_t b)
{
    fprintf (stderr,"calloc called with %zu\n", c);
    return (void *)__real_calloc (c, b);
}


void *
__wrap_malloc (size_t c)
{
    fprintf (stderr,"malloc called with %zu\n", c);
    return (void *)__real_malloc (c);
}

*/

extern gdsl_element_t 
alloc_pssm_entry (void *pssm_entry)
{
    /*
    pssm_entry_t *entry = (pssm_entry_t *) pssm_entry;
    pssm_entry_t *value = (pssm_entry_t *) malloc (sizeof (pssm_entry_t));

    assert(value != NULL);

    memcpy(value, entry, sizeof(pssm_entry_t));
    */

    return (gdsl_element_t) pssm_entry;
}

extern void 
free_pssm_entry (gdsl_element_t e)
{
    //free (e);
}

#define bfs

extern long int
compare_pssm_entries (gdsl_element_t e1, void* e2)
{
    float val1, val2;

#ifdef bfs
        val1 =  ((pssm_entry_t *) e1)->score_offset;
        val2 =  ((pssm_entry_t *) e2)->score_offset;
#else
        val1 =  ((pssm_entry_t *) e1)->pssm_score;
        val2 =  ((pssm_entry_t *) e2)->pssm_score;
#endif

    if (val1 < val2)
        return -1;
    else {
        if (val2 < val1)
            return 1;
    }
    return 0;

#ifdef bfs
        val1 =  ((pssm_entry_t *) e1)->pssm_score;
        val2 =  ((pssm_entry_t *) e2)->pssm_score;
#else
        val1 =  ((pssm_entry_t *) e1)->score_offset;
        val2 =  ((pssm_entry_t *) e2)->score_offset;
#endif

    //fprintf(stderr, "val1: %f val2: %f\n", val1, val2);
    if (val1 < val2)
        return -1;
    else {
        if (val2 < val1)
            return 1;
    }

    return 0;
}



