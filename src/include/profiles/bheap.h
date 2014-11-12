#include <core/structs.h>

/*
 * initbh
 *   Initialize heap
 *
 * @arg heap_struct* bh
 *   A pointer to the heap data structure
 */
void initbh(heap_struct* bh, int size);

/*
 * destroybh
 *   Free and destroy heap
 *
 * @arg heap_struct* bh
 *   A pointer to the heap data structure
 */
void destroybh(heap_struct* bh);

/*
 * insertbh
 *   Insert an element into the heap
 *
 * @arg heap_struct* bh
 *   A pointer to the heap data structure
 *
 * @arg int element
 *   New element inserted in the heap
 */
void insertbh(heap_struct* bh, alignment_struct* element);

/*
 * deletebh
 *   Delete and return the top priority element in the heap
 *
 * @arg heap_struct* bh
 *   A pointer to the heap data structure
 *
 * @return
 *   The top priority element in the heap
 */
alignment_struct* deletebh(heap_struct* bh);
