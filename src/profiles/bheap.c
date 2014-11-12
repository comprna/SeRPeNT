#include <profiles/bheap.h>

/* initbh
 *
 * @see include/profiles/bheap.h
 */
void initbh(heap_struct* bh, int size)
{
  bh->heap_size = 0;
  bh->heap = (int*) malloc((size + 1) * sizeof(int));
  bh->alignments = (alignment_struct**) malloc((size + 1) * sizeof(alignment_struct*));
  bh->heap[0] = -INT_MAX;
}

/*
 * destroybh
 *
 * @see include/profiles/bheap.h
 */
void destroybh(heap_struct* bh)
{
  free(bh->heap);
  free(bh->alignments);
}

/*
 * insertbh
 *
 * @see include/profiles/bheap.h
 */
void insertbh(heap_struct* bh, alignment_struct* element)
{
  bh->heap_size++;

  // Insert in the last place
  bh->heap[bh->heap_size] = element->start;
  bh->alignments[bh->heap_size] = element;

  // Adjust its position
  int now = bh->heap_size;
  while(bh->heap[now/2] > element->start) {
    bh->heap[now] = bh->heap[now/2];
    bh->alignments[now] = bh->alignments[now/2];
    now /= 2;
  }
  bh->heap[now] = element->start;

  // Add pointer
  bh->alignments[now] = element;
}

/*
 * deletebh
 *
 * @see include/profiles/bheap.h
 */
alignment_struct* deletebh(heap_struct* bh)
{
  int min_element, last_element, child, now;
  alignment_struct* first_alignment;
  alignment_struct* last_alignment;

  min_element = bh->heap[1];
  first_alignment = bh->alignments[1];
  last_element = bh->heap[bh->heap_size--];
  last_alignment = bh->alignments[bh->heap_size + 1];

  for(now = 1; now * 2 <= bh->heap_size; now = child) {
    child = now * 2;

    if((child != bh->heap_size) && (bh->heap[child+1] < bh->heap[child]))
      child++;

    if(last_element > bh->heap[child]) {
      bh->heap[now] = bh->heap[child];
      bh->alignments[now] = bh->alignments[child];
    }
    else
      break;
  }

  bh->heap[now] = last_element;
  bh->alignments[now] = last_alignment;

  return first_alignment;
}
