#include <annotate/profilemap.h>

/*
 * map_search
 *
 * @see include/annotate/profilemap.h
 */
int map_search(map_struct* map, char* chrom, int strand)
{
  int imax = map->size;
  int imin = 0;
  char* id = (char*) malloc(MAX_FEATURE * sizeof(char));

  // Build identifier
  strncpy(id, chrom, MAX_FEATURE);
  if (strand == FWD_STRAND)
    strcat(id, "+\0");
  else
    strcat(id, "-\0");

  imax--;
  while (imax >= imin) {
    int imid = imin + ((imax - imin) / 2);

    int result = strcmp(map->elements[imid].identifier, id);
    if(result == 0) {
      free(id);
      return imid;
    }
    else if (result < 0)
      imin = imid + 1;
    else         
      imax = imid - 1;
  }

  free(id);
  return -1;
}

/*
 * insert_sort
 *
 * @see include/annotate/profilemap.h
 */
void map_sort(map_struct* map)
{
  int i;
  int j;

  for(i = 0; i < map->size; i++) {
    char* id_x = (char*) malloc(sizeof(char) * MAX_FEATURE);
    itnode_struct* p_x;
    strncpy(id_x, map->elements[i].identifier, MAX_FEATURE);
    p_x = map->elements[i].root;

    j = i;
    while(j >= 0 && strcmp(map->elements[j - 1].identifier, id_x) > 0) {
      strncpy(map->elements[j].identifier, map->elements[j - 1].identifier, MAX_FEATURE);
      map->elements[j].root = map->elements[j - 1].root;
      j--;
    }

    strncpy(map->elements[j].identifier, id_x, MAX_FEATURE);
    map->elements[j].root = p_x;
    free(id_x);
  }
}

/*
 * map_init
 *
 * @see include/annotation/profilemap.h
 */
void map_init(map_struct* map)
{
  map->elements = NULL;
  map->size = 0;
}

/*
 * map_profile
 *
 * @see include/annotate/profilemap.h
 */
void map_add_profile(map_struct* map, profile_struct_annotation* profile)
{
  char* id = (char*) malloc(MAX_FEATURE * sizeof(char));

  // Build identifier
  strncpy(id, profile->chromosome, MAX_FEATURE);
  if (profile->strand == FWD_STRAND)
    strcat(id, "+\0");
  else
    strcat(id, "-\0");

  // Map is empty
  if (map->size == 0) {
    map->elements = (map_element_struct*) malloc(sizeof(map_element_struct) * (map->size + 1));
    strncpy(map->elements[map->size].identifier, id, MAX_FEATURE);
    map->elements[map->size].root = insert_itnode(NULL, profile->start, profile->end, profile);
    map->size++;
  }

  // Map is not empty
  else {
    int mapidx = map_search(map, profile->chromosome, profile->strand);

    // Element does not exist
    if (mapidx < 0) {
      map->elements = (map_element_struct*) realloc(map->elements, sizeof(map_element_struct) * (map->size + 1));
      strncpy(map->elements[map->size].identifier, id, MAX_FEATURE);
      map->elements[map->size].root = insert_itnode(NULL, profile->start, profile->end, profile);
      map->size++;
      map_sort(map);
    }

    // Element does exist
    else
      map->elements[mapidx].root = insert_itnode(map->elements[mapidx].root, profile->start, profile->end, profile);
  }

  free(id);
}

/*
 * map_annotate
 *
 * @see include/annotate/profilemap.h
 */
void map_annotate(map_struct* map, char* chrom, int start, int end, int strand, char* feature)
{
  int position;

  if ((position = map_search(map, chrom, strand)) < 0)
    return;

  search_itnode(map->elements[position].root, start, end, feature);
}

/*
 * map_destroy
 *
 * @see include/annotate/profilemap.h
 */
void map_destroy(map_struct* map) {
  int i;

  for (i = 0; i < map->size; i++)
    destroy_itnode(map->elements[i].root);

  free(map->elements);
}
