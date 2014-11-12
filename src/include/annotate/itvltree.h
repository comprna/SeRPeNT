#include <core/structs.h>


/*
 * insert_itnode
 *   A utility function to insert a new Interval Search Tree Node
 *
 * @arg itnode_struct *node
 *   A pointer to the root node of the interval tree
 * @arg int low
 *   Lower bound of the new interval
 * @arg int high
 *   Upper bound of the new interval
 *
 * @return
 *   The root node of the interval tree
 */
itnode_struct* insert_itnode(itnode_struct* node, int low, int high, profile_struct* profile);


/*
 * search_itnode
 *   Aearches a given interval i in a given interval tree
 *
 * @arg itnode_struct *node
 *   A pointer to the root node of the interval tree
 * @arg int low
 *   Lower bound of the new interval to be found
 * @arg int high
 *   Upper bound of the new interval to be found
 *
 * @return
 *   The index of the profile that overlaps with the interval.
 *   -1 if no profiles overlaps with the interval.
 */
void search_itnode(itnode_struct* root, int low, int high, char* annotation);


/*
 * destroy_itnode
 *   Frees all the nodes that are rooted from the given node
 *
 * @arg itnode_struct *node
 *   A pointer to the root node of the interval tree
 */
void destroy_itnode(itnode_struct* root);
