#include <annotate/itvltree.h>

/*
 * max
 *   A utility function to get maximum of two integers
 */
int max(int a, int b)
{
  return (a > b)? a : b;
}

/*
 * height
 *   A utility function to get height of the tree
 */
int height(itnode_struct* root)
{
  if (root == NULL)
    return 0;

  return root->height;
}

/*
 * do_overlap
 *   A utility function to check if given two intervals overlap
 */
int do_overlap(int low_a, int high_a, int low_b, int high_b)
{
  if (low_a <= high_b && low_b <= high_a) return 1;
  return 0;
}

/*
 * new_itnode
 *   Helper function that allocates a new node with the given key and
 *   NULL left and right pointers
 */
itnode_struct* new_itnode(int low, int high, profile_struct_annotation* profile)
{
  itnode_struct* node = (itnode_struct*) malloc(sizeof(struct itnode_struct));

  node->low = low;
  node->high = high;
  node->profile = profile;
  node->left   = NULL;
  node->right  = NULL;
  node->height = 1;  // new node is initially added at leaf

  return(node);
}

/*
 * right_rotate
 *   A utility function to right rotate subtree rooted with y
 */
itnode_struct* right_rotate(itnode_struct* y)
{
    itnode_struct* x = y->left;
    itnode_struct* T2 = x->right;

    // Perform rotation
    x->right = y;
    y->left = T2;
 
    // Update heights
    y->height = max(height(y->left), height(y->right)) + 1;
    x->height = max(height(x->left), height(x->right)) + 1;
 
    // Return new root
    return x;
}

/*
 * left_rotate
 *   A utility function to left rotate subtree rooted with x
 */
itnode_struct* left_rotate(itnode_struct* x)
{
  itnode_struct* y = x->right;
  itnode_struct* T2 = y->left;
 
  // Perform rotation
  y->left = x;
  x->right = T2;

  //  Update heights
  x->height = max(height(x->left), height(x->right)) + 1;
  y->height = max(height(y->left), height(y->right)) + 1;

  // Return new root
  return y;
}

/*
 * get_balance
 *   Get Balance factor of a node
 */
int get_balance(itnode_struct* root)
{
  if (root == NULL)
    return 0;

  return height(root->left) - height(root->right);
}

/*
 * insert_itnode
 *
 * @see include/annotate/itvltree.h
 */
itnode_struct* insert_itnode(itnode_struct* node, int low, int high, profile_struct_annotation* profile)
{
  // Perform the normal BST rotation
  if (node == NULL)
    return(new_itnode(low, high, profile));

  /* SEEK AND DESTROY */
  if (high < node->high)
    node->left  = insert_itnode(node->left, low, high, profile);
  else
    node->right = insert_itnode(node->right, low, high, profile);

  // Update height of this ancestor node
  node->height = max(height(node->left), height(node->right)) + 1;

  // Get the balance factor of this ancestor node to check whether
  // this node became unbalanced
  int balance = get_balance(node);
 
  // If this node becomes unbalanced, then there are 4 cases
  // Left Left Case
  if (balance > 1 && high < node->left->high)
    return right_rotate(node);

  // Right Right Case
  if (balance < -1 && high > node->right->high)
    return left_rotate(node);

  // Left Right Case
  if (balance > 1 && high > node->left->high) {
    node->left =  left_rotate(node->left);
    return right_rotate(node);
  }

  // Right Left Case
  if (balance < -1 && high < node->right->high) {
    node->right = right_rotate(node->right);
    return left_rotate(node);
  }

  // return the (unchanged) node pointer
  return node;
}

/*
 * search_itnode
 * 
 * @see include/annotate/itvltree.h
 */
void search_itnode(itnode_struct* root, int low, int high, char* annotation)
{
  // Base Case, tree is empty
  if (root == NULL)
    return;

  // If given interval overlaps with root
  // Check annotation score in order to overwrite annotation
  if (do_overlap(root->low, root->high, low, high)) {
    int score;

    if ((root->low > low) && (root->high < high))
      score = (root->high - root->low + 1) - (root->low - low) - (high - root->high);
    else if ((root->low > low) && (root->high >= high))
      score = (high - root->low + 1) - (root->low - low);
    else if ((root->low >= low) && (root->high < high))
      score = (root->high - low + 1) - (high - root->high);
    else
      score = high - low + 1;

    if (score > root->profile->anscore) { 
      strncpy(root->profile->annotation, annotation, MAX_FEATURE);
      root->profile->anscore = score;
    }
  }

  // If left child of root is present and max of left child is
  // greater than or equal to given interval, then i may
  // overlap with an interval is left subtree
  if (low < root->low)
    search_itnode(root->left, low, high, annotation);

  // Else interval can only overlap with right subtree
  if (high > root->high)
    search_itnode(root->right, low, high, annotation);
}

/*
 * destroy_itnode
 *
 * @see include/annotate/itvltree.h
 */
void destroy_itnode(itnode_struct* root)
{
  // Left child
  if (root->left != NULL) 
    destroy_itnode(root->left);
  
  // Recursive case 2
  if (root->right != NULL)
    destroy_itnode(root->right);

  // Base case
  free(root);
  root = NULL;
}
