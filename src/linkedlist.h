typedef struct Node {
  double val;
  int col;
  struct Node *next;
} Node;

void add_to_val(Node *head, int col, double val);
bool is_col_in_list(Node *head, int col);
void free_list(Node *head);
void add_to_list(Node **head, int col, double val);
