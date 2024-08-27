/*linkedlist.c -- linked list to assemble arrays
 *
 * Written on Thursday, 15 August 2024.
 */
#include "linkedlist.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

bool is_col_in_list(Node *head, int col) {
  Node *current = head;
  while (current != NULL) {
    if (current->col == col) {
      return true;
    }
    current = current->next;
  }
  return false;
}

void free_list(Node *head) {
  Node *current = head;
  Node *next_node;
  while (current != NULL) {
    next_node = current->next;
    free(current);
    current = next_node;
  }
}

void add_to_val(Node *head, int col, double val) {
  Node *current = head;
  while (current != NULL) {
    if (current->col == col) {
      current->val += val;
      break;
    }
    current = current->next;
  }
}

void add_to_list(Node **head, int col, double val) {
  Node *new_node = (Node *)malloc(sizeof(Node));
  new_node->val = val;
  new_node->col = col;
  new_node->next = *head;
  *head = new_node;
}
