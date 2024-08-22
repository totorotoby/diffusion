/*linkedlist.c -- linked list to assemble arrays
 *
 * Written on Thursday, 15 August 2024.
  */
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include "linkedlist.h"

bool is_col_in_list(Node *head, int col)
{
  Node *current = head;
  while (current != NULL)
    {
      if (current->col == col)
	{
	  return true;
	}
      current = current->next;
    }
  return false;
}

void free_list(Node *head)
{
  Node *current = head;
  Node *next_node;
  while (current != NULL)
    {
      next_node = current->next;
      free(current);            
      current = next_node;  
    }
}

void add_to_val(Node *head, int col, double mval, double sval)
{
  Node *current = head;
  while (current != NULL)
    {
      if (current->col == col)
	{
	  current->mval += mval;
	  current->sval += sval;
	  break;
	}
      current = current->next;
    }
}

void add_to_list(Node **head, int col, double mval, double sval)
{
  Node *new_node = (Node *)malloc(sizeof(Node));
  new_node->mval = mval;
  new_node->sval = sval;
  new_node->col = col;
  new_node->next = *head;
  *head = new_node;
}
