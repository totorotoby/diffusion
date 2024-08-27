/*mesh.c -- Getting mesh from Tmesh
 *
 * Written on Friday, 26 July 2024.
 */
#include "mesh.h"
#include "matrix.h"
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BUFFER_SIZE 1024
#define NV_PER_ELEM 3

int min(const int a, const int b) { return (a < b) ? a : b; }

int max(const int a, const int b) { return (a > b) ? a : b; }

int get_mesh(const char *filename, int *num_nodes, int *num_elements,
             double **vx, double **vy, int **EToV) {

  FILE *mesh_file = fopen(filename, "r");
  if (!mesh_file) {
    perror("Failed to open file");
    return 1;
  }

  char buffer[BUFFER_SIZE];
  long offset = 0;
  size_t len = strlen("NODES");
  char found = 0;

  while (fgets(buffer, BUFFER_SIZE, mesh_file)) {
    char *pos = strstr(buffer, "NODES");
    if (pos) {
      offset += (pos - buffer) + len;
      found = 1;
      break;
    }
    offset += strlen(buffer);
  }

  if (found) {
    fseek(mesh_file, offset, SEEK_SET);
    fscanf(mesh_file, "%d", num_nodes);

    *vx = (double *)malloc(sizeof(double) * *num_nodes);
    *vy = (double *)malloc(sizeof(double) * *num_nodes);
    int n = 0;

    double xcoord;
    double ycoord;
    while (n < *num_nodes - 1) {
      fscanf(mesh_file, "%d", &n);
      fscanf(mesh_file, "%lf %lf", *vx + n, *vy + n);
    }
    fscanf(mesh_file, "%*s %d", num_elements);
    *EToV = (int *)malloc(sizeof(int) * *num_elements * NV_PER_ELEM);
    n = 0;
    while (n < *num_elements - 1) {
      fscanf(mesh_file, "%d", &n);
      // if ever elements are different than triangles
      // this needs to be changed to a loop
      fscanf(mesh_file, "%d %d %d", *EToV + NV_PER_ELEM * n,
             *EToV + NV_PER_ELEM * n + 1, *EToV + NV_PER_ELEM * n + 2);
    }
  } else {
    printf("Something is wrong with mesh file format\n");
    return 1;
  }

  // write_EToV_array(*EToV, *num_elements, "plotting/EToV_array");

  fclose(mesh_file);

  return 0;
}

// may only work for triangular meshes
int get_mesh_degrees(const int *EToV, int **degree, int *num_edges,
                     const int num_elements, const int num_nodes) {

  *degree = (int *)calloc(num_nodes, sizeof(int));
  // could be parallelized
#pragma omp parallel for
  for (int i = 0; i < NV_PER_ELEM * num_elements; i++) {
    int vertex = EToV[i];
#pragma omp atomic
    (*(*degree + vertex))++;
  }

  int sum = 0;
#pragma omp parallel for reduction(+ : sum)
  for (int v = 0; v < num_nodes; v++) {
    sum += (*degree)[v];
  }

  *num_edges = (sum + num_nodes) / 2;

  return 0;
}

int tri_mesh_csr_adj_matrix(csr **adj_matrix, int **boundary_nodes,
                            int *num_boundary_nodes, const int *EToV,
                            const int num_elements, const int num_edges,
                            int num_nodes) {

  int error;
  error = allocate_csr_novals(adj_matrix, 2 * num_edges, num_nodes, num_nodes);

  // hashing to check if edge has already been added using a * N + b, where a =
  // min(v, u), and b = max(v, u) this probably won't work if we get too big...
  // O(n^2) storage... can convert to linked-list or hashed array of linked list
  int *edge_check = (int *)calloc(num_nodes * num_nodes, sizeof(int));

  for (int e = 0; e < num_elements; e++) {
    int v1 = EToV[3 * e];
    int v2 = EToV[3 * e + 1];
    int v3 = EToV[3 * e + 2];

    int edges[3][2] = {{v1, v2}, {v2, v3}, {v3, v1}};

    for (int edge = 0; edge < 3; edge++) {
      int a = min(edges[edge][0], edges[edge][1]);
      int b = max(edges[edge][0], edges[edge][1]);
      // problem is here
      if (!edge_check[a * num_nodes + b]) {
        (*adj_matrix)->row_ptr[a + 1]++;
        (*adj_matrix)->row_ptr[b + 1]++;
      }
      edge_check[a * num_nodes + b]++;
    }
  }

  // getting number of boundary nodes
  *num_boundary_nodes = 0;
  for (int e = 0; e < num_elements; e++) {
    int v1 = EToV[3 * e];
    int v2 = EToV[3 * e + 1];
    int v3 = EToV[3 * e + 2];

    int edges[3][2] = {{v1, v2}, {v2, v3}, {v3, v1}};

    for (int edge = 0; edge < 3; edge++) {
      int a = min(edges[edge][0], edges[edge][1]);
      int b = max(edges[edge][0], edges[edge][1]);
      if (edge_check[a * num_nodes + b] == 1) {
        (*num_boundary_nodes)++;
      }
    }
  }

  printf("# boundary nodes: %d\n", *num_boundary_nodes);
  *boundary_nodes = (int *)malloc(sizeof(int) * *num_boundary_nodes);

  // convert row_ptr from histogram to cumulative sum
  for (int v = 1; v <= num_nodes; v++) {
    (*adj_matrix)->row_ptr[v] += (*adj_matrix)->row_ptr[v - 1];
  }

  int *pos = (int *)malloc(num_nodes * sizeof(int));
  for (int v = 0; v < num_nodes; v++) {
    pos[v] = (*adj_matrix)->row_ptr[v];
  }

  int b_index = 0;
  for (int e = 0; e < num_elements; e++) {
    int v1 = EToV[3 * e];
    int v2 = EToV[3 * e + 1];
    int v3 = EToV[3 * e + 2];

    int edges[3][2] = {{v1, v2}, {v2, v3}, {v3, v1}};

    for (int edge = 0; edge < 3; edge++) {
      int a = min(edges[edge][0], edges[edge][1]);
      int b = max(edges[edge][0], edges[edge][1]);

      if (edge_check[a * num_nodes + b] == 1) {
        int found1 = 0;
        int found2 = 0;
        for (int vi = 0; vi < *num_boundary_nodes; vi++) {
          if (a == (*boundary_nodes)[vi])
            found1 = 1;
          if (b == (*boundary_nodes)[vi])
            found2 = 1;
        }
        if (!found1)
          (*boundary_nodes)[b_index++] = a;
        if (!found2)
          (*boundary_nodes)[b_index++] = b;
      }
      if (edge_check[a * num_nodes + b]) {
        (*adj_matrix)->cols[pos[a]++] = b;
        (*adj_matrix)->cols[pos[b]++] = a;
        edge_check[a * num_nodes + b] = 0;
      }
    }
  }

  write_int_array((*boundary_nodes), (*num_boundary_nodes),
                  "plotting/boundary_nodes");

  free(pos);
  free(edge_check);
  return error;
}

int degree_comp(const void *v1, const void *v2, const void *degree) {
  return ((int *)degree)[*((int *)v1)] - ((int *)degree)[*((int *)v2)];
}

int cuthill_mckee(const csr *adj_matrix, const int *degree, int **EToV,
                  double **vx, double **vy, int **boundary_nodes,
                  const int num_boundary_nodes, const int num_nodes,
                  const int num_elements) {

  // initalize permutation matrix
  int *perm = (int *)malloc(num_nodes * sizeof(int));
  for (int v = 0; v < num_nodes; v++)
    perm[v] = -1;

  // getting a starting node
  int current_degree = INT_MAX;
  int current;
  int next;

  for (int v = 0; v < num_nodes; v++) {
    if (degree[v] <= current_degree) {
      current = v;
      current_degree = degree[v];
    }
  }
  perm[current] = 0;

  int *node_ptr = adj_matrix->row_ptr;
  int *neighbors = adj_matrix->cols;
  // queue to hold order to visit nodes
  int *queue = (int *)malloc(num_nodes * sizeof(int));
  queue[0] = current;
  int front = 0;
  int back = 1;
  int visited = 1;

  while (front < back) {

    int current = queue[front++];
    int begin = node_ptr[current];
    int end = node_ptr[current + 1];

    // sort by degree before adding to permutation matrix
    qsort_r(neighbors + begin, end - begin, sizeof(int), (void *)degree_comp,
            degree);

    for (int i = begin; i < end; i++) {
      if (perm[neighbors[i]] == -1) {
        queue[back++] = neighbors[i];
        perm[neighbors[i]] = visited;
        visited++;
      }
    }
  }

  // permuting orginal arrays
  double *nvx = (double *)malloc(num_nodes * sizeof(double));
  double *nvy = (double *)malloc(num_nodes * sizeof(double));
  int *nboundary_nodes = (int *)malloc(num_boundary_nodes * sizeof(int));
  int *nEToV = (int *)malloc(NV_PER_ELEM * num_elements * sizeof(int));

  for (int v = 0; v < num_boundary_nodes; v++)
    nboundary_nodes[v] = perm[(*boundary_nodes)[v]];

#pragma omp parallel for
  for (int v = 0; v < num_nodes; v++) {
    nvx[perm[v]] = (*vx)[v];
    nvy[perm[v]] = (*vy)[v];
  }

#pragma omp parallel for
  for (int e = 0; e < num_elements; e++) {
    for (int v = 0; v < NV_PER_ELEM; v++) {
      nEToV[3 * e + v] = perm[(*EToV)[3 * e + v]];
    }
  }

  write_csr_novals(adj_matrix, "plotting/test_row_ptr", "plotting/test_col");
  write_int_array(perm, num_nodes, "plotting/perm");
  write_double_array(*vx, num_nodes, "plotting/vx");
  write_double_array(*vy, num_nodes, "plotting/vy");

  free(*EToV);
  free(*vx);
  free(*vy);
  free(*boundary_nodes);

  *EToV = nEToV;
  *vx = nvx;
  *vy = nvy;
  *boundary_nodes = nboundary_nodes;

  free(queue);
  free(perm);

  return 0;
}
