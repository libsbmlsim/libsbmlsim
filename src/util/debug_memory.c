/**
 * <!--------------------------------------------------------------------------
 * This file is part of libSBMLSim.  Please visit
 * http://fun.bio.keio.ac.jp/software/libsbmlsim/ for more
 * information about libSBMLSim and its latest version.
 *
 * Copyright (C) 2011-2017 by the Keio University, Yokohama, Japan
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation.  A copy of the license agreement is provided
 * in the file named "LICENSE.txt" included with this software distribution.
 * ---------------------------------------------------------------------- -->*/
#define C_DEBUG_MEMORY_C  /* Have to be defined before including debug_memory.h */

#include "debug_memory.h"

debug_node_t* create_node(Site* site) {
#ifdef DEBUG_MEMORY_DEBUG
  printf("create_node for %s %4d [%p]\n", site->s.file, site->s.line, site);
#endif
  debug_node_t* head = malloc(sizeof(debug_node_t));
  head->site = site;
  head->next = NULL;
  return head;
}

void add_node(Site* site) {
#ifdef DEBUG_MEMORY_DEBUG
  printf("   add_node for %s %4d [%p]\n", site->s.file, site->s.line, site);
#endif
  debug_node_t* current;
  if (debug_root_node == NULL) {
    debug_root_node = create_node(site);
  } else {
    current = debug_root_node;
    while (current->next != NULL) {
      current = current->next;
    }
    current->next = create_node(site);
  }
}

void remove_node(Site* site) {
#ifdef DEBUG_MEMORY_DEBUG
  printf("remove_node for %s %4d [%p]\n", site->s.file, site->s.line, site);
#endif
  debug_node_t* tmp_node;
  debug_node_t* current;
  if (debug_root_node->site == site) {
    tmp_node = debug_root_node;
    debug_root_node = debug_root_node->next;
    free(tmp_node);
    return;
  }
  current = debug_root_node;
  while (current->next != NULL) {
    if (current->next->site == site) {
      tmp_node = current->next;
      current->next = tmp_node->next;
      free(tmp_node);
      return;
    }
    current = current->next;
  }
}

void print_debug_node(debug_node_t* node) {
  printf("address : %p\n", node->site);
  printf("size    : %zu bytes\n", node->site->s.n);
  printf("file    : %s\n", node->site->s.file);
  printf("line    : %d\n", node->site->s.line);
  printf("----------------------------------------\n");
}

void print_allocated_memory(void) {
  debug_node_t *current = debug_root_node;
  if (current == NULL) return;
  printf("=== Allocated Memory (v%s) ==========\n", DEBUG_MEMORY_DOTTED_VERSION);
  while (current != NULL) {
    print_debug_node(current);
    current = current->next;
  }
  printf("Total   : %d bytes\n", total_allocated);
  printf("========================================\n");
}

void* debug_malloc(size_t n, char *file, int line) {
  char *rp;
  rp = (char*)malloc(sizeof(Site)+n);
  total_allocated += n;
  ((Site*)rp)->s.n = n;
  ((Site*)rp)->s.file = file;
  ((Site*)rp)->s.line = line;
  add_node((Site*)rp);
  return (void*)(rp + sizeof(Site));
}

void* debug_calloc(size_t c, size_t n, char *file, int line) {
  char *rp;
  rp = (char*)calloc(sizeof(Site)+n*c, sizeof(char));
  total_allocated += n*c;
  ((Site*)rp)->s.n = n*c;
  ((Site*)rp)->s.file = file;
  ((Site*)rp)->s.line = line;
  add_node((Site*)rp);
  return (void*)(rp + sizeof(Site));
}

void* debug_realloc(void *ptr, size_t n, char *file, int line) {
  char *rp;
  char *tmp;
  tmp = ((char*)ptr) - sizeof(Site);
  total_allocated -= ((Site*)tmp)->s.n;
  rp = (char*)realloc(tmp, sizeof(Site)+n);
  total_allocated += n;
  ((Site*)rp)->s.n = n;
  ((Site*)rp)->s.file = file;
  ((Site*)rp)->s.line = line;
  add_node((Site*)rp);
  return (void*)(rp + sizeof(Site));
}

void debug_free(void *p, char *file, int line) {
  char *rp;
  rp = ((char*)p) - sizeof(Site);
  total_allocated -= ((Site*)rp)->s.n;
  ((Site*)rp)->s.file = file;
  ((Site*)rp)->s.line = line;
  remove_node((Site*)rp);
  free(rp);
}

char* debug_strdup(const char* str, char *file, int line) {
  char *rp = NULL;
  if (str) {
    size_t n = strlen(str) + 1;
    rp = malloc(sizeof(Site)+n);
    if (rp) {
#ifdef _MSC_VER
      strcpy_s(rp + sizeof(Site), n, str);
#else
      strcpy(rp + sizeof(Site), str);
#endif
      total_allocated += n;
      ((Site*)rp)->s.n = n;
      ((Site*)rp)->s.file = file;
      ((Site*)rp)->s.line = line;
      add_node((Site*)rp);
    }
  }
  return rp + sizeof(Site);
}
