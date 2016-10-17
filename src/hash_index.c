/*
 *  Copyright (C) 2006  Hin-Tak Leung
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  A copy of the GNU General Public License is available via WWW at
 *  http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
 *  writing to the Free Software Foundation, Inc., 51 Franklin Street
 *  Fifth Floor, Boston, MA 02110-1301  USA.
 */

/* 
   The design decision for this hash index implementation is
   somewhat different from the usual: all we want is for
   the look up to be very fast. We don't care about the usual obsessions 
   of hash designers about dispersion - cache locality actually works
   in our favour, as genetic data is likely to be partially sorted,
   so we *do* want similar strings to be similiarly located,
   unlike the usual design.

   For look up to be fast, 
   (1) we pick a very simple hash function,
   (2) and have all the checking during insertion,
   (3) over-allocate so that the link-lists are almost all single-member,
   (4) use bit-masking instead of division.
*/

/* MAX_HASH_SIZE is the size of the table.
   There is no reason why it can't be 2^31,
   but 2^31 is probably a bit unrealistic. 
   Just here to keep it from growing out of bound 
*/

#include <stdlib.h>
/* glibc requires this to pull in strdup */
#ifndef __USE_BSD
#define __USE_BSD
#endif

#include <string.h>
#include "hash_index.h"

#define MAX_HASH_SIZE 1000000

/* Dan J Berstein's hash function */
static unsigned int hash(const char *s, int bitmask) {
  unsigned int h;
  for (h = 5381; *s; s++) {
    h *= 33;
    h += *s;
  }
  return (h & bitmask);
}

index_db index_create(int size) {
  index_db result = (index_db)calloc(1, sizeof(struct index_db_struct));
  if (result) {
    int hash_size = 1;
    while ((hash_size < size) && (hash_size < MAX_HASH_SIZE))
      hash_size <<=1;
    result->nodelist =(t_node *)calloc(hash_size, sizeof(t_node));
    result->bitmask = hash_size - 1;
  }
  return result;
}

/* check every thing at insert, so that we don't need to 
   check at lookup */
int index_insert(index_db db, const char *name, int value) {
  if ((strlen(name) < MAX_ID) && (index_lookup(db, name) < 0) && (value >=0)) {
    t_node this_node = (t_node) calloc(1, sizeof(struct node));
    if ((this_node)  && strcpy(this_node->name, name)) {
      this_node->value = value;
      int idx = hash(name, db->bitmask);
      this_node->next = db->nodelist[idx];
      db->nodelist[idx] = this_node;
      return 0;
    }
  }

  return -1;
}

/* an case-insensitive version of insert:
   convert to all uppercase, and convert to
   all lowercase and insert both of those.
*/
int index_insert_case_independent(index_db db, const char *name, int value) {
  
  char *dup_up = strdup(name);
  char *dup_low = strdup(name);
  char *dup_up_ptr  = dup_up;
  char *dup_low_ptr = dup_low;
  /* poor man's lower case and uppercase */
  while (*dup_up) {
    *dup_up = (*dup_up) | 0x20;
    dup_up++;
  }
  while (*dup_low) {
    *dup_low = (*dup_low) & 0xDF;
    dup_low++;
  }
  /* insert the all uppercase or all lowercase versions if they are different */
  if (strcmp(name, dup_up_ptr)) {
    index_insert(db, dup_up_ptr, value);
  }
  if (strcmp(name, dup_low_ptr) && strcmp(dup_up_ptr, dup_low_ptr)) {
    index_insert(db, dup_low_ptr, value);
  }
  free(dup_up_ptr);
  free(dup_low_ptr);
  return index_insert(db, name, value);
}

int index_lookup(index_db db, const char *name) {
  t_node np;
  for(np=db->nodelist[hash(name, db->bitmask)]; np != NULL; np=np->next)
    if(!strcmp(name,np->name))
      return(np->value);

  return -1;
}

void index_destroy(index_db db) {
  if(!db)
    return;

  int i =0;
  for (i = 0; i <= db->bitmask ; i++) {
    t_node node = db->nodelist[i]; 
    while (node) {
      t_node this_node = node;
      node = node->next; /* have to copy first before destroying */
      free(this_node);
    }
  }
  free(db->nodelist);
  free(db);
}
