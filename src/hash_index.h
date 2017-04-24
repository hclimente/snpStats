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

/* the maximum length of strings the look-up table */
#define MAX_ID 128

/* for MAX_ID < 4 or 8, name[] is smaller than 
   *name and also have better locality */
typedef struct node {
  struct node *next;
  char name[MAX_ID];
  int value;
} *t_node;

typedef struct index_db_struct {
  t_node *nodelist;
  int bitmask;
} *index_db;

/*
  This module provides 4 routines - creating a lookup table,
  insert a new entry, lookup an old entry, and destroy
  the lookup table.

  The create routine take a size argument, for the expected
  number of entries - does not need to be accurate (over=wastage,
  under=slower-lookup) - in the worse case of size=1, 
  look-up becomes a linear search of all the entries.
  Returns an "index_db" type for the newly created table, or NULL
  from failure (to allocate the required memory).

  The insert routine takes an index_db created earlier,
  a name and a non-negative value pair. return 0 for 
  successful, or -1 for failure. value should be non-negative.

  The lookup routine takes an index_db created earlier,
  and a name. Returns an earlier inserted value,
  or -1 if not found.

  The destroy routine takes an index_db created earlier,
  and de-allocate all resources associate with it. Returns void.

*/

index_db index_create(int size);
int index_insert(index_db db, const char *name, int value);
int index_insert_case_independent(index_db db, const char *name, int value);
int index_lookup(index_db db, const char *name);
void index_destroy(index_db db);
