#include <stdio.h>
#include <stdbool.h>
#include <assert.h>

#include "test.h"
#include "itree/iv.h"
#include "itree/search.h"

bool test_all()
{
  return test_iv() && test_search();
}



// ------------------------------------------------------------------
// Test code
// ------------------------------------------------------------------

bool test_iv()
{
  IV *iv = init_IV(2);
  Interval *a = init_Interval(12, 22);
  Interval *b = init_Interval(23, 33);
  Interval *c = init_Interval(34, 44);

  printf("iv: vector realloc when size is passed\n");
  add_IV(iv, *a);
  add_IV(iv, *b);
  assert(iv->available == 2);
  add_IV(iv, *c);
  assert(iv->available == 4);

  printf("iv: correct addition of elements\n");
  assert(iv->v[0].start == 12);
  assert(iv->v[1].start == 23);
  assert(iv->v[2].start == 34);

  free(a);
  free(b);
  free(c);
  free_IV(iv);

  return true;
}

bool test_search()
{
  IntervalResult *res;
  IA *ia = init_set_IA(2);
  Interval a = {.start = 10,.stop = 20 };
  Interval b = {.start = 30,.stop = 40 };
  ia->v[0] = a;
  ia->v[1] = b;
  struct IntervalTree *tree = build_tree(ia);

  printf("search: count point overlaps in itree\n");
  assert(count_point_overlaps(5, tree) == 0);
  assert(count_point_overlaps(25, tree) == 0);
  assert(count_point_overlaps(45, tree) == 0);
  assert(count_point_overlaps(10, tree) == 1);
  assert(count_point_overlaps(11, tree) == 1);
  assert(count_point_overlaps(20, tree) == 1);

  printf("search: retrieve intervals overlapping itree\n");
  res = get_point_overlaps(5, tree);
  assert(res->inbetween == true);
  free_IntervalResult(res);

  res = get_point_overlaps(12, tree);
  assert(res->inbetween == false);
  assert(res->iv->v[0].start == 10);
  free_IntervalResult(res);

  res = get_point_overlaps(25, tree);
  assert(res->inbetween == true);
  free_IntervalResult(res);

  Interval c = {.start = 1,.stop = 5 };
  Interval d = {.start = 5,.stop = 15 };
  Interval e = {.start = 22,.stop = 28 };

  printf("search: count interval overlaps in itree\n");
  assert(count_interval_overlaps(&c, tree) == 0);
  assert(count_interval_overlaps(&d, tree) == 1);
  assert(count_interval_overlaps(&e, tree) == 0);

  printf("search: retrieve intervals overlapping itree\n");
  res = get_interval_overlaps(&c, tree);
  assert(res->inbetween == true);
  free_IntervalResult(res);

  res = get_interval_overlaps(&d, tree);
  assert(res->inbetween == false);
  assert(res->iv->v[0].start == 10);
  free_IntervalResult(res);

  res = get_interval_overlaps(&e, tree);
  assert(res->inbetween == true);
  free_IntervalResult(res);

  free_IntervalTree(tree);
  return true;
}
