#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>

struct list_el_i {
	int val ;
	struct list_el_i * next ;
} ;

typedef struct list_el_i item_i ;

struct list_el_d {
	double val ;
	struct list_el_d * next ;
} ;

typedef struct list_el_d item_d ;
item_i * llist_add_el_i(int val_i, item_i* cur) ;
item_d * llist_add_el_d(double val_d, item_d* cur) ;
int llist_get_el_i(int ind, item_i *cur) ;
double llist_get_el_d(int ind, item_d *cur) ;
item_i *llist_delete_el_i(int ind, item_i *cur) ;
item_d *llist_delete_el_d(int ind, item_d *cur) ;
int llist_get_ind_i(int val, item_i * cur) ;
int llist_min(item_i *cur) ;
int llist_destroy_i(item_i * cur ) ;
int llist_destroy_d(item_d * cur ) ;
int llist_print(item_i *cur) ;

