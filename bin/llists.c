#include "llists.h"

item_i * llist_add_el_i(int val_i, item_i* cur)
{
	item_i *ptr, *head ;
	head = cur ; 
	ptr = Calloc(1,item_i) ;
	ptr->val = val_i ;
	ptr->next = NULL ;
	
	if (cur != NULL)
	{
		// existing list
		while (cur->next!=NULL)
		{
			cur = cur->next ;
		}
		cur->next=ptr ;
		cur = head ;
	}
	else
	{
		cur = ptr ;
	}
	
	
	return cur ;
}

int llist_destroy_i(item_i * cur ) 
{
	item_i * del ;
	while (cur!=NULL)
	{
		del=cur->next ;
		Free(cur) ;
		cur=del;
	}
	
	return 0;
}

int llist_destroy_d(item_d * cur ) 
{
	item_d * del ;
	while (cur!=NULL)
	{
		del=cur->next ;
		Free(cur) ;
		cur=del;
	}
	return 0 ;
}

item_d * llist_add_el_d(double val_d, item_d* cur)
{
	item_d *ptr , * head;
	head = cur ;
	ptr = Calloc(1,item_d) ;
	ptr->val = val_d ;
	ptr->next = NULL ;
	if (cur != NULL)
	{
		// existing list
		while (cur->next!=NULL)
		{
			cur = cur->next ;
		}
		cur->next=ptr ;
		cur = head ;
	}
	else
	{
		cur = ptr ;
	}
	
	return cur ;
}

int llist_get_el_i(int ind, item_i *cur)
{
	item_i* ptr = cur;
	int i, val ;
	for (i=0 ; i<ind ; i++)
	{
	
		ptr = cur->next ;
		cur = ptr ;
	} 
	val = ptr->val ;
	return val ;
	
}

double llist_get_el_d(int ind, item_d *cur)
{
	item_d* ptr = cur;
	int i;
	double  val ;
	for (i=0 ; i<ind ; i++)
	{
		ptr = cur->next ;
		cur = ptr ;
	} 
	val = ptr->val ;
	return val ;
	
}

int llist_get_ind_i(int val, item_i * cur)
{
	item_i * ptr = cur ;
	int i=0 ;
	while (ptr!=NULL)
	{
		if (ptr->val == val) 
		{
			ptr=NULL ;
		}
		else
		{
			i++ ;
			ptr=ptr->next ;
		}
	}
	return i ;
}

item_i *llist_delete_el_i(int ind, item_i *cur)
{
	item_i *ptr = cur, *del, *head = ptr;
	int i, val ;
	
	// if ind > 0
	for (i=0 ; i<(ind-1) ; i++)
	{
		ptr = cur->next ;
		cur = ptr ;
	} 
	
	// ptr points to entry before ind
	if (ind > 0)
	{
		del=ptr->next ;
		ptr->next = del->next ;
	}
	else
	{
		head = ptr->next ;
		del = ptr ;
	}
	
	Free(del) ;
	del = NULL ;
	return head ;
}

item_d *llist_delete_el_d(int ind, item_d *cur)
{
	item_d* ptr = cur, *del, *head = ptr;
	int i;
	
	// if ind > 0
	for (i=0 ; i<(ind-1) ; i++)
	{
		ptr = cur->next ;
		cur = ptr ;
	} 
	
	// ptr points to entry before ind
	if (ind > 0)
	{
		del=ptr->next ;
		ptr->next = del->next ;
	}
	else
	{
		head = ptr->next ;
		del = ptr ;
	}
	
	Free(del) ;
	del = NULL ;
	return head ;
}

int llist_min(item_i *cur)
{
	item_i * ptr = cur ;
	int m ;
	
	m=(ptr==NULL) ? 0 : cur->val ;
	
	while (ptr!=NULL)
	{
		if (ptr->val< m) {m=ptr->val ;} 
		ptr = ptr->next ;
		 
	}
	return m ;
}

int llist_print(item_i *cur)
{
	item_i * ptr = cur ;
	int m ;
	
	m=(ptr==NULL) ? 0 : cur->val ;
	
	while (ptr!=NULL)
	{
		Rprintf("%i, ", ptr->val) ; 
		ptr = ptr->next ;
		 
	}
	return 0 ;
}
