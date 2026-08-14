/*************************************************************************
    > File Name: kqueue.c
    > Author: xlzh
    > Mail: xiaolongzhang2015@163.com
    > Created Time: 2022年09月08日 星期四 16时02分22秒
 ************************************************************************/

#include "kqueue.h"


/*! @macro
  * @abstract  whether the queue is empty
  @param  queue      the pointer to the queue object
  @return            empty: 1  otherwise: 0
 */
#define kqueue_empty(kqueue) ((kqueue)->n_item > 0 ? 0 : 1)


/*! @macro
  * @abstract  whether the queue is full
  @param  queue      the pointer to the queue object
  @return            full: 1  otherwise: 0
 */
#define kqueue_full(kqueue) ((kqueue)->n_item+1 == (kqueue)->capacity ? 1 : 0)



kqueue_t *kqueue_init(int size)
{
    kqueue_t *kqueue;

    if (size <= 0) {
        fprintf(stderr, "[SysErr:kqueue_init] the size (%d) of queue must be larger than 0!", size);
        exit(-1);
    }

    kqueue = (kqueue_t *) calloc(1, sizeof(kqueue_t));
    if (kqueue == NULL) {
        fprintf(stderr, "[SysErr:kqueue_init] failed to allocate memory for the queue!\n");
        exit(-1);
    }

    kqueue->capacity = size + 1;
    kqueue->items = (void **) calloc(kqueue->capacity, sizeof(void *));
    if (kqueue->items == NULL) {
        fprintf(stderr, "[SysErr:kqueue_init] failed to allocate memory for the items of the queue!\n");
        exit(-1);
    }

    pthread_mutex_init(&kqueue->lock, NULL);
    pthread_cond_init(&kqueue->not_empty, NULL);
    pthread_cond_init(&kqueue->not_full, NULL);

    return kqueue;
}


void kqueue_get_access(kqueue_t *kqueue)
{
    pthread_mutex_lock(&kqueue->lock);

    while (kqueue_full(kqueue) && !kqueue->is_finish) {
        pthread_cond_wait(&kqueue->not_full, &kqueue->lock);
    }
    pthread_mutex_unlock(&kqueue->lock);
}


int kqueue_push(kqueue_t *kqueue, void *data)
{
    pthread_mutex_lock(&kqueue->lock);

    /* need to call kqueue_get_access to check whether there has space */
    kqueue->items[kqueue->rear] = data;
    kqueue->rear = (kqueue->rear + 1) % kqueue->capacity;
    kqueue->n_item++;

    /* send signal to kqueue_pop that there has available data */
    pthread_cond_signal(&kqueue->not_empty);
    pthread_mutex_unlock(&kqueue->lock);

    return 0;
}


void *kqueue_get_front(kqueue_t *kqueue)
{
    void *front_data;
    pthread_mutex_lock(&kqueue->lock);

    while (kqueue_empty(kqueue) && !kqueue->is_finish) {
        pthread_cond_wait(&kqueue->not_empty, &kqueue->lock);
    }
    front_data = kqueue_empty(kqueue) ? NULL : kqueue->items[kqueue->front];
    pthread_mutex_unlock(&kqueue->lock);

    return front_data;
}


void *kqueue_pop(kqueue_t *kqueue)
{
    void *front_data;
    pthread_mutex_lock(&kqueue->lock);

    /* need to call kqueue_get_front to check whether there has items */
    front_data = kqueue->items[kqueue->front];
    kqueue->front = (kqueue->front + 1) % kqueue->capacity;
    kqueue->n_item--;

    /* send signal to kqueue_push that there has available space */
    pthread_cond_signal(&kqueue->not_full);
    pthread_mutex_unlock(&kqueue->lock);

    return front_data;
}


void kqueue_set_finish(kqueue_t *kqueue)
{
    pthread_mutex_lock(&kqueue->lock);
    kqueue->is_finish = 1;
    pthread_cond_signal(&kqueue->not_empty);
    pthread_mutex_unlock(&kqueue->lock);
}


void kqueue_set_restart(kqueue_t *kqueue)
{
    pthread_mutex_lock(&kqueue->lock);

    if (kqueue->n_item > 0) {
        fprintf(stderr, "[SysErr:kqueue_set_restart] the tasks in the queue are not finished!\n");
        exit(-1);
    }

    kqueue->front = kqueue->rear = 0;
    kqueue->n_item = 0;
    kqueue->is_finish = 0;

    pthread_mutex_unlock(&kqueue->lock);
}


void kqueue_destroy(kqueue_t *kqueue)
{
    if (kqueue == NULL)  /* the kqueue has been released */
        return ;

    if (kqueue->items == NULL)  /* the items of the queue has been released */
        return ;

    free(kqueue->items);
    free(kqueue);
}


void kqueue_foreach(kqueue_t *kqueue, kqueue_access_func func, void *func_args)
{
    pthread_mutex_lock(&kqueue->lock);

    int front = kqueue->front;
    int n_item = kqueue->n_item;

    for (int i=0; i < n_item; i++) {
        if (func(kqueue->items[front], func_args) < 0)
            break;

        front = (front + 1) % kqueue->capacity;
    }

    pthread_mutex_unlock(&kqueue->lock);
}
