/*************************************************************************
    > File Name: kqueue.h
    > Author: xlzh
    > Mail: xiaolongzhang2015@163.com
    > Created Time: 2022年09月08日 星期四 15时55分12秒
 ************************************************************************/

#ifndef QUEUE_KQUEUE_H
#define QUEUE_KQUEUE_H

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>


/*! @typedef kqueue_t
 @abstract struct of the single-direction circular queue
 @field front, rear       the front and rear index of the circular queue
 @field capacity         the capacity of the queue (one more than the maximum number of items, to distinguish the full/empty state)
 @field n_item            the number of items in the queue
 @field is_finish        1: the queue is finished, no more items will be pushed
 @field items             the pointer to the item array of the queue
 @field lock              the mutex used to protect the queue
 @field not_empty        the condition variable, signaled when an item is pushed
 @field not_full         the condition variable, signaled when an item is popped
*/
typedef struct kqueue_t {
    int front;
    int rear;
    int capacity;
    int n_item;
    int is_finish;
    void **items;
    pthread_mutex_t lock;
    pthread_cond_t not_empty;
    pthread_cond_t not_full;
} kqueue_t;


/*! @function
  * @abstract  initialize the circular queue
  @param  size        the maximum number of items allowed
  @return             a circular queue object
 */
kqueue_t *kqueue_init(int size);


/*! @function
  * @abstract  get access to the kqueue before pushing an item into it
  @param  kqueue     the pointer to the queue object
  @return            void
 */
void kqueue_get_access(kqueue_t *kqueue);


/*! @function
  * @abstract  add a new item to the rear of the queue
  @param  kqueue      the pointer to the queue object
  @param  data        the user defined data field
  @return             0: success  -1:failed
 */
int kqueue_push(kqueue_t *kqueue, void *data);


/*! @function
  * @abstract  get the front data to decide whether it is needed to pop the data
  @param  kqueue     the pointer to the queue object
  @return            (void*)data: success  NULL: no more items in the queue
 */
void *kqueue_get_front(kqueue_t *kqueue);


/*! @function
  * @abstract  remove the item in the front of the queue
  @param  kqueue     the pointer to the queue object
  @return            (void*)data: success  NULL: no more items in the queue
 */
void *kqueue_pop(kqueue_t *kqueue);


/*! @function
  * @abstract  set the finish signal for the queue
  @param  kqueue     the pointer to the queue object
  @return            void
 */
void kqueue_set_finish(kqueue_t *kqueue);


/*! @function
  * @abstract  restart the queue (keep the queue memory untouched)
  @param  kqueue     the pointer to the queue object
  @return            void
 */
void kqueue_set_restart(kqueue_t *kqueue);


/*! @function
  * @abstract  callback function to access the data of the kqueue
  @param  item       the pointer to the queue data
  @param  args       the args of user provided
  @return            0: continue iterating the kqueue; other: stop the iteration and return
 */
typedef int (*kqueue_access_func)(void *item, void *args);


/*! @function
  * @abstract  iterate the items in the queue
  @param  kqueue     the pointer to the queue object
  @param  func       the callback function to operate the kqueue data
  @param  func_args  the args for the user provided function
  @return            void
 */
void kqueue_foreach(kqueue_t *kqueue, kqueue_access_func func, void *func_args);


/*! @function
  * @abstract  check whether the tasks of the queue are finished
  @param  queue      the pointer to the queue object
  @return            0:not yet; 1:finish
 */
#define kqueue_is_finish(_kqueue) ((_kqueue)->is_finish && ((_kqueue)->n_item==0))


/*! @function
  * @abstract  destroy the queue
  @param  kqueue     the pointer to the queue object
  @return            void
 */
void kqueue_destroy(kqueue_t *kqueue);


/*! @macro
  * @abstract  the number of items in the queue
  @param  kqueue     the pointer to the queue object
  @return            the number of items
 */
#define kqueue_size(_kqueue) ((_kqueue)->n_item)


#endif //QUEUE_KQUEUE_H
