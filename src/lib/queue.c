#include "queue.h"

#include <assert.h>
#include <stdlib.h>

struct queue make_queue(void) {
    struct queue q;
    struct queue_node *dummy =
        (struct queue_node *)malloc(sizeof(struct queue_node));
    dummy->data = NULL;
    dummy->next = NULL;
    q.front = dummy;
    q.back = dummy;
    return q;
}

void free_queue(struct queue *q) {
    while (!queue_empty(q)) {
        dequeue(q);
    }
    free(q->front);
    q->front = NULL;
    q->back = NULL;
}

int queue_empty(const struct queue *q) { return q->front == q->back; }

void enqueue(struct queue *q, void *entry) {
    q->back->data = entry;
    q->back->next = (struct queue_node *)malloc(sizeof(struct queue_node));
}

void *dequeue(struct queue *q) {
    assert(!queue_empty(q));
    struct queue_node *old_node = q->front;
    void *data = old_node->data;
    q->front = old_node->next;
    free(old_node);
    return data;
}