/**
 * \file Simple C queue.
 */

#ifndef QUEUE_H
#define QUEUE_H

struct queue_node {
    void *data;
    struct queue_node *next;
};

struct queue {
    struct queue_node *front, *back;
};

struct queue make_queue(void);

void free_queue(struct queue *q);

int queue_empty(const struct queue *q);

void enqueue(struct queue *q, void *entry);

void *dequeue(struct queue *q);

#endif  // QUEUE_H
