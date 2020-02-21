# treetime/taskqueue

This module implements multi-producer-multii-consumer message queue that
serializes incoming tasks in order to balance the load across workers.

Task producer (api) pushes tasks into the queue. Task consumers (workers) pop 
tasks from the queue.

This may be used as is (e.g. locally, in dev mode) or may be replaced with
AWS Batch or other cloud solutions which include message queues.
