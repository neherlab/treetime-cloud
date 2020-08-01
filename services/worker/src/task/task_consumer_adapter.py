import json
import sys

from pika.adapters.blocking_connection import BlockingChannel
from pika.spec import Basic, BasicProperties

from .types import Task, TaskConsumerCallable, TaskResult


class TaskConsumerAdapter:
  """
  Wraps Task consumer, allowing to consume raw messages from Pika queue.
  Deserializes raw messages into Tasks, validates and consumes these tasks
  using the wrapper consumer callable.
  """

  def __init__(self, consumer: TaskConsumerCallable):
    super().__init__()
    self._consumer = consumer

  def _deserialize_task(self, body: bytes) -> Task:
    """
    Deserializes task from raw bytes into a Task object
    """
    try:
      task = json.loads(body)["data"]
      task_id = task["taskId"]
      input_filenames = task["inputFilenames"]
    except KeyError as error:
      # TODO: report failure back to api
      sys.stderr.write(f"Error: task {task_id}: data is corrupt. Aborting.\n")
      raise error

    return Task(task_id, input_filenames)

  def __call__(self, _1: BlockingChannel, _2: Basic.Deliver,
               _3: BasicProperties, body: bytes) -> TaskResult:
    """
    Implements call interface, so that we can pass this object as a callblack
    to Pika.
    Receives raw Pika message, converts it into a Task and consumes this task
    using the given consumer
    """

    # TODO: introduce a schema, proper deserialization and validation
    task = self._deserialize_task(body)
    return self._consumer(task)
