import json
import sys

from typing import Callable, Dict

import pika
from pika.adapters.blocking_connection import BlockingChannel
from pika.spec import Basic, BasicProperties

from file_store import FileStore


class Task:
  """
  Describes a task
  """

  task_id: str
  input_filenames: Dict[str, str]


# Type alias for raw Pika consumer callback
PikaConsumerCallable = Callable[
    [BlockingChannel, Basic.Deliver, BasicProperties, bytes], str]

# Type alias for Task consumer callback
TaskConsumerCallable = Callable[[Task], None]


class MessageQueueClient:
  """
  Implements a message queue client.
  Connects to the message queue service, retrieves and consumes messages using
  a given consumer.
  """

  def __init__(self, host: str, queue: str, consumer: PikaConsumerCallable):
    super().__init__()
    self._connection = pika.BlockingConnection(pika.ConnectionParameters(host))
    self._channel = self._connection.channel()
    self._channel.queue_declare(queue)
    self._channel.basic_consume(queue,
                                on_message_callback=consumer,
                                auto_ack=True)

  def start_consuming(self) -> None:
    try:
      print(" [*] Waiting for messages. To exit press CTRL+C")
      self._channel.start_consuming()
    except KeyboardInterrupt:
      self._channel.stop_consuming()
    finally:
      self._connection.close()


class TaskConsumerAdapter(PikaConsumerCallable):
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
      print(type(input_filenames))
    except KeyError as error:
      # TODO: report failure back to api
      sys.stderr.write(f"Error: task {task_id}: data is corrupt. Aborting.\n")
      raise error

    return Task(task_id, input_filenames)

  def __call__(self, _1: pika.BlockingConnection, _2: pika.spec.Basic.Deliver,
               _3: pika.spec.BasicProperties, body: bytes) -> None:
    """
    Implements call interface, so that we can pass this object as a callblack
    to Pika.
    Receives raw Pika message, converts it into a Task and consumes this task
    using the given consumer
    """

    # TODO: introduce a schema, proper deserialization and validation
    task = self._deserialize_task(body)
    self._consumer(task)


class TaskConsumer(TaskConsumerCallable):
  """ Consumes Task objects """

  def __init__(self, file_store):
    self._file_store = file_store

  def __call__(self, task: Task) -> None:
    """
    Implements call interface, so that we can pass this object as a callblack
    """
    self._consume_task(task)

  def _consume_task(self, task: Task) -> None:
    """
    Consumes a task by running a useful job on it
    """
    task_id = task.task_id
    input_filenames = task.input_filenames

    for _, filename in input_filenames.items():
      self._file_store.download_input_file(task_id, filename)

    # TODO: run the actual algorithm
    # root = ".data/"
    # cfg = tree_time_config.treetime_webconfig
    # ttw = TreeTimeWeb(root, cfg)
    # ttw.run()

    # TODO: upload output files, not inputs
    for _, filename in input_filenames.items():
      self._file_store.upload_output_file(task_id, filename)


def main():
  file_store = FileStore(
      endpoint_url="http://treetime-dev-filestore:9000",
      aws_access_key_id="minioadmin",
      aws_secret_access_key="minioadmin",
  )

  consumer = TaskConsumerAdapter(TaskConsumer(file_store=file_store))

  queue = MessageQueueClient(host="treetime-dev-taskqueue",
                             queue="tasks",
                             consumer=consumer)

  queue.start_consuming()


if __name__ == "__main__":
  main()
