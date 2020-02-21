from typing import Callable, Dict, NamedTuple

from pika.adapters.blocking_connection import BlockingChannel
from pika.spec import Basic, BasicProperties


class Task(NamedTuple):
  task_id: str
  input_filenames: Dict[str, str]


# Describes raw Pika consumer callback
PikaConsumerCallable = Callable[
    [BlockingChannel, Basic.Deliver, BasicProperties, bytes], None]

# Describes Task consumer callback
TaskConsumerCallable = Callable[[Task], None]
