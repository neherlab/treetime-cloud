import json
import os
import requests
from pprint import pprint

import pika

from .types import PikaConsumerCallable

API_HOST = os.environ['API_HOST']
API_PORT = os.environ['API_PORT']
API_ADDRESS = f"{API_HOST}:{API_PORT}"

class MessageQueueClient:
  """
  Implements a message queue client.
  Connects to the message queue service, retrieves and consumes messages using
  a given consumer.
  """

  def __init__(self, host: str, queue: str, consumer: PikaConsumerCallable):
    super().__init__()
    self._consumer = consumer
    self._connection_tasks = pika.BlockingConnection(pika.ConnectionParameters(host))
    self._channel_tasks = self._connection_tasks.channel()
    self._channel_tasks.queue_declare(queue, durable=False)
    self._channel_tasks.basic_consume(
        queue,
        on_message_callback=self.on_message,
        auto_ack=True,
    )

  def on_message(self, channel, method, properties, body):
    result = self._consumer(channel, method, properties, body)
    data = dict(result._asdict())
    data = json.loads(json.dumps(data))
    # print(type(data))
    pprint(data)
    requests.post(f"http://{API_ADDRESS}/api/v1/taskResults", data=data)

  def start_consuming(self) -> None:
    try:
      print(" [*] Waiting for messages. To exit press CTRL+C")
      self._channel_tasks.start_consuming()
    except KeyboardInterrupt:
      self._channel_tasks.stop_consuming()
    finally:
      self._connection_tasks.close()
