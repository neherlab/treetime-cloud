import pika

from .types import PikaConsumerCallable


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
