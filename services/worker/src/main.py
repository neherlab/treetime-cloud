import os

from filestore import FileStore
from task import MessageQueueClient, TaskConsumer, TaskConsumerAdapter

FILESTORE_HOST = os.environ['FILESTORE_HOST']
FILESTORE_PORT = os.environ['FILESTORE_PORT']
FILESTORE_ADDRESS = f'{FILESTORE_HOST}:{FILESTORE_PORT}'

TASK_QUEUE_HOST = os.environ['TASK_QUEUE_HOST']


def main() -> None:
  file_store = FileStore(
      endpoint_url=FILESTORE_ADDRESS,
      aws_access_key_id="minioadmin",
      aws_secret_access_key="minioadmin",
  )

  consumer = TaskConsumerAdapter(TaskConsumer(file_store=file_store))

  queue = MessageQueueClient(
      host=TASK_QUEUE_HOST,
      queue="tasks",
      consumer=consumer,
  )

  queue.start_consuming()


if __name__ == "__main__":
  main()
