from filestore import FileStore
from task import MessageQueueClient, TaskConsumer, TaskConsumerAdapter


def main() -> None:
  file_store = FileStore(
      endpoint_url="http://treetime-dev-filestore:9000",
      aws_access_key_id="minioadmin",
      aws_secret_access_key="minioadmin",
  )

  consumer = TaskConsumerAdapter(TaskConsumer(file_store=file_store))

  queue = MessageQueueClient(
      host="treetime-dev-taskqueue",
      queue="tasks",
      consumer=consumer,
  )

  queue.start_consuming()


if __name__ == "__main__":
  main()
