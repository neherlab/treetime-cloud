import json
import sys

import pika

from file_store import FileStore

import tree_time_config

if __name__ == "__main__":
  root = ".data/"

  cfg = tree_time_config.treetime_webconfig

  connection = pika.BlockingConnection(
      pika.ConnectionParameters(host="treetime-dev-taskqueue"))

  channel = connection.channel()

  channel.queue_declare(queue="tasks")

  filestore = FileStore(
      endpoint_url="http://treetime-dev-filestore:9000",
      aws_access_key_id="minioadmin",
      aws_secret_access_key="minioadmin",
  )

  def callback(_1, _2, _3, body):
    try:
      task = json.loads(body)["data"]
      task_id = task["taskId"]
      input_filenames = task["inputFilenames"]
      print(input_filenames)
    except KeyError:
      # TODO: report failure back to api
      sys.stderr.write(f"Error: task {task_id}: data is corrupt. Aborting.")

    for _, filename in input_filenames.items():
      filestore.download_input_file(task_id, filename)

    for _, filename in input_filenames.items():
      filestore.upload_output_file(task_id, filename)

    # ttw = TreeTimeWeb(root, cfg)
    # ttw.run()

  channel.basic_consume(queue="tasks",
                        on_message_callback=callback,
                        auto_ack=True)

  print(" [*] Waiting for messages. To exit press CTRL+C")
  channel.start_consuming()
