import json
import os
import sys

from botocore.client import Config
import boto3
import pika

import tree_time_config


if __name__ == "__main__":
    root = "../data/"

    cfg = tree_time_config.treetime_webconfig

    connection = pika.BlockingConnection(
        pika.ConnectionParameters(host="treetime-dev-taskqueue")
    )

    channel = connection.channel()

    channel.queue_declare(queue="tasks")

    def callback(ch, method, properties, body):
        try:
            task = json.loads(body)["data"]
            task_id = task["taskId"]
            input_filepaths = task["inputFilepaths"]
            print(input_filepaths)
        except:
            # TODO: report failure back to api
            sys.stderr.write(
                f"Error: task {task_id}: data is corrupt. Aborting."
            )

        session = boto3.session.Session()
        s3 = boto3.client(
            service_name="s3",
            endpoint_url="http://treetime-dev-filestore:9000",
            aws_access_key_id="minioadmin",
            aws_secret_access_key="minioadmin",
            config=Config(signature_version="s3v4"),
        )

        def download_from_s3(path: str):
            out_path = os.path.join("..", "data", path)
            os.makedirs(os.path.dirname(out_path), exist_ok=True)
            s3.download_file("treetime", path, out_path)

        for _, filepath in input_filepaths.items():
            download_from_s3(filepath)

        # ttw = TreeTimeWeb(root, cfg)
        # ttw.run()

    channel.basic_consume(
        queue="tasks", on_message_callback=callback, auto_ack=True
    )

    print(" [*] Waiting for messages. To exit press CTRL+C")
    channel.start_consuming()
