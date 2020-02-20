import os

from botocore.client import Config
import boto3


class FileStore:

  def __init__(
      self,
      endpoint_url: str,
      aws_access_key_id: str,
      aws_secret_access_key: str,
  ):
    boto3.session.Session()
    self.s3 = boto3.client(
        service_name="s3",
        endpoint_url=endpoint_url,
        aws_access_key_id=aws_access_key_id,
        aws_secret_access_key=aws_secret_access_key,
        config=Config(signature_version="s3v4"),
    )

  def download_input_file(self, prefix: str, filename: str):
    remote_path = f"{prefix}/input/{filename}"
    local_path = os.path.join(".data", prefix, "input", filename)

    os.makedirs(os.path.dirname(local_path), exist_ok=True)

    self.s3.download_file("treetime", remote_path, local_path)

  def upload_output_file(self, prefix: str, filename: str):
    # TODO: don't upload input, get real output instead
    local_path = os.path.join(".data", prefix, "input", filename)
    remote_path = f"{prefix}/output/{filename}"

    self.s3.upload_file(local_path, "treetime", remote_path)
