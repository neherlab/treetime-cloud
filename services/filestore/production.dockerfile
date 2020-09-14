FROM minio/minio:RELEASE.2020-02-07T23-28-16Z

CMD ["minio", "server", "--quiet", "/tmp"]
