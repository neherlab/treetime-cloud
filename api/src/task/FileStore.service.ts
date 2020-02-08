import { Injectable } from '@nestjs/common'

import S3 from 'aws-sdk/clients/s3'

@Injectable()
export class FileStoreService {
  // TODO: replace hardcoded values with dynamically configurable variables
  private readonly s3 = new S3({
    accessKeyId: 'minioadmin',
    secretAccessKey: 'minioadmin',
    endpoint: 'http://treetime-dev-file_store:9000',
    s3ForcePathStyle: true,
    signatureVersion: 'v4',
    region: 'eu-central-1',
  })

  // TODO: replace hardcoded values with dynamically configurable variables
  private readonly BUCKET_NAME = 'treetime'

  public async uploadInputFiles(prefix: string, files: Express.Multer.File[]) {
    return Promise.all(files.map(file => this.uploadInputFile(prefix, file)))
  }

  public async uploadInputFile(prefix: string, file: Express.Multer.File) {
    const filename = file?.originalname ?? file?.filename
    const filepath = `${prefix}/${filename}`
    const data = file.buffer
    return this.uploadFileToS3(filepath, data)
  }

  public async uploadFileToS3(filepath: string, data: Buffer) {
    return this.s3
      .putObject({
        Bucket: this.BUCKET_NAME,
        Key: filepath,
        Body: data,
      })
      .promise()
  }
}
