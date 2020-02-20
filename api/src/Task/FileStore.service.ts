import { Injectable } from '@nestjs/common'

import S3 from 'aws-sdk/clients/s3'

export interface UploadedFileData {
  DATES: Express.Multer.File[]
  FASTA: Express.Multer.File[]
  NWK: Express.Multer.File[]
}

export interface UploadedFiles {
  DATES?: Express.Multer.File
  FASTA?: Express.Multer.File
  NWK?: Express.Multer.File
}

export interface UploadedFilepaths {
  DATES?: string
  FASTA?: string
  NWK?: string
}

@Injectable()
export class FileStoreService {
  // TODO: replace hardcoded values with dynamically configurable variables
  private readonly s3 = new S3({
    accessKeyId: 'minioadmin',
    secretAccessKey: 'minioadmin',
    endpoint: 'http://treetime-dev-filestore:9000',
    s3ForcePathStyle: true,
    signatureVersion: 'v4',
    region: 'eu-central-1',
  })

  // TODO: replace hardcoded values with dynamically configurable variables
  private readonly BUCKET_NAME = 'treetime'

  // TODO: persist this in a database instead
  private readonly filenames: Map<string, UploadedFilepaths> = new Map<string, UploadedFilepaths>() // prettier-ignore

  public async uploadInputFiles(prefix: string, filedata: UploadedFileData) {
    const files: UploadedFiles = {
      DATES: filedata?.DATES?.[0],
      FASTA: filedata?.FASTA?.[0],
      NWK: filedata?.NWK?.[0],
    }

    return Promise.all([
      Object.entries(files).map(async ([type, file]) => {
        if (!file) {
          return
        }

        const filename = file?.originalname ?? file?.filename
        const filepath = `${prefix}/input/${filename}`
        const data = file.buffer

        // TODO: these 2 operations should be executed atomically:
        // if one operation fails, all changes should be reverted.
        await this.addFilenameToTask(type, prefix, filename)
        await this.uploadFileToS3(filepath, data)
      }),
    ])
  }

  public async getFilenamesForTask(taskId: string) {
    return this.filenames.get(taskId)
  }

  // prettier-ignore
  private async addFilenameToTask(type: string, prefix: string, filename: string) {
    const existingFiles = this.filenames.get(prefix) ?? {}
    this.filenames.set(prefix, { ...existingFiles, [type]: filename })
  }

  private async uploadFileToS3(filepath: string, data: Buffer) {
    return this.s3
      .putObject({
        Bucket: this.BUCKET_NAME,
        Key: filepath,
        Body: data,
      })
      .promise()
  }
}
