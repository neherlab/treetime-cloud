import { BadRequestException, Injectable } from '@nestjs/common'

import S3 from 'aws-sdk/clients/s3'
import { concurrent } from 'fasy'
import { isEqual } from 'lodash'

export enum FileType {
  DATES = 'DATES',
  FASTA = 'FASTA',
  NWK = 'NWK',
}

export interface UploadedFileData {
  [key: string]: Express.Multer.File[]
}

export interface UploadedFilepaths {
  [key: string]: string | undefined
}

export interface Upload {
  type: string
  file: Express.Multer.File
}

// export async function upload(prefix: string) {
//   return async ({ type, file }: Upload) =>
// }

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
    const expectedFileTypes = Object.keys(FileType)

    const fileTypes = Object.keys(filedata).sort()
    if (!isEqual(fileTypes, expectedFileTypes)) {
      const numExpectedFileTypes = expectedFileTypes.length
      const fileTypesList = fileTypes.join(', ')
      const expectedFileTypesList = expectedFileTypes.join(', ')

      throw new BadRequestException(
        `Expected ${numExpectedFileTypes} files to be uploaded of types: "${expectedFileTypesList}", but got "${fileTypesList}"'`,
      )
    }

    const files: Upload[] = Object.entries(filedata).map(([type, files]) => ({ type, file: files[0] }))
    return concurrent.map(({ type, file }: Upload) => {
      if (!file) {
        throw new BadRequestException(`Unable to read file of type ${type}`)
      }

      const filename = file.originalname ?? file.filename
      if (!filename) {
        throw new BadRequestException(`File with type ${type} did not have a filename`)
      }

      const filepath = `${prefix}/input/${filename}`
      const data = file.buffer

      // TODO: these 2 operations should be executed atomically:
      // if one operation fails, all changes should be reverted.
      return Promise.all([this.addFilenameToTask(type, prefix, filename), this.uploadFileToS3(filepath, data)])
    }, files)
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
