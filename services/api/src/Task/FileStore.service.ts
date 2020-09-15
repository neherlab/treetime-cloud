import { BadRequestException, Injectable } from '@nestjs/common'
import S3 from 'aws-sdk/clients/s3'
import { concurrent } from 'fasy'
import { isEqual } from 'lodash'

import { getenv } from '../../lib/getenv'

const FILESTORE_HOST = getenv('FILESTORE_HOST')
const FILESTORE_PORT = getenv('FILESTORE_PORT')
const FILESTORE_ADDRESS = `${FILESTORE_HOST}:${FILESTORE_PORT}`

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

@Injectable()
export class FileStoreService {
  // TODO: replace hardcoded values with dynamically configurable variables
  private readonly s3 = new S3({
    accessKeyId: 'minioadmin',
    secretAccessKey: 'minioadmin',
    endpoint: FILESTORE_ADDRESS,
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
    await this.ensureBucket()

    return this.s3
      .putObject({
        Bucket: this.BUCKET_NAME,
        Key: filepath,
        Body: data,
      })
      .promise()
  }

  private async ensureBucket() {
    try {
      await this.s3
        .createBucket({
          Bucket: this.BUCKET_NAME,
          CreateBucketConfiguration: {
            LocationConstraint: 'eu-central-1',
          },
        })
        .promise()
    } catch (error) {
      // This is a bug in aws-sdk that prevents from implementing this check in a type-safe way
      // See: https://github.com/aws/aws-sdk-js/issues/2611
      // eslint-disable-next-line @typescript-eslint/no-unsafe-member-access
      if (/* !(error instanceof AWSError) || */ error?.code !== 'BucketAlreadyOwnedByYou') {
        throw error
      }
    }
  }
}
