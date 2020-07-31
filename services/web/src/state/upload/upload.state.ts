import { FileType } from './upload.types'

export interface UploadState {
  files: Map<FileType, File>
  pending: number
  error: Error | null
}

export const uploadDefaultState: UploadState = {
  files: new Map<FileType, File>(),
  pending: 0, // number of async updates in-flight
  error: null,
}
