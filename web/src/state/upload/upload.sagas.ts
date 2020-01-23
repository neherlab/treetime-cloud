import { call, takeLatest } from 'redux-saga/effects'

import axios from 'axios'

import {
  triggerUploadFiles,
  uploadFilesAsync,
  UploadFilesPayload,
} from './upload.actions'

import fsaSaga from '../util/fsaSaga'

import { FileType } from './upload.types'

interface UploadFileParams {
  type: FileType
  file: File
}

function uploadFilesApi({ files }: UploadFilesPayload) {
  // TODO: upload in parallel
  return [...files.entries()].map(([type, file]) => {
    const fileType = type.toLowerCase()
    const formData = new FormData()
    formData.set('taskId', 'TODO-123-456') // TODO: get real task id
    formData.append('file', file)
    return axios.post(
      `http://localhost:5000/api/v1/upload/${fileType}`,
      formData,
      {
        headers: {
          'Accept': 'application/json',
          'Content-Type': 'multipart/form-data',
        },
      },
    )
  })
}

export function* uploadFileWorker(params: UploadFilesPayload) {
  yield call(uploadFilesApi, params)
}

export const uploadSaga = takeLatest(
  triggerUploadFiles,
  fsaSaga(uploadFilesAsync, uploadFileWorker),
)

export default [uploadSaga]
