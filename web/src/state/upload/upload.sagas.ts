import { takeLatest } from 'redux-saga/effects'

import axios from 'axios'

import {
  triggerUploadFiles,
  uploadFilesAsync,
  UploadFilesPayload,
} from './upload.actions'

import fsaSaga from '../util/fsaSaga'

function uploadFilesApi({ files, taskId }: UploadFilesPayload) {
  const formData = new FormData()

  formData.set('taskId', taskId)

  files.forEach((file, type) => {
    const typeString = type.toString()
    formData.append(typeString, file)
  })

  return axios.post(`/api/v1/upload`, formData, {
    headers: {
      'Accept': 'application/json',
      'Content-Type': 'multipart/form-data',
    },
  })
}

export const uploadSaga = takeLatest(
  triggerUploadFiles,
  fsaSaga(uploadFilesAsync, uploadFilesApi),
)

export default [uploadSaga]
