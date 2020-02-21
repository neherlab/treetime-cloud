import { takeLatest } from 'redux-saga/effects'

import axios from 'axios'

import {
  asyncPostTask,
  getTaskIdAsync,
  PostTaskPayload,
  triggerGetTaskId,
  triggerPostTask,
} from './task.actions'

import fsaSaga from '../util/fsaSaga'

function apiGetTaskId() {
  return axios.get(`/api/v1/taskId`, {
    headers: {
      'Accept': 'application/json',
      'Content-Type': 'multipart/form-data',
    },
  })
}

export const sagaGetTaskId = takeLatest(
  triggerGetTaskId,
  fsaSaga(getTaskIdAsync, apiGetTaskId),
)

function apiPostTask(payload: PostTaskPayload) {
  return axios.post(`/api/v1/task`, { payload })
}

export const sagaPostTask = takeLatest(
  triggerPostTask,
  fsaSaga(asyncPostTask, apiPostTask),
)
