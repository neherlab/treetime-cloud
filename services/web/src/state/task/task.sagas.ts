import axios from 'axios'
import { takeLatest } from 'redux-saga/effects'

import fsaSaga from 'src/state/util/fsaSaga'

import { postTaskAsync, getTaskIdAsync, PostTaskPayload, getTaskIdTrigger, postTaskTrigger } from './task.actions'

function apiGetTaskId() {
  return axios.get(`/api/v1/taskId`, {
    headers: {
      'Accept': 'application/json',
      'Content-Type': 'multipart/form-data',
    },
  })
}

export const sagaGetTaskId = takeLatest(getTaskIdTrigger, fsaSaga(getTaskIdAsync, apiGetTaskId))

function apiPostTask(payload: PostTaskPayload) {
  return axios.post(`/api/v1/task`, { payload })
}

export const sagaPostTask = takeLatest(postTaskTrigger, fsaSaga(postTaskAsync, apiPostTask))

export default [sagaGetTaskId, sagaPostTask]
